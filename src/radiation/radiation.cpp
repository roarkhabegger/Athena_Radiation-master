//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation
//======================================================================================


#include <sstream>  // msg
#include <iostream>  // cout
#include <stdexcept> // runtime erro
#include <stdio.h>  // fopen and fwrite


// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp" 
#include "radiation.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "integrators/rad_integrators.hpp"

// constructor, initializes data structures and parameters

// The default opacity function.
// Do nothing. Keep the opacity as the initial value
inline void DefaultOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
{
  
}


Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin):
    pmy_block(pmb), ir(pmb->ncells3,pmb->ncells2,pmb->ncells1,pmb->nfre_ang),
    ir1(pmb->ncells3,pmb->ncells2,pmb->ncells1,pmb->nfre_ang),
// constructor overload resolution of non-aggregate class type AthenaArray<Real>
    flux{ {pmb->ncells3, pmb->ncells2, pmb->ncells1+1, pmb->nfre_ang},
      {pmb->ncells3, pmb->ncells2+1, pmb->ncells1, pmb->nfre_ang,
       (pmb->pmy_mesh->f2 ? AthenaArray<Real>::DataStatus::allocated :
        AthenaArray<Real>::DataStatus::empty)},
      {pmb->ncells3+1, pmb->ncells2, pmb->ncells1, pmb->nfre_ang,
       (pmb->pmy_mesh->f3 ? AthenaArray<Real>::DataStatus::allocated :
        AthenaArray<Real>::DataStatus::empty)}},
    coarse_ir_(pmb->ncc3, pmb->ncc2, pmb->ncc1,pmb->nfre_ang,
             (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
              AthenaArray<Real>::DataStatus::empty)),
    rad_bvar(pmb, &ir, &coarse_ir_, flux){
  // read in the parameters
  int nmu = pin->GetInteger("radiation","nmu");
  // total number of polar angles covering 0 to pi/2
  nzeta = pin->GetOrAddInteger("radiation","nzeta",0); 
  // total number of azimuthal angles covering 0 to pi
  npsi = pin->GetOrAddInteger("radiation","npsi",0); 
  angle_flag = pin->GetOrAddInteger("radiation","angle_flag",0);
  prat = pin->GetReal("radiation","Prat");
  crat = pin->GetReal("radiation","Crat");
  rotate_theta=pin->GetOrAddInteger("radiation","rotate_theta",0);
  rotate_phi=pin->GetOrAddInteger("radiation","rotate_phi",0);
  reduced_c  = crat * pin->GetOrAddReal("radiation","reduced_factor",1.0);
  nfreq = pin->GetOrAddInteger("radiation","n_frequency",1);
  vmax = pin->GetOrAddReal("radiation","vmax",0.9);
  tunit = pin->GetOrAddReal("radiation","Tunit",1.e7);
  t_floor_ = pin->GetOrAddReal("radiation", "tfloor", TINY_NUMBER);

  Mesh *pm = pmb->pmy_mesh;

//  ir_output=pin->GetOrAddInteger("radiation","ir_output",0);
  
  set_source_flag = 0;

  // equivalent temperature for electron
  telectron = 5.94065e9;
  telectron /= tunit;

  // number of cells for three dimensions
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;  
  // calculate noct based on dimension
  int ndim = 1;
  if(nc2 > 1) ndim = 2;
  if(nc3 > 1) ndim = 3;
  
   
  int n_ang; // number of angles per octant and number of octant
  // total calculate total number of angles based on dimensions
  if(angle_flag == 1){

    if(ndim == 1){
      noct = 2;
      n_ang = nzeta;
    }else if(ndim == 2){
      if(npsi <= 1){
        n_ang = nzeta;
      }else if(nzeta == 0){
        n_ang = npsi/2;
      }else{
        n_ang = nzeta*npsi;
      }
      noct = 4;
    }else if(ndim == 3){
      n_ang = nzeta*npsi/2;
      noct = 8;
    }

  }else{ 

    if(ndim == 1){
      n_ang = nmu;
      noct = 2;
    }else if(ndim == 2){
      noct = 4;
      if(angle_flag == 0){
        n_ang = nmu * (nmu + 1)/2;
      }else if(angle_flag == 10){
        n_ang = nmu;
      }
    }else if(ndim == 3){
      noct = 8;
      if(angle_flag == 0){
        n_ang = nmu * (nmu + 1)/2;
      }else if(angle_flag == 10){
        n_ang = nmu * nmu/2;
      }
    }// end 3D
  }
  
  nang = n_ang * noct;
  
  n_fre_ang = nang * nfreq;

//  if(ir_output > n_fre_ang){

//    std::stringstream msg;
//    msg << "### FATAL ERROR in Radiation Class" << std::endl
//        << "number of output specific intensity is too large!";
//    throw std::runtime_error(msg.str().c_str());
//  }
  
//  if(ir_output > 0){
//    ir_index.NewAthenaArray(ir_output);
//    dump_ir.NewAthenaArray(ir_output,nc3,nc2,nc1);
//  }
  
  pmb->RegisterMeshBlockData(ir);

  // If user-requested time integrator is type 3S*, allocate additional memory registers
  std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");
  if (integrator == "ssprk5_4" || STS_ENABLED) {
    // future extension may add "int nregister" to Hydro class
    ir2.NewAthenaArray(nc3, nc2, nc1, n_fre_ang);
  }
  
 // Do not add to cell-centered refinement, as 
// radiation variables need to be done in different order      
//  if(pm->multilevel){
//    refinement_idx = pmy_block->pmr->AddToRefinement(&ir, &coarse_ir_);
//  }
  
  rad_mom.NewAthenaArray(13,nc3,nc2,nc1);
  rad_mom_cm.NewAthenaArray(4,nc3,nc2,nc1);
  sigma_s.NewAthenaArray(nc3,nc2,nc1,nfreq);
  sigma_a.NewAthenaArray(nc3,nc2,nc1,nfreq);
  sigma_ae.NewAthenaArray(nc3,nc2,nc1,nfreq);
  sigma_planck.NewAthenaArray(nc3,nc2,nc1,nfreq);
  
  grey_sigma.NewAthenaArray(3,nc3,nc2,nc1);

  
  mu.NewAthenaArray(3,nc3,nc2,nc1,nang);
  wmu.NewAthenaArray(nang);

  wfreq.NewAthenaArray(nfreq);
  
  if(angle_flag == 1)
    AngularGrid(angle_flag, nzeta, npsi);
  else
    AngularGrid(angle_flag, nmu);    

  
  // Initialize the frequency weight
  FrequencyGrid();
  
  // set a default opacity function
  UpdateOpacity = DefaultOpacity;
  
  
  pradintegrator = new RadIntegrator(this, pin);
  
// enroll radiation boundary value object
  rad_bvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&rad_bvar);
  pmb->pbval->bvars_main_int.push_back(&rad_bvar);

  
  // dump the angular grid and radiation parameters in a file
  if(Globals::my_rank ==0){
    FILE *pfile;
    std::stringstream msg;
    if((pfile = fopen("Rad_angles.txt","w")) == NULL){
        msg << "### FATAL ERROR in Radiation Class" << std::endl
            << "Output file Rad_angles.txt could not be opened";
        throw std::runtime_error(msg.str().c_str());
    }
      // damp the angular grid in one cell

    fprintf(pfile,"Prat          %4.2e \n",prat);
    fprintf(pfile,"Crat          %4.2e \n",crat);
    fprintf(pfile,"reduced_c     %4.2e \n",reduced_c);
    fprintf(pfile,"Vmax          %4.2e \n",vmax);
    fprintf(pfile,"Tunit         %4.2e \n",tunit);
    fprintf(pfile,"Compt         %d  \n",pradintegrator->compton_flag_);
    fprintf(pfile,"Planck        %d  \n",pradintegrator->planck_flag_);
    fprintf(pfile,"Tfloor        %4.2e \n",t_floor_);
    fprintf(pfile,"rotate_theta  %d  \n",rotate_theta);
    fprintf(pfile,"rotate_phi    %d  \n",rotate_phi);
    fprintf(pfile,"adv_flag:     %d  \n",pradintegrator->adv_flag_);
    fprintf(pfile,"nzeta:        %d  \n",nzeta);
    fprintf(pfile,"npsi:         %d  \n",npsi);
    
    for(int n=0; n<nang; ++n){
      fprintf(pfile,"%2d   %e   %e   %e    %e\n",n,mu(0,0,0,0,n),mu(1,0,0,0,n),
             mu(2,0,0,0,n), wmu(n));
    }

    
    fclose(pfile);
  
  }
  
  

}

// destructor
// destructor not used
//Radiation::~Radiation()
//{
//  ir.DeleteAthenaArray();
//  ir1.DeleteAthenaArray();
//  rad_mom.DeleteAthenaArray();
//  rad_mom_cm.DeleteAthenaArray();
//  sigma_s.DeleteAthenaArray();
//  sigma_a.DeleteAthenaArray();
//  sigma_ae.DeleteAthenaArray();
//  sigma_planck.DeleteAthenaArray();
//  grey_sigma.DeleteAthenaArray();
  
//  if(ir_output > 0){
//    ir_index.DeleteAthenaArray();
//    dump_ir.DeleteAthenaArray();
//  }
  
///  mu.DeleteAthenaArray();
//  wmu.DeleteAthenaArray();
//  wfreq.DeleteAthenaArray();
  
//  flux[X1DIR].DeleteAthenaArray();
//  if(pmy_block->block_size.nx2 > 1) flux[X2DIR].DeleteAthenaArray();
//  if(pmy_block->block_size.nx3 > 1) flux[X3DIR].DeleteAthenaArray();
  
//  delete pradintegrator;
  
//}


//Enrol the function to update opacity

void Radiation::EnrollOpacityFunction(OpacityFunc MyOpacityFunction)
{
  UpdateOpacity = MyOpacityFunction;
  
}




