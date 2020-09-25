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
//! \file rad_integrators.cpp
//  \brief implementation of radiation integrators
//======================================================================================

#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../radiation.hpp"
#include "rad_integrators.hpp"
#include "../../coordinates/coordinates.hpp"


RadIntegrator::RadIntegrator(Radiation *prad, ParameterInput *pin)
{

  pmy_rad = prad;

  MeshBlock *pmb = prad->pmy_block;
  Coordinates *pco=pmb->pcoord;

  int nang=prad->nang;
  int nfreq=prad->nfreq;

  rad_xorder = pin->GetOrAddInteger("time","rad_xorder",2);
  if (rad_xorder == 3) {
    if (NGHOST < 3){ 
      std::stringstream msg;
      msg << "### FATAL ERROR in radiation reconstruction constructor" << std::endl
          << "rad_xorder=" << rad_xorder <<
          " (PPM) reconstruction selected, but nghost=" << NGHOST << std::endl
          << "Reconfigure with --nghost=3  " <<std::endl;
      ATHENA_ERROR(msg);
    }
  }

  
      // factor to separate the diffusion and advection part
  taufact_ = pin->GetOrAddInteger("radiation","taucell",5);
  compton_flag_=pin->GetOrAddInteger("radiation","Compton",0);
  planck_flag_=pin->GetOrAddInteger("radiation","Planck",0);
  adv_flag_=pin->GetOrAddInteger("radiation","Advection",0);
  flux_correct_flag_ = pin->GetOrAddInteger("radiation","CorrectFlux",0);
  tau_limit_ =  pin->GetOrAddReal("radiation","tau_limit",0);



  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2, 
  ncells3 = pmb->ncells3; 

 
  
  x1face_area_.NewAthenaArray(ncells1+1);
  if(ncells2 > 1) {
    x2face_area_.NewAthenaArray(ncells1);
    x2face_area_p1_.NewAthenaArray(ncells1);
  }
  if(ncells3 > 1) {
    x3face_area_.NewAthenaArray(ncells1);
    x3face_area_p1_.NewAthenaArray(ncells1);
  }
  cell_volume_.NewAthenaArray(ncells1);


  cwidth2_.NewAthenaArray(ncells1);
  cwidth3_.NewAthenaArray(ncells1);

  dflx_.NewAthenaArray(ncells1,prad->n_fre_ang);

  // arrays for spatial recontruction 
  il_.NewAthenaArray(ncells1,prad->n_fre_ang);
  ilb_.NewAthenaArray(ncells1,prad->n_fre_ang);

  ir_.NewAthenaArray(ncells1,prad->n_fre_ang);
  
  if(adv_flag_ > 0){
    temp_i1_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
    temp_i2_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  }

  vel_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  velx_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  vely_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  velz_.NewAthenaArray(ncells3,ncells2,ncells1,prad->n_fre_ang);
  
  vncsigma_.NewAthenaArray(nang);
  vncsigma2_.NewAthenaArray(nang);
  wmu_cm_.NewAthenaArray(nang);
  tran_coef_.NewAthenaArray(nang);
  cm_to_lab_.NewAthenaArray(nang);
  ir_cm_.NewAthenaArray(prad->n_fre_ang);

  if(prad->angle_flag == 1){
    int &nzeta = prad->nzeta;
    int &npsi = prad->npsi;
    if(nzeta > 0){
      g_zeta_.NewAthenaArray(2*nzeta+1);
      q_zeta_.NewAthenaArray(2*nzeta+2*NGHOST);
      ql_zeta_.NewAthenaArray(2*nzeta+2*NGHOST);
      qr_zeta_.NewAthenaArray(2*nzeta+2*NGHOST);
      

      if(npsi > 0){
        zeta_flux_.NewAthenaArray(ncells3,ncells2,ncells1,(2*nzeta+1)*2*npsi);
        zeta_area_.NewAthenaArray(2*npsi,2*nzeta+1);
      }
      else{
        zeta_flux_.NewAthenaArray(ncells3,ncells2,ncells1,2*nzeta+1);
        zeta_area_.NewAthenaArray(2*nzeta+1);
      }

      pco->ZetaArea(prad, zeta_area_);
    }

    if(npsi > 0){
      g_psi_.NewAthenaArray(2*npsi+1);
      q_psi_.NewAthenaArray(2*npsi+2*NGHOST);
      ql_psi_.NewAthenaArray(2*npsi+2*NGHOST);
      qr_psi_.NewAthenaArray(2*npsi+2*NGHOST);      


      if(nzeta > 0){
        psi_flux_.NewAthenaArray(ncells3,ncells2,ncells1,2*nzeta*(2*npsi+1));
        psi_area_.NewAthenaArray(2*nzeta,2*npsi+1);
      }
      else{
        psi_flux_.NewAthenaArray(ncells3,ncells2,ncells1,2*npsi+1); 
        psi_area_.NewAthenaArray(2*npsi+1);
      }

      pco->PsiArea(prad, psi_area_); 
    }

    dflx_ang_.NewAthenaArray(nang);
    ang_vol_.NewAthenaArray(nang);
    pco->AngularVol(prad, ang_vol_);
  }

  // calculate the advection velocity at the cell faces

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;
  jl = js, ju=je, kl=ks, ku=ke;

  if(ncells2 > 1)
  {
    if(ncells3 == 1){
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    }else{
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
 
  }

  // calculate velx_
  for (int k=kl; k<=ku; ++k){
    for (int j=jl; j<=ju; ++j){
        // get the velocity at the interface
      for(int i=is-1; i<=ie+1; ++i){
        Real dxl = pco->x1f(i)-pco->x1v(i-1);
        Real dxr = pco->x1v(i) - pco->x1f(i);
        Real factl = dxr/(dxl+dxr);
        Real factr = dxl/(dxl+dxr);
        for(int ifr=0; ifr<nfreq; ++ifr){
          Real *cosx = &(prad->mu(0,k,j,i-1,0));
          Real *cosx1 = &(prad->mu(0,k,j,i,0));
          Real *veln = &(velx_(k,j,i,ifr*nang));
#pragma omp simd aligned(cosx,cosx1,veln:ALI_LEN)
          for(int n=0; n<nang; ++n){
            // linear intepolation between x1v(i-1), x1f(i), x1v(i)
            veln[n] = prad->reduced_c *
                                (factl * cosx[n] + factr * cosx1[n]);                       
          }// end n
        }// end ifr
      }// end i
    }
  }// end k

  // calculate vely_
  if(ncells2 > 1){
    il = is-1, iu = ie+1, kl = ks, ku = ke;
    if (ncells3 >  1) // 2D
      kl = ks-1, ku = ke+1;

    for (int k=kl; k<=ku; ++k){
      for (int j=js; j<=je+1; ++j){
        // get the velocity at the interface
        for(int i=il; i<=iu; ++i){
          Real dxl = pco->x2f(j)-pco->x2v(j-1);
          Real dxr = pco->x2v(j) - pco->x2f(j);
          Real factl = dxr/(dxl+dxr);
          Real factr = dxl/(dxl+dxr);
          for(int ifr=0; ifr<nfreq; ++ifr){
            Real *cosy = &(prad->mu(1,k,j-1,i,0));
            Real *cosy1 = &(prad->mu(1,k,j,i,0));
            Real *veln = &(vely_(k,j,i,ifr*nang));
#pragma omp simd aligned(cosy,cosy1,veln:ALI_LEN)
            for(int n=0; n<nang; ++n){
            // linear intepolation between x2v(j-1), x2f(j), x2v(j)
              veln[n] = prad->reduced_c *
                          (factl * cosy[n] + factr * cosy1[n]);
            }


          }// end ifr
        }// end i
      }
    }// end k 
  }// ncells2

  // calculate vely_
  if(ncells3 > 1){
    il =is-1, iu=ie+1, jl=js-1, ju=je+1;

    for (int k=ks; k<=ke+1; ++k){
      for (int j=jl; j<=ju; ++j){
        // get the velocity at the interface
        for(int i=il; i<=iu; ++i){
          Real dxl = pco->x3f(k) - pco->x3v(k-1);
          Real dxr = pco->x3v(k) - pco->x3f(k);
          Real factl = dxr/(dxl+dxr);
          Real factr = dxl/(dxl+dxr);
          for(int ifr=0; ifr<nfreq; ++ifr){
            Real *cosz = &(prad->mu(2,k-1,j,i,0));
            Real *cosz1 = &(prad->mu(2,k,j,i,0));
            Real *veln = &(velz_(k,j,i,ifr*nang));
#pragma omp simd aligned(cosz,cosz1,veln:ALI_LEN)
            for(int n=0; n<nang; ++n){
            // linear intepolation between x2v(j-1), x2f(j), x2v(j)
              veln[n] = prad->reduced_c *
                          (factl * cosz[n] + factr * cosz1[n]);
            }
          }// end ifr
        }// end i
      }// end j
    }// end k 
  }// ncells3

}
// destructor

RadIntegrator::~RadIntegrator()
{
 
  x1face_area_.DeleteAthenaArray();
  if(pmy_rad->pmy_block->ncells2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if(pmy_rad->pmy_block->ncells3 > 1) {
    x3face_area_.DeleteAthenaArray();
    x3face_area_p1_.DeleteAthenaArray();
  }
  cell_volume_.DeleteAthenaArray();

  cwidth2_.DeleteAthenaArray();
  cwidth3_.DeleteAthenaArray();

  dflx_.DeleteAthenaArray();

  il_.DeleteAthenaArray();
  ilb_.DeleteAthenaArray();

  ir_.DeleteAthenaArray();

  if(adv_flag_ > 0){
    temp_i1_.DeleteAthenaArray();
    temp_i2_.DeleteAthenaArray();
  }

  vel_.DeleteAthenaArray();
  velx_.DeleteAthenaArray();
  vely_.DeleteAthenaArray();
  velz_.DeleteAthenaArray();

  
  vncsigma_.DeleteAthenaArray();
  vncsigma2_.DeleteAthenaArray();
  wmu_cm_.DeleteAthenaArray();
  tran_coef_.DeleteAthenaArray();
  cm_to_lab_.DeleteAthenaArray();
  ir_cm_.DeleteAthenaArray();

  if(pmy_rad->angle_flag == 1){
    int &nzeta = pmy_rad->nzeta;
    int &npsi = pmy_rad->npsi;
    if(nzeta > 0){
      g_zeta_.DeleteAthenaArray();
      q_zeta_.DeleteAthenaArray();
      ql_zeta_.DeleteAthenaArray();
      qr_zeta_.DeleteAthenaArray();
      zeta_flux_.DeleteAthenaArray();
      zeta_area_.DeleteAthenaArray();
    }
    if(npsi > 0){
      g_psi_.DeleteAthenaArray();
      q_psi_.DeleteAthenaArray();
      ql_psi_.DeleteAthenaArray();
      qr_psi_.DeleteAthenaArray();
      psi_flux_.DeleteAthenaArray();     
      psi_area_.DeleteAthenaArray(); 
    }
    dflx_ang_.DeleteAthenaArray();
    ang_vol_.DeleteAthenaArray();
  }  
  
}

