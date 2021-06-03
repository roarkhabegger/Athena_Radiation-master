//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C++ headers
#define _USE_MATH_DEFINES
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"


//======================================================================================
/*! \file beam.cpp
 *  Dynamic diffusion test for cosmic rays
 * Compare the numerical solution with analytic solution
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================

Real dens0, pres0;
Real alpha;
Real beta;
Real angle;
Real ampSin;
Real ampLin;
Real lambSin; 
Real lambLin;

Real sigma;

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);

void ProjectCROuterX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh);

void ProjectCRInnerX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh);

Real pres(Real x1, Real x2);
Real GradP(Real x1, Real x2);
Real dens(Real x1, Real x2);


//======================================================================================
// HYDRO PROFILES
//======================================================================================
Real pres(Real x1, Real x2) {
  Real myVal = 0.0;
  Real r = pow(pow(x2,2)+pow(x1,2),0.5);         
  Real theta = atan2(x2,x1);
  Real dist = r * cos(0.5*M_PI - angle - theta);
  myVal = pres0*(1+ampLin*x2/lambLin + ampSin*sin(x2/lambSin));

  return myVal;
 
}
Real GradP(Real x1, Real x2) {
  Real myVal = 0.0;
  Real r = pow(pow(x2,2)+pow(x1,2),0.5);         
  Real theta = atan2(x2,x1);
  Real dist = r * cos(0.5*M_PI - angle - theta);
  myVal = pres0*(ampLin/lambLin + ampSin/lambSin*cos(dist/lambSin));

  return myVal;
 
}
Real dens(Real x1, Real x2) {
  Real myVal = 0.0;
  Real r = pow(pow(x2,2)+pow(x1,2),0.5);         
  Real theta = atan2(x2,x1);
  Real dist = r * cos(0.5*M_PI - angle - theta);
  myVal = dens0;//dens0*(1+ampLin*dist/lambLin + ampSin*sin(dist/lambSin));

  return myVal;
 
}



//======================================================================================
// INITUSER_____DATA
//======================================================================================
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  if(CR_ENABLED){
    pcr->EnrollOpacityFunction(Diffusion);
  }


}
void Mesh::InitUserMeshData(ParameterInput *pin){

  if (pin->GetString("mesh","ix2_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x2, ProjectPressureInnerX2);
  }
  if (pin->GetString("mesh","ox2_bc")=="user"){
    EnrollUserBoundaryFunction(outer_x2, ProjectPressureOuterX2);
  }
  
  if(CR_ENABLED){
    EnrollUserCRBoundaryFunction(inner_x2, ProjectCRInnerX2);
    EnrollUserCRBoundaryFunction(outer_x2, ProjectCROuterX2);
  }
  

}

//======================================================================================
// PROB GENERATOR
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
  Real gamma = peos->GetGamma();

  dens0 = pin->GetReal("problem","Dens");
  pres0 = pin->GetReal("problem","Pres");
  ampSin = pin->GetReal("problem","sinAmp");
  ampLin = pin->GetReal("problem","linAmp");
  lambSin = pin->GetReal("problem","sinLambda");
  lambLin = pin->GetReal("problem","linLambda");


  if (MAGNETIC_FIELDS_ENABLED) {
    alpha = pin->GetOrAddReal("problem","alpha",0.0);
    angle = pin->GetReal("problem","angle");
  }
  if(CR_ENABLED){
    beta = pin->GetOrAddReal("problem","beta",0.0);
    sigma = pin->GetReal("cr","sigma");
  }
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x2 = pcoord->x2v(j);
        Real x1 = pcoord->x1v(i);
        phydro->u(IDN,k,j,i) = dens(x1,x2);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pres(x1,x2)/(gamma-1);
        if(CR_ENABLED){
          Real r = pow(pow(x2,2)+pow(x1,2),0.5);         
          Real theta = atan2(x2,x1);
          Real dist = r * cos(0.5*M_PI - angle - theta);
          Real crp = beta*pres(x1,x2);
          Real B0  = sqrt(alpha*pres0);
          Real vs = B0/pow(4*PI*dens0,0.5)*cos(angle)*(-1)
                    *(ampLin/lambLin+ampSin/lambSin*cos(dist/lambSin))
                    /abs(ampLin/lambLin+ampSin/lambSin*cos(dist/lambSin));

          
          pcr->u_cr(CRE,k,j,i) = 3.0*crp;
          pcr->u_cr(CRF1,k,j,i) = 0.0; 
          pcr->u_cr(CRF2,k,j,i) = 0.0; //4.0*crp*vs;
          pcr->u_cr(CRF3,k,j,i) = 0.0;
        }
      }// end i
    }
  }    
    // initialize uniform interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0   = pow(alpha*pres0,0.5);
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          Real x2 = pcoord->x2v(j);

          pfield->b.x1f(k,j,i) = 0.0;
          
          
        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          Real x2 = pcoord->x2v(j);
          Real x1 = pcoord->x1v(i);
          pfield->b.x2f(k,j,i) = B0;
        }
      }
    }
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x3f(k,j,i) = 0.0;
        }
      }
    }
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            phydro->u(IEN,k,j,i) += 0.5*( pow(pfield->b.x1f(k,j,i),2.0)
                                         +pow(pfield->b.x2f(k,j,i),2.0)
                                         +pow(pfield->b.x3f(k,j,i),2.0)
                                         );
          }
        }
      }
    }
  }
  if(CR_ENABLED){

  // Default values are 1/3
    int nz1 = block_size.nx1 + 2*(NGHOST);
    int nz2 = block_size.nx2;
    if(nz2 > 1) nz2 += 2*(NGHOST);
    int nz3 = block_size.nx3;
    if(nz3 > 1) nz3 += 2*(NGHOST);
    for(int k=0; k<nz3; ++k){
      for(int j=0; j<nz2; ++j){
        for(int i=0; i<nz1; ++i){
          pcr->sigma_diff(0,k,j,i) = sigma;
          pcr->sigma_diff(1,k,j,i) = sigma;
          pcr->sigma_diff(2,k,j,i) = sigma;
        }
      }
    }// end k,j,i

  }// End CR
  return;
}

//======================================================================================
// CR Boundary Conditions
//======================================================================================
void ProjectCRInnerX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh)
{
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
	  Real dx1 = (pco->x2v(js-j) - pco->x2v(js-j+1));
          Real dx2 = (pco->x2v(js-j+1) - pco->x2v(js-j+2));
          u_cr(CRE,k,js-j,i) = dx1/dx2*( u_cr(CRE,k,js-j+1,i) - u_cr(CRE,k,js-j+2,i))
			       + u_cr(CRE,k,js-j+1,i);
          u_cr(CRF1,k,js-j,i) = 0.0; //u_cr(CRF1,k,j,is);
          u_cr(CRF2,k,js-j,i) = dx1/dx2*( u_cr(CRF2,k,js-j+1,i) - u_cr(CRF2,k,js-j+2,i))
			       + u_cr(CRF2,k,js-j+1,i);
          //u_cr(CRF2,k,js-j,i) = 0.0; //u_cr(CRF2,k,j,is);
          u_cr(CRF3,k,js-j,i) = 0.0; //u_cr(CRF3,k,j,is);

        }
      }
    }
  }
}
void ProjectCROuterX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh)
{
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
	  Real dx1 = (pco->x2v(je+j) - pco->x2v(je+j-1));
          Real dx2 = (pco->x2v(je+j-1) - pco->x2v(je+j-2));
          u_cr(CRE,k,je+j,i) = dx1/dx2*( u_cr(CRE,k,je+j-1,i) - u_cr(CRE,k,je+j-2,i))
			       + u_cr(CRE,k,je+j-1,i);
          u_cr(CRF1,k,je+j,i) = 0.0; //u_cr(CRF1,k,j,is);
          u_cr(CRF2,k,je+j,i) = dx1/dx2*( u_cr(CRF2,k,je+j-1,i) - u_cr(CRF2,k,je+j-2,i))
			       + u_cr(CRF2,k,je+j-1,i);
          //u_cr(CRF2,k,je+j,i) = 0.0; //u_cr(CRF2,k,j,is);
          u_cr(CRF3,k,je+j,i) = 0.0; //u_cr(CRF3,k,j,is);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
	  Real dx1 = (pco->x2v(jl-j) - pco->x2v(jl-j+1));
          Real dx2 = (pco->x2v(jl-j+1) - pco->x2v(jl-j+2));
          prim(IDN,k,jl-j,i) = dx1/dx2*( prim(IDN,k,jl-j+1,i) - prim(IDN,k,jl-j+2,i))
			       + prim(IDN,k,jl-j+1,i);
          prim(IVX,k,jl-j,i) = 0.0;
          prim(IVY,k,jl-j,i) =  prim(IVY,k,jl+j-1,i);  // reflect 2-velocity
          prim(IVZ,k,jl-j,i) = 0.0;
          prim(IPR,k,jl-j,i) = dx1/dx2*( prim(IPR,k,jl-j+1,i) - prim(IPR,k,jl-j+2,i))
			       + prim(IPR,k,jl-j+1,i);
        }
      }
    }
  //}

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0   = pow(alpha*pres0,0.5);
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(jl-j),i) =  0.0;
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(jl-j),i) = 1.0* b.x2f(k,jl+j-1,i);// std::pow(alpha*pressure,0.5);  // reflect 2-field
        }
      }
    }
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(jl-j),i) = 0.0;// b.x3f(k,(jl+j-1),i);
        }
      }
    }
  }
}
//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureOuterX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
	  Real dx1 = (pco->x2v(ju+j) - pco->x2v(ju+j-1));
          Real dx2 = (pco->x2v(ju+j-1) - pco->x2v(ju+j-2));
          prim(IDN,k,ju+j,i) = dx1/dx2*( prim(IDN,k,ju+j-1,i) - prim(IDN,k,ju+j-2,i))
	           	       + prim(IDN,k,ju+j-1,i);
          prim(IVX,k,ju+j,i) = 0.0;
          prim(IVY,k,ju+j,i) = prim(IVY,k,ju-j+1,i);  // reflect 2-velocity
          prim(IVZ,k,ju+j,i) = 0.0;
          prim(IPR,k,ju+j,i) = dx1/dx2*( prim(IPR,k,ju+j-1,i) - prim(IPR,k,ju+j-2,i))
	           	       + prim(IPR,k,ju+j-1,i);
        }
      }
    }
  //}

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0   = pow(alpha*pres0,0.5);
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(ju+j),i) =  0.0;
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(ju+j),i) = b.x2f(k,ju-j+1) ;// std::pow(alpha*pressure,0.5);  // reflect 2-field
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(ju+j),i) =  0.0;
        }
      }
    }

  }
}


//======================================================================================
void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc)
{ 


  // set the default opacity to be a large value in the default hydro case
  CosmicRay *pcr=pmb->pcr;
  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is-1, iu=pmb->ie+1;
  if(pmb->block_size.nx2 > 1){
    jl -= 1;
    ju += 1;
  }
  if(pmb->block_size.nx3 > 1){
    kl -= 1;
    ku += 1;
  }

  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){

        pcr->sigma_diff(0,k,j,i) = sigma;
        pcr->sigma_diff(1,k,j,i) = sigma;
        pcr->sigma_diff(2,k,j,i) = sigma;  

      }
    }
  }

  Real invlim=1.0/pcr->vmax;

  // The information stored in the array
  // b_angle is
  // b_angle[0]=sin_theta_b
  // b_angle[1]=cos_theta_b
  // b_angle[2]=sin_phi_b
  // b_angle[3]=cos_phi_b




  if(MAGNETIC_FIELDS_ENABLED){
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
    // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                         + pcr->cwidth(i);
          Real dprdx=(u_cr(CRE,k,j,i+1) - u_cr(CRE,k,j,i-1))/3.0;
          dprdx /= distance;
          pcr->b_grad_pc(k,j,i) = bcc(IB1,k,j,i) * dprdx;
        }
    //y component
        pmb->pcoord->CenterWidth2(k,j-1,il,iu,pcr->cwidth1);       
        pmb->pcoord->CenterWidth2(k,j,il,iu,pcr->cwidth);
        pmb->pcoord->CenterWidth2(k,j+1,il,iu,pcr->cwidth2);

        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                         + pcr->cwidth(i);
          Real dprdy=(u_cr(CRE,k,j+1,i) - u_cr(CRE,k,j-1,i))/3.0;
          dprdy /= distance;
          pcr->b_grad_pc(k,j,i) += bcc(IB2,k,j,i) * dprdy;

        }
    // z component
        pmb->pcoord->CenterWidth3(k-1,j,il,iu,pcr->cwidth1);       
        pmb->pcoord->CenterWidth3(k,j,il,iu,pcr->cwidth);
        pmb->pcoord->CenterWidth3(k+1,j,il,iu,pcr->cwidth2);

        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                          + pcr->cwidth(i);
          Real dprdz=(u_cr(CRE,k+1,j,i) - u_cr(CRE,k-1,j,i))/3.0;
          dprdz /= distance;
          pcr->b_grad_pc(k,j,i) += bcc(IB3,k,j,i) * dprdz;

          // now only get the sign
//          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) = 1.0;
//          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) 
//            = -1.0;
//          else pcr->b_grad_pc(k,j,i) = 0.0;
        }


      // now calculate the streaming velocity
      // streaming velocity is calculated with respect to the current coordinate 
      //  system
      // diffusion coefficient is calculated with respect to B direction
        for(int i=il; i<=iu; ++i){
          Real pb= bcc(IB1,k,j,i)*bcc(IB1,k,j,i)
                  +bcc(IB2,k,j,i)*bcc(IB2,k,j,i)
                  +bcc(IB3,k,j,i)*bcc(IB3,k,j,i);
          Real inv_sqrt_rho = 1.0/sqrt(prim(IDN,k,j,i));
          Real va1 = bcc(IB1,k,j,i)*inv_sqrt_rho;
          Real va2 = bcc(IB2,k,j,i)*inv_sqrt_rho;
          Real va3 = bcc(IB3,k,j,i)*inv_sqrt_rho;

          Real va = sqrt(pb/prim(IDN,k,j,i));

          Real dpc_sign = 0.0;
          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = 1.0;
          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = -1.0;
          
          pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
          pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
          pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;

          // now the diffusion coefficient

          if(va < TINY_NUMBER){
            pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
          }else{
            pcr->sigma_adv(0,k,j,i) = fabs(pcr->b_grad_pc(k,j,i))
                          /(sqrt(pb)* va * (1.0 + 1.0/3.0) 
                                    * invlim * u_cr(CRE,k,j,i)); 
          }

          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;  

          // Now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(pb);
          if(btot > TINY_NUMBER){
            pcr->b_angle(0,k,j,i) = bxby/btot;
            pcr->b_angle(1,k,j,i) = bcc(IB3,k,j,i)/btot;
          }else{
            pcr->b_angle(0,k,j,i) = 1.0;
            pcr->b_angle(1,k,j,i) = 0.0;
          }
          if(bxby > TINY_NUMBER){
            pcr->b_angle(2,k,j,i) = bcc(IB2,k,j,i)/bxby;
            pcr->b_angle(3,k,j,i) = bcc(IB1,k,j,i)/bxby;
          }else{
            pcr->b_angle(2,k,j,i) = 0.0;
            pcr->b_angle(3,k,j,i) = 1.0;            
          }

        }//        

      }// end j
    }// end k

  }// End MHD  
  else{
  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
  // x component
      pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
      for(int i=il; i<=iu; ++i){
         Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                        + pcr->cwidth(i);
         Real grad_pr=(u_cr(CRE,k,j,i+1) - u_cr(CRE,k,j,i-1))/3.0;
         grad_pr /= distance;

         Real va = 0.0;

         if(va < TINY_NUMBER){
           pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
           pcr->v_adv(0,k,j,i) = 0.0;
         }else{
           Real sigma2 = fabs(grad_pr)/(va * (1.0 + 1.0/3.0) 
                             * invlim * u_cr(CRE,k,j,i)); 
           if(fabs(grad_pr) < TINY_NUMBER){
             pcr->sigma_adv(0,k,j,i) = 0.0;
             pcr->v_adv(0,k,j,i) = 0.0;
           }else{
             pcr->sigma_adv(0,k,j,i) = sigma2;
             pcr->v_adv(0,k,j,i) = -va * grad_pr/fabs(grad_pr);     
           }
        }

        pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
        pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
       
        pcr->v_adv(1,k,j,i) = 0.0;
        pcr->v_adv(2,k,j,i) = 0.0;




      }

    }
  }

  }
}
