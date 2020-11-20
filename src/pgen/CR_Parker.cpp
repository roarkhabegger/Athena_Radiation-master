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

static Real sigma=1.e8;
static Real vx = 0.0;
static Real vy = 0.0;
static Real vz = 0.0;
static int direction =0;
Real Pcr = 0.0;
Real PcrBkg = 0.0;//.000001;
Real crFlux = 0.0;

Real dens0 = 1.0;
Real B0    = 1.0;
Real E0    = 1.0;
Real P0    = 1.0;

Real crRad = 0.0;

Real CRtLim = 0.0;
Real H = 1.0;
Real gH = 0.01;

Real g0 = -1.0;

Real LogMeshSpacingX2(Real x, RegionSize rs);

void myGravity(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons);

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void FixMHDLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh);
void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh);

void FixMHDRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh);
void FixCRsourceRight(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh);


void FixMHDBottom(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh);
void FixCRsourceBottom(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh);


Real LogMeshSpacingX2(Real x, RegionSize rs){
  Real xf, xrat;
  xrat = pow(rs.x2max/rs.x2min,1.0/((Real) rs.nx2));
  xf = rs.x2min*pow(xrat,x*rs.nx2); 
  return xf;
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 
 
    
}

void Mesh::InitUserMeshData(ParameterInput *pin){

  Real x2rat = pin->GetOrAddReal("mesh","x2rat",0.0);
  if (x2rat< 0.0) {
    EnrollUserMeshGenerator(X2DIR,LogMeshSpacingX2);
  }

  if (pin->GetString("mesh","ix2_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x2, FixMHDBottom);
  }
  //if (pin->GetString("mesh","ix1_bc")=="user"){
  //  EnrollUserBoundaryFunction(inner_x1, FixMHDLeft);
  //}
  //if (pin->GetString("mesh","ox1_bc")=="user"){
  //  EnrollUserBoundaryFunction(outer_x1, FixMHDRight);
  //}
  EnrollUserExplicitSourceFunction(myGravity);
  if(CR_ENABLED){  
    EnrollUserCRBoundaryFunction(inner_x1, FixCRsourceLeft);
    //EnrollUserCRBoundaryFunction(outer_x1, FixCRsourceRight);
    //EnrollUserCRBoundaryFunction(inner_x2, FixCRsourceBottom);
  }
}



void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  if(CR_ENABLED){
    pcr->EnrollOpacityFunction(Diffusion);
  }


}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
  Real gamma = peos->GetGamma();

  Real vx=pin->GetReal("problem","xVel");
  Pcr = pin->GetReal("problem","CRPres");
  PcrBkg = pin->GetReal("problem","BackgroundCRPres");

  B0   = pin->GetReal("problem","Bx");
  dens0 = pin->GetReal("problem","Dens");
  H = pin->GetReal("problem","ScaleHeight");
  gH = pin->GetReal("problem","gScaleHeight");
  P0 = pin->GetReal("problem","Pres");
  E0 = P0/(gamma-1)+PcrBkg*3.0;
  g0 = pin->GetReal("problem","grav");

  crFlux = 3.0*B0*Pcr/(sqrt(4.0*PI*dens0)) ;
 
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x2 = pcoord->x2v(j);
         
        Real prof = pow(cosh(x2/gH),-1.0*gH/H);
        Real density = dens0*prof;  
        Real pressure = P0*prof;
        Real vel = vx;
        Real crPress = PcrBkg*prof;
        Real B =  B0*sqrt(prof);

        Real gProf = g0*tanh(x2/(gH));

        phydro->u(IDN,k,j,i) = density;

        phydro->u(IM1,k,j,i) = vx*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){
          phydro->u(IEN,k,j,i) = (0.5*density*vx*vx+pressure/(gamma-1.0));
        }
        
        if(CR_ENABLED){
            pcr->u_cr(CRE,k,j,i) = 3.0*crPress;
            pcr->u_cr(CRF1,k,j,i) = 3.0*B*crPress/sqrt(4.0*PI*density);
            pcr->u_cr(CRF2,k,j,i) = 0.0;
            pcr->u_cr(CRF3,k,j,i) = 0.0;
            phydro->u(IEN,k,j,i) += 3.0*crPress;
        }
      }// end i
    }
  }
  //Need to set opactiy sigma in the ghost zones
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

    // Add horizontal magnetic field lines, to show streaming and diffusion 
  // along magnetic field ines
  if(MAGNETIC_FIELDS_ENABLED){

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          Real x2 = pcoord->x2v(j);
          //Real prof = exp(-x2/H);
          Real prof = pow(cosh(x2/gH),-1.0*gH/H);
          pfield->b.x1f(k,j,i) = B0*sqrt(prof);
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = 0.0;
          }
        }
      }

    }

    if(block_size.nx3 > 1){

      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }// end nx3

    // set cell centerd magnetic field
    // Add magnetic energy density to the total energy
    pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);
    E0 += 0.5*SQR(B0);
    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));
      
        }
      }
    }

  }// end MHD
    
  return;
}

void myGravity(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons) 
{  
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x2 = pmb->pcoord->x2v(j);
        Real gProf = g0*tanh(x2/(gH));
        Real src = dt*prim(IDN,k,j,i)*gProf;
        cons(IM2,k,j,i) += src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);
      }
    }
  }


  return;
}


void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc){ 
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
         Real va = 1.0/sqrt(prim(IDN,k,j,i));
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


void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr, 
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie, 
    int js, int je, int ks, int ke, int ngh){
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {
          Real x2 = abs(pco->x2v(j));
          //std::cout << "(" << i << "," << j << "," << k << ") up=" << upBound << ", lo=" << loBound <<std::endl; 
          if ((time<CRtLim) && (x2 <= crRad)) {
            u_cr(CRE,k,j,is-i) = 3.0*Pcr;
            u_cr(CRF1,k,j,is-i) = crFlux;//u_cr(CRF1,k,j,is);
            u_cr(CRF2,k,j,is-i) = u_cr(CRF2,k,j,is);
            u_cr(CRF3,k,j,is-i) = u_cr(CRF3,k,j,is);
          } else {
            Real prof = exp(-1.0*abs(x2)/H);
            u_cr(CRE,k,j,is-i) = 3.0*PcrBkg*prof;
            u_cr(CRF1,k,j,is-i) = 3.0*B0*PcrBkg*pow(prof,1.5)/sqrt(4.0*PI*dens0);//u_cr(CRF1,k,j,is);
            u_cr(CRF2,k,j,is-i) = 0.0;//u_cr(CRF2,k,j,is);
            u_cr(CRF3,k,j,is-i) = 0.0;//u_cr(CRF3,k,j,is);
          }
        }
      }
    }
  }
}



void FixMHDBottom(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh){
  // copy hydro variables into ghost zones, reflecting v1
  for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
      for (int j=1; j<=ngh; ++j) {
        Real x2 = pco->x2v(js-j);
        Real prof = pow(cosh(x2/gH),-1.0*gH/H);
        //Real prof = exp(-1.0*x2/H);
        //std::cout << "j="<< js-j << " dens0=" << dens0 <<std::endl;
        prim(IDN,k,js-j,i) = dens0;
        prim(IVX,k,js-j,i) = 0.0; // reflect 1-velocity
        prim(IVY,k,js-j,i) = 0.0;
        prim(IVZ,k,js-j,i) = 0.0;
        if(NON_BAROTROPIC_EOS)
          prim(IEN,k,js-j,i) = P0;
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones, reflecting b1
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) { 
    for (int i=is; i<=ie; ++i) { 
#pragma simd
      for (int j=1; j<=(NGHOST); ++j) { 
          b.x2f(k,js-j,i) =  0.0;
      } 
    }}
    if(ie > is){ 
     for (int k=ks; k<=ke; ++k) {
     for (int i=is; i<=ie+1; ++i) {
#pragma simd
      for (int j=1; j<=(NGHOST); ++j) {
        Real x2 = pco->x2v(js-j);
        Real prof = pow(cosh(x2/gH),-1.0*gH/H);
        //Real prof = exp(-1.0*x2/H);
        b.x1f(k,js-j,i) =  B0*sqrt(prof);
      }
     }}  
    }
    if(ke > ks){        
     for (int k=ks; k<=ke+1; ++k) {
      for (int i=is; i<=ie; ++i) {
#pragma simd
       for (int j=1; j<=(NGHOST); ++j) {
         b.x3f(k,js-j,i) = 0.0;
       }
      }}
    }
  }
}


