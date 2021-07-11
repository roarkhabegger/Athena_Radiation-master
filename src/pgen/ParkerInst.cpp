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
#include <cfloat>      // FLT_MAX
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
//
//Hydrostatic Equilibrium variables
Real dens0, pres0, vx;
Real H, nGrav;
Real g0;
Real alpha;
Real beta;
int cooling;

//Perturbation variables
Real crPertCenterX;
Real crPertCenterY;
Real crPertRadius;
Real crPertAmp;
Real crPertStartTime;
Real crPertEndTime;
Real crPertSteep;

//Profile functions
Real densProfile(Real x1, Real x2);
Real presProfile(Real x1, Real x2);
Real gravProfile(Real x1, Real x2);
Real CRSourceProfile(Real x1, Real x2);

//For logarithmic height spacing
Real LogMeshSpacingX2(Real x, RegionSize rs);

//Gravity function tanh
void mySource(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons);

//CR Source term function
void myCRSource(MeshBlock *pmb, const Real time, const Real dt,
                const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
                AthenaArray<Real> &u_cr);


//Boundary condtions
//Hydro conditions
//X1 boundaries to match HSE
void ProfilesOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void ProfilesInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
//x2 boundaries to extend HSE
void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
//CR boundary conditions
//X1 match HSE on x1 bounds
void ProfilesCROuterX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

void ProfilesCRInnerX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);
//Extend HSE at x2 bounds
void ProjectCROuterX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

void ProjectCRInnerX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

//cr Diffusion
Real sigma;
void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);


Real densProfile(Real x1, Real x2){
  Real rho = dens0;
  Real prof = pow(cosh(x2/(nGrav*H)),-1.0*nGrav);
  rho *= prof;

  return rho;
}
Real presProfile(Real x1, Real x2){
  Real pres = pres0;
  Real rho = densProfile(x1,x2);
  pres *= rho/dens0;

  return pres;
}

Real gravProfile(Real x1, Real x2){
  Real g = g0 * tanh(x2/(nGrav*H));
  return g;
}

Real LogMeshSpacingX2(Real x, RegionSize rs){
  Real xf, xrat;
  xrat = pow(rs.x2max/rs.x2min,1.0/((Real) rs.nx2));
  xf = rs.x2min*pow(xrat,x*rs.nx2);
  return xf;
}

Real CRSourceProfile(Real x1, Real x2){
  Real crE = 0.0;

  return crE;
}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{

  if(CR_ENABLED){
   pcr->EnrollOpacityFunction(Diffusion);
   // pcr->EnrollUserCRSource(myCRSource);
  }


}
void Mesh::InitUserMeshData(ParameterInput *pin){

//  Real x2rat = pin->GetOrAddReal("mesh","x2rat",0.0);
//  if (x2rat< 0.0) {
//    EnrollUserMeshGenerator(X2DIR,LogMeshSpacingX2);
//  }
  if (pin->GetString("mesh","ix2_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x2, ProjectPressureInnerX2);
  }
  if (pin->GetString("mesh","ox2_bc")=="user"){
    EnrollUserBoundaryFunction(outer_x2, ProjectPressureOuterX2);
  }
  if (pin->GetString("mesh","ix1_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x1, ProfilesInnerX1);
  }
  if (pin->GetString("mesh","ox1_bc")=="user"){
    EnrollUserBoundaryFunction(outer_x1, ProfilesOuterX1);
  }

  EnrollUserExplicitSourceFunction(mySource);
  cooling = pin->GetOrAddInteger("problem","cooling",0);

  if(CR_ENABLED){
    if (pin->GetString("mesh","ix2_bc")=="user")
      EnrollUserCRBoundaryFunction(inner_x2, ProjectCRInnerX2);
    if (pin->GetString("mesh","ox2_bc")=="user")
      EnrollUserCRBoundaryFunction(outer_x2, ProjectCROuterX2);
    if (pin->GetString("mesh","ix1_bc")=="user")
      EnrollUserCRBoundaryFunction(inner_x1, ProfilesCRInnerX1);
    if (pin->GetString("mesh","ox1_bc")=="user")
      EnrollUserCRBoundaryFunction(outer_x1, ProfilesCROuterX1);
  //  EnrollUserTimeStepFunction(CRSourceTimeStep);
  }
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Real gamma = peos->GetGamma();

  vx=pin->GetReal("problem","xVel");
  //H = pin->GetReal("problem","ScaleHeight");
  nGrav = pin->GetReal("problem","GravNumScaleHeight");
  beta = pin->GetOrAddReal("problem","beta",0.0);
  alpha = pin->GetOrAddReal("problem","alpha",0.0);

  dens0 = pin->GetReal("problem","Dens");
  g0 = pin->GetReal("problem","grav");
  pres0 = pin->GetReal("problem","Pres");
  H = pres0/dens0/(-1.0*g0)*(1+alpha+beta);
  Real PPertAmp = 0.0; //pin->GetReal("problem","pertPresAmplitude");

  if(CR_ENABLED){
    crPertCenterX = pin->GetReal("problem","pertX");
    crPertCenterY = pin->GetReal("problem","pertY");
    //crPertRadius = pin->GetReal("problem","pertR");
    crPertAmp = pin->GetReal("problem","pertAmplitude");
    //crPertStartTime = pin->GetReal("problem","pertStart");
    //crPertEndTime = pin->GetReal("problem","pertEnd");
    crPertSteep = pin->GetReal("problem","pertDx");
    PPertAmp = pin->GetOrAddReal("problem","pertPresAmplitude",0.0);
  }
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {

      Real x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real density = densProfile(x1,x2);
        Real pressure = presProfile(x1,x2);

        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = vx*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pressure/(gamma-1) + 0.5*density*pow(vx,2.0);
        if(CR_ENABLED){
          Real crp = beta*pressure; //*(1+ampCR*(1-x2/centerCR));
          sigma = pin->GetReal("cr","sigma");
          Real dist = pow(pow(x1-crPertCenterX,2.0)+pow(x2-crPertCenterY,2.0),0.5);
          Real myVal = pres0*crPertAmp/(crPertSteep*pow(2*M_PI,0.5))*exp(-0.5*pow(dist/crPertSteep,2.0));
                       //3.0*pres0*beta*crPertAmp*0.5*(1 - tanh((dist-crPertRadius)/crPertSteep));

          pcr->u_cr(CRE,k,j,i) = 3.0*crp+3.0*beta*myVal;// exp(-40.0*(dist_sq));
          pcr->u_cr(CRF1,k,j,i) = myVal*beta*(2.0*(x1-crPertCenterX)/pow(crPertSteep,2.0))/sigma + vx*4.0*crp;
          pcr->u_cr(CRF2,k,j,i) = 0.0;//sin(atan2(x2,x1))*pcr->u_cr(CRE,k,j,i)*4.0/3.0;
          pcr->u_cr(CRF3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) += PPertAmp*myVal/(gamma-1) ;
        }
      }// end i
    }
  }
    // initialize uniform interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0   = pow(2*alpha*pres0,0.5);
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          Real x2 = pcoord->x2v(j);
          Real x1 = pcoord->x1f(i);
          Real pressure= presProfile(x1,x2);

          Real B = B0*sqrt(pressure/pres0);

          pfield->b.x1f(k,j,i) = B;


        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          Real x2 = pcoord->x2f(j);
          Real x1 = pcoord->x1v(i);
          Real B = B0;// 0.1*B0;
          pfield->b.x2f(k,j,i) = 0.0;
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
            Real x2 = pcoord->x2v(j);
            Real x1 = pcoord->x1v(i);
            Real pressure= presProfile(x1,x2);

            Real B = B0*sqrt(pressure/pres0);

            phydro->u(IEN,k,j,i) += 0.5*( pow(B,2.0));
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

void mySource(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &cons)
{
  Real T0 = 9.773;
  Real T1 = 0.1238;
  Real T2 = 0.007594;
  Real A = 651790.0;
  Real B = 0.705;
  Real C = 3.0;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x2 = pmb->pcoord->x2v(j);
        Real x1 = pmb->pcoord->x1v(i);
        Real gravity = gravProfile(x1,x2);
        Real src = dt*prim(IDN,k,j,i)*gravity;
        cons(IM2,k,j,i) += src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);

        if (cooling == 1) {
          Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i);
          Real n = prim(IDN,k,j,i);
          cons(IEN,k,j,i) -= pow(n,2.0)*(A*exp(-1*T0/(temp-T1))+B*exp(-1*T2/temp))-n*C;
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn ProfilesInnerX1()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm
void ProfilesInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  Real edge = pco->x1v(is);
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
          Real x2 = pco->x2v(j);
          Real x1 = pco->x1v(is-i);
          Real density = densProfile(x1,x2);
          Real pressure = presProfile(x1,x2);
          prim(IDN,k,j,is-i) = density+(prim(IDN,k,j,is) - densProfile(edge,x2));
          prim(IVX,k,j,is-i) = vx+prim(IVX,k,j,is);
          prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
          prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
          prim(IEN,k,j,is-i) = pressure+(prim(IEN,k,j,is) - presProfile(edge,x2));

      }
    }
  }
   // set magnetic field in inlet ghost zones
  edge = pco->x1f(is);
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0   = pow(2*alpha*pres0,0.5);
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
#pragma omp simd
      for(int i=1; i<=ngh; ++i){
          Real x2 = pco->x2v(j);
          Real x1 = pco->x1f(is-i);
          Real pressure= presProfile(x1,x2);

          Real B = B0*sqrt(pressure/pres0);
        b.x1f(k,j,is-i) = B + (b.x1f(k,j,is) - B0*sqrt(presProfile(edge,x2)/pres0));
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
#pragma omp simd
      for(int i=1; i<=ngh; ++i){
        b.x2f(k,j,is-i) = b.x2f(k,j,is);
      }
    }}

    for(int k=ks; k<=ke+1; ++k){
    for(int j=js; j<=je; ++j){
#pragma omp simd
      for(int i=1; i<=ngh; ++i){
        b.x3f(k,j,is-i) = b.x3f(k,j,is);
      }
    }}

  }


  return;
}
//----------------------------------------------------------------------------------------
//! \fn ProfilesOuterX1()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm
void ProfilesOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  Real edge = pco->x1v(ie);
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
          Real x2 = pco->x2v(j);
          Real x1 = pco->x1v(ie+i);
          Real density = densProfile(x1,x2);
          Real pressure = presProfile(x1,x2);
          prim(IDN,k,j,ie+i) = density+(prim(IDN,k,j,ie) - densProfile(edge,x2));
          prim(IVX,k,j,ie+i) = vx+prim(IVX,k,j,ie);
          prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
          prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
          prim(IEN,k,j,ie+i) = pressure+(prim(IEN,k,j,ie) - presProfile(edge,x2));

      }
    }
  }
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    Real fedge = pco->x1f(ie+1);
    Real B0   = pow(2*alpha*pres0,0.5);
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
#pragma omp simd
      for(int i=2; i<=ngh+1; ++i){
          Real x2 = pco->x2v(j);
          Real x1 = pco->x1f(ie+i);
          Real pressure= presProfile(x1,x2);

          Real B = B0*sqrt(pressure/pres0);
        b.x1f(k,j,ie+i) = B+(b.x1f(k,j,ie+1)-B0*sqrt(presProfile(fedge,x2)/pres0));
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
#pragma omp simd
      for(int i=1; i<=ngh; ++i){
        b.x2f(k,j,ie+i) = b.x2f(k,j,ie);
      }
    }}

    for(int k=ks; k<=ke+1; ++k){
    for(int j=js; j<=je; ++j){
#pragma omp simd
      for(int i=1; i<=ngh; ++i){
        b.x3f(k,j,ie+i) = b.x3f(k,j,ie);
      }
    }}

  }


  return;
}
//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {

    Real edge = pco->x2v(jl);
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real x1 = pco->x1v(i);
          Real x2 = pco->x2v(jl-j);

          prim(IDN,k,jl-j,i) = densProfile(x1,x2)+(prim(IDN,k,jl,i) - densProfile(x1,edge));
			       //dx1/dx2*( prim(IDN,k,jl-j+1,i) - prim(IDN,k,jl-j+2,i))
			       //+ prim(IDN,k,jl-j+1,i);
          prim(IVX,k,jl-j,i) = vx+prim(IVX,k,jl,i);
          prim(IVY,k,jl-j,i) = prim(IVY,k,jl,i);  // reflect 2-velocity
          prim(IVZ,k,jl-j,i) = prim(IVZ,k,jl,i);
          prim(IPR,k,jl-j,i) = presProfile(x1,x2)+(prim(IEN,k,jl,i) - presProfile(x1,edge));
			       //dx1/dx2*( prim(IPR,k,jl-j+1,i) - prim(IPR,k,jl-j+2,i))
			       //+ prim(IPR,k,jl-j+1,i);
        }
      }
    }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real fedge = pco->x2v(jl);
    Real B0   = pow(2*alpha*pres0,0.5);
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          Real x2 = pco->x2v(jl-j);
          Real x1 = pco->x1f(i);
          Real pressure= presProfile(x1,x2);

          Real B = B0*sqrt(pressure/pres0);
          b.x1f(k,(jl-j),i) =  B + (b.x1f(k,jl,i) - B0*sqrt(presProfile(x1,fedge)/pres0));
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(jl-j),i) = b.x2f(k,jl,i);// std::pow(alpha*pressure,0.5);  // reflect 2-field
        }
      }
    }
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(jl-j),i) = b.x3f(k,jl,i);// b.x3f(k,(jl+j-1),i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureOuterX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
    Real edge = pco->x2v(ju);
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real x1 = pco->x1v(i);
          Real x2 = pco->x2v(ju+j);

          prim(IDN,k,ju+j,i) = densProfile(x1,x2)+(prim(IDN,k,ju,i) - densProfile(x1,edge));
          prim(IVX,k,ju+j,i) = vx+prim(IVX,k,ju,i);
          prim(IVY,k,ju+j,i) = prim(IVY,k,ju,i);  // reflect 2-velocity
          prim(IVZ,k,ju+j,i) = prim(IVZ,k,ju,i);
          prim(IPR,k,ju+j,i) = presProfile(x1,x2)+(prim(IPR,k,ju,i) - presProfile(x1,edge));
        }
      }
    }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real fedge = pco->x2v(ju);
    Real B0   = pow(2*alpha*pres0,0.5);
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          Real x2 = pco->x2v(ju+j);
          Real x1 = pco->x1f(i);
          Real pressure= presProfile(x1,x2);

          Real B = B0*sqrt(pressure/pres0);
          b.x1f(k,(ju+j),i) =  B + (b.x1f(k,ju,i)-B0*sqrt(presProfile(x1,fedge)/pres0));
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=2; j<=ngh+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(ju+j),i) = b.x2f(k,ju+1,i) ;// std::pow(alpha*pressure,0.5);  // reflect 2-field
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(ju+j),i) =  b.x3f(k,ju,i);
        }
      }
    }

  }
  return;
}
//======================================================================================
// CR Boundary Conditions
//======================================================================================
void ProfilesCRInnerX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {
          Real x1 = pco->x1v(is-i);
          Real x2 = pco->x2v(j);
          u_cr(CRE,k,j,is-i) = 3.0*beta*(presProfile(x1,x2) - presProfile(pco->x1v(is),x2)) + u_cr(CRE,k,j,is);
          u_cr(CRF1,k,j,is-i) = vx*4.0*beta*presProfile(pco->x1v(is),x2) + u_cr(CRF1,k,j,is); //u_cr(CRF1,k,j,is);
          u_cr(CRF2,k,j,is-i) = u_cr(CRF2,k,j,is);
          u_cr(CRF3,k,j,is-i) = u_cr(CRF3,k,j,is); //u_cr(CRF3,k,j,is);

        }
      }
    }
  }
  return;
}
void ProfilesCROuterX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=1; i<=ngh; ++i) {
          Real x1 = pco->x1v(ie+i);
          Real x2 = pco->x2v(j);
          u_cr(CRE,k,j,ie+i) = 3.0*beta*(presProfile(x1,x2) - presProfile(pco->x1v(ie),x2)) + u_cr(CRE,k,j,ie);
          u_cr(CRF1,k,j,ie+i) = vx*4.0*beta*presProfile(pco->x1v(ie),x2) + u_cr(CRF1,k,j,ie); //u_cr(CRF1,k,j,is);
          u_cr(CRF2,k,j,ie+i) = u_cr(CRF2,k,j,ie);

          u_cr(CRF3,k,j,ie+i) = u_cr(CRF3,k,j,ie); //u_cr(CRF3,k,j,is);

        }
      }
    }
  }
  return;
}
void ProjectCRInnerX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real x1 = pco->x1v(i);
          Real x2 = pco->x2v(js-j);

          u_cr(CRE,k,js-j,i) = 3.0*beta*(presProfile(x1,x2)-presProfile(x1,pco->x2v(js))) + u_cr(CRE,k,js,i);
          u_cr(CRF1,k,js-j,i) = vx*4.0*beta*presProfile(x1,pco->x2v(js)) + u_cr(CRF1,k,js,i); //u_cr(CRF1,k,j,is);
          u_cr(CRF2,k,js-j,i) = u_cr(CRF2,k,js,i);
          u_cr(CRF3,k,js-j,i) = u_cr(CRF3,k,js,i); //u_cr(CRF3,k,j,is);

        }
      }
    }
  }
  return;
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
          Real x1 = pco->x1v(i);
          Real x2 = pco->x2v(je+j);

          u_cr(CRE,k,je+j,i) = 3.0*beta*(presProfile(x1,x2)-presProfile(x1,pco->x2v(je))) + u_cr(CRE,k,je,i);
          u_cr(CRF2,k,je+j,i) = vx*4.0*beta*presProfile(x1,pco->x2v(je)) + u_cr(CRF2,k,je,i); //u_cr(CRF1,k,j,is);
          u_cr(CRF1,k,je+j,i) = u_cr(CRF1,k,je,i);
          u_cr(CRF3,k,je+j,i) = u_cr(CRF3,k,je,i); //u_cr(CRF3,k,j,is);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
void myCRSource(MeshBlock *pmb, const Real time, const Real dt,
                const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
                AthenaArray<Real> &u_cr){
  //CR Source function
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;

  if (CR_ENABLED) {
    //if ((time >= crPertStartTime)&&(time<crPertEndTime)){
      Real myVal = 0.0;
      Real t1 = crPertStartTime;
      Real t2 = crPertEndTime;
      Real timeDep = (1-exp(-1*time/t1))*exp(-1*time/t2)*t1/t2*pow(t2/t1+1,t1/t2+1);
      for(int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            Real x1 = pmb->pcoord->x1v(i);
            Real x2 = pmb->pcoord->x2v(j);
            myVal = 0.0;

            Real dist = pow(pow(x1-crPertCenterX,2.0)+pow(x2-crPertCenterY,2.0),0.5);
            //Real dist = fabs(x2-crPertCenterY);

            myVal = 3.0*pres0*beta*crPertAmp*0.5*(1 - tanh(dist/crPertSteep));

            u_cr(CRE,k,j,i) += myVal*timeDep;
          }
        }
      }
    //}
  }
  return;
}



//----------------------------------------------------------------------------------------
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
