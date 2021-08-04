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
#include "../scalars/scalars.hpp"


//======================================================================================
/*! \file ParkerInst.cpp
 *  Magneto-hydrostatic equilibrium for local sim of a galactic disk (vertical structure)
 *  Perturbation is gaussian profile in CR Pressure
 *====================================================================================*/

//Hydrostatic Equilibrium variables
Real dens0, pres0, vx;
Real H, nGrav;
Real g0;
Real alpha;
Real beta;
Real dfloor, pfloor;
const Real float_min{std::numeric_limits<float>::min()};

//Cooling and perturbation scalar
int cooling, tracking;

//Perturbation variables
Real crPertCenterX;
Real crPertCenterY;
Real crPertAmp;
Real crPertSteep;

//Profile functions
Real densProfile(Real x1, Real x2);
Real presProfile(Real x1, Real x2);
Real gravProfile(Real x1, Real x2);
Real pertProfile(Real x1, Real x2);

//For logarithmic height spacing
Real LogMeshSpacingX2(Real x, RegionSize rs);

//Gravity function tanh, and cooling
void mySource(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons);

//CR LOCAL Source term function
//void myCRSource(MeshBlock *pmb, const Real time, const Real dt,
//                const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
//                AthenaArray<Real> &u_cr);


//Boundary condtions
//Hydro conditions
//X1 boundaries to match HSE
// void ProfilesOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
// void ProfilesInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

//x2 boundaries to extend HSE
void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//CR boundary conditions
//X1 match HSE on x1 bounds
// void ProfilesCROuterX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
//     const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//     AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
//     int js, int je, int ks, int ke, int ngh);
//
// void ProfilesCRInnerX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
//     const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//     AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
//     int js, int je, int ks, int ke, int ngh);

//Extend HSE at x2 bounds
void ProjectCROuterX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

void ProjectCRInnerX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

//cr Diffusion variables and function
Real sigma;
void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);


//Implement functions
Real densProfile(Real x1, Real x2)
{
  Real rho = dens0;
  Real prof = pow(cosh(x2/(nGrav*H)),-1.0*nGrav);
  rho *= prof;
  return rho;
}

Real presProfile(Real x1, Real x2)
{
  Real pres = pres0;
  Real rho = densProfile(x1,x2);
  pres *= rho/dens0;
  return pres;
}

Real gravProfile(Real x1, Real x2)
{
  Real g = g0 * tanh(x2/(nGrav*H));
  return g;
}

Real pertProfile(Real x1, Real x2)
{
  Real dist = pow(pow(x1-crPertCenterX,2.0)+pow(x2-crPertCenterY,2.0),0.5);
  Real p = crPertAmp/(crPertSteep*pow(2*M_PI,0.5))*exp(-0.5*pow(dist/crPertSteep,2.0));
  return p;
}

Real LogMeshSpacingX2(Real x, RegionSize rs)
{
  Real xf, xrat;
  xrat = pow(rs.x2max/rs.x2min,1.0/((Real) rs.nx2));
  xf = rs.x2min*pow(xrat,x*rs.nx2);
  return xf;
}

// Set local (each cell independently) opacity funciton
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  if(CR_ENABLED){
   pcr->EnrollOpacityFunction(Diffusion);
  }
}

//Set up initial MESH data
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  //Real x2rat = pin->GetOrAddReal("mesh","x2rat",0.0);
  //if (x2rat< 0.0) {
  //  EnrollUserMeshGenerator(X2DIR,LogMeshSpacingX2);
  //}
  // MHD boundary conditions
  if (pin->GetString("mesh","ix2_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x2, ProjectPressureInnerX2);
  }
  if (pin->GetString("mesh","ox2_bc")=="user"){
    EnrollUserBoundaryFunction(outer_x2, ProjectPressureOuterX2);
  }
  // if (pin->GetString("mesh","ix1_bc")=="user"){
  //   EnrollUserBoundaryFunction(inner_x1, ProfilesInnerX1);
  // }
  // if (pin->GetString("mesh","ox1_bc")=="user"){
  //   EnrollUserBoundaryFunction(outer_x1, ProfilesOuterX1);
  // }

  // Source Functions
  EnrollUserExplicitSourceFunction(mySource);
  cooling = pin->GetOrAddInteger("problem","cooling",0);

  if(CR_ENABLED){
    //CR Boundary conditions
    if (pin->GetString("mesh","ix2_bc")=="user")
      EnrollUserCRBoundaryFunction(inner_x2, ProjectCRInnerX2);
    if (pin->GetString("mesh","ox2_bc")=="user")
      EnrollUserCRBoundaryFunction(outer_x2, ProjectCROuterX2);
    // if (pin->GetString("mesh","ix1_bc")=="user")
    //   EnrollUserCRBoundaryFunction(inner_x1, ProfilesCRInnerX1);
    // if (pin->GetString("mesh","ox1_bc")=="user")
    //   EnrollUserCRBoundaryFunction(outer_x1, ProfilesCROuterX1);
  }
}

//Setup initial mesh and variables
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Real gamma = peos->GetGamma();

  // Load variables
  vx=pin->GetReal("problem","xVel");
  nGrav = pin->GetReal("problem","GravNumScaleHeight");
  beta = pin->GetOrAddReal("problem","beta",0.0);
  alpha = pin->GetOrAddReal("problem","alpha",0.0);
  dens0 = pin->GetReal("problem","Dens");
  g0 = pin->GetReal("problem","grav");
  pres0 = pin->GetReal("problem","Pres");
  dfloor = pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*float_min)) ;
  pfloor = pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*float_min)) ;


  // Derived variables
  H = pres0/dens0/(-1.0*g0)*(1+alpha+beta);

  if(CR_ENABLED){
    //Load CR Variables
    crPertCenterX = pin->GetReal("problem","pertX");
    crPertCenterY = pin->GetReal("problem","pertY");
    sigma = pin->GetReal("cr","sigma");
    crPertAmp = pin->GetReal("problem","pertAmplitude");
    crPertSteep = pin->GetReal("problem","pertDx");
  }
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real density = densProfile(x1,x2);
        Real pressure = presProfile(x1,x2);

        //set hydro variables
        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = vx*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pressure/(gamma-1) + 0.5*density*pow(vx,2.0);

        if(CR_ENABLED){
          // get CR parameters
          Real crp = beta*pressure; //*(1+ampCR*(1-x2/centerCR));
          Real pertVal = pertProfile(x1,x2);
          Real pertValNoX = pertProfile(crPertCenterX,x2);
          Real dPcdz = -1.0*beta*presProfile(x1,x2)/H*tanh(x2/(nGrav*H));

          // set CR variables
          pcr->u_cr(CRE,k,j,i) = 3.0*crp+3.0*beta*pres0*pertVal;
          pcr->u_cr(CRF1,k,j,i) = vx*4.0*crp;
          pcr->u_cr(CRF2,k,j,i) = -1.0*dPcdz/sigma;
          pcr->u_cr(CRF3,k,j,i) = 0.0;

          // Setup scalar tracker for perturbation
          if ((NSCALARS > 0) ) {
            pscalars->s(0,k,j,i) = pertVal/crPertAmp;
            pscalars->s(1,k,j,i) = pertValNoX/crPertAmp;
          }
        }
      }// end i
    }
  }
    // initialize uniform interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0   = pow(2*alpha,0.5);
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          Real x2 = pcoord->x2v(j);
          Real x1 = pcoord->x1f(i);
          Real pressure= presProfile(x1,x2);

          Real B = B0*sqrt(pressure);

          pfield->b.x1f(k,j,i) = B;


        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          Real x2 = pcoord->x2f(j);
          Real x1 = pcoord->x1v(i);
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

            phydro->u(IEN,k,j,i) += 0.5*( pow(pfield->b.x1f(k,j,i),2.0)
                                         +pow(pfield->b.x2f(k,j,i),2.0)
                                         +pow(pfield->b.x3f(k,j,i),2.0));
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
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2v(j);
          pcr->sigma_diff(0,k,j,i) = sigma;
          pcr->sigma_diff(1,k,j,i) = sigma;
          pcr->sigma_diff(2,k,j,i) = sigma;
        }
      }
    }// end k,j,i

  }// End CR
  return;
}

//Source function with gravity and cooling
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
        //GRAVITY
        Real x2 = pmb->pcoord->x2v(j);
        Real x1 = pmb->pcoord->x1v(i);
        Real gravity = gravProfile(x1,x2);
        Real src = dt*prim(IDN,k,j,i)*gravity;

        cons(IM2,k,j,i) += src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);

        //COOLING
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


// //----------------------------------------------------------------------------------------
// //! \fn ProfilesInnerX1()
// //  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm
// void ProfilesInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
// {
//   for (int k=ks; k<=ke; ++k) {
//     for (int j=js; j<=je; ++j) {
// #pragma omp simd
//       for (int i=1; i<=ngh; ++i) {
//           Real x2 = pco->x2v(j);
//           Real x1 = pco->x1v(is-i);
//           Real density = densProfile(x1,x2);
//           Real pressure = presProfile(x1,x2);
//
//           Real Diff = prim(IDN,k,j,is-1+i) - densProfile(pco->x1v(is-1+i),x2);
//           prim(IDN,k,j,is-i) = density - Diff;
//
//           Diff = prim(IVX,k,j,is-1+i) - vx;
//           prim(IVX,k,j,is-i) = vx - Diff;
//           prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
//           prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
//
//           Diff = prim(IPR,k,j,is-1+i) - presProfile(pco->x1v(is-1+i),x2);
//           prim(IPR,k,j,is-i) = pressure-Diff;
//
//       }
//     }
//   }
//    // set magnetic field in inlet ghost zones
//   if (MAGNETIC_FIELDS_ENABLED) {
//     Real B0   = pow(2*alpha,0.5);
//     for(int k=ks; k<=ke; ++k){
//     for(int j=js; j<=je; ++j){
// #pragma omp simd
//       for(int i=1; i<=ngh; ++i){
//         Real x2 = pco->x2v(j);
//         Real x1 = pco->x1f(is-i);
//
//         Real B = B0 * sqrt(presProfile(x1,x2));
//         Real Diff = b.x1f(k,j,is-1+i) - B0*sqrt(presProfile(pco->x1v(is-1+i),x2));
//
//         b.x1f(k,j,is-i) = B - Diff;
//       }
//     }}
//
//     for(int k=ks; k<=ke; ++k){
//     for(int j=js; j<=je+1; ++j){
// #pragma omp simd
//       for(int i=1; i<=ngh; ++i){
//         b.x2f(k,j,is-i) = b.x2f(k,j,is-1+i);
//       }
//     }}
//
//     for(int k=ks; k<=ke+1; ++k){
//     for(int j=js; j<=je; ++j){
// #pragma omp simd
//       for(int i=1; i<=ngh; ++i){
//         b.x3f(k,j,is-i) = b.x3f(k,j,is-1+i);
//       }
//     }}
//
//   }
//
//
//   return;
// }
//
// //----------------------------------------------------------------------------------------
// //! \fn ProfilesOuterX1()
// //  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm
// void ProfilesOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
//       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
// {
//   for (int k=ks; k<=ke; ++k) {
//     for (int j=js; j<=je; ++j) {
// #pragma omp simd
//       for (int i=1; i<=ngh; ++i) {
//           Real x2 = pco->x2v(j);
//           Real x1 = pco->x1v(ie+i);
//           Real density = densProfile(x1,x2);
//           Real pressure = presProfile(x1,x2);
//
//           Real Diff = prim(IDN,k,j,ie+1-i) - densProfile(pco->x1v(ie+1-i),x2);
//           prim(IDN,k,j,ie+i) = density - Diff;
//
//           Diff = prim(IVX,k,j,ie+1-i) - vx;
//           prim(IVX,k,j,ie+i) = vx - Diff;
//           prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
//           prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
//
//           Diff = prim(IPR,k,j,ie+1-i) - presProfile(pco->x1v(ie+1-i),x2);
//           prim(IEN,k,j,ie+i) = pressure - Diff;
//
//       }
//     }
//   }
//    // set magnetic field in inlet ghost zones
//   if (MAGNETIC_FIELDS_ENABLED) {
//     Real B0   = pow(2*alpha,0.5);
//     for(int k=ks; k<=ke; ++k){
//     for(int j=js; j<=je; ++j){
// #pragma omp simd
//       for(int i=2; i<=ngh+1; ++i){
//         Real x2 = pco->x2v(j);
//         Real x1 = pco->x1f(ie+i);
//         Real B = B0 * sqrt(presProfile(x1,x2));
//         Real Diff = b.x1f(k,j,ie+2-i) - B0*sqrt(presProfile(pco->x1v(ie+2-i),x2));
//         b.x1f(k,j,ie+i) = B-Diff;
//       }
//     }}
//
//     for(int k=ks; k<=ke; ++k){
//     for(int j=js; j<=je+1; ++j){
// #pragma omp simd
//       for(int i=1; i<=ngh; ++i){
//         b.x2f(k,j,ie+i) = b.x2f(k,j,ie+1-i);
//       }
//     }}
//
//     for(int k=ks; k<=ke+1; ++k){
//     for(int j=js; j<=je; ++j){
// #pragma omp simd
//       for(int i=1; i<=ngh; ++i){
//         b.x3f(k,j,ie+i) = b.x3f(k,j,ie+1 - i);
//       }
//     }}
//
//   }
//
//
//   return;
// }
//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real dx = pco->x2v(jl+1) - pco->x2v(jl);
        Real delta = pco->x2v(jl-j) - pco->x2v(jl);
        Real dy = prim(IDN,k,jl+1,i) - prim(IDN,k,jl,i);
        Real myVal = dy/dx*delta + prim(IDN,k,jl,i);
        if (myVal < dfloor) {
          prim(IDN,k,jl-j,i) = dfloor;
        } else {
          prim(IDN,k,jl-j,i) = myVal;
        }

        prim(IVX,k,jl-j,i) = 0.0;//vx- Diff;
        prim(IVY,k,jl-j,i) = prim(IVY,k,jl,i);//prim(IVY,k,jl-1+j,i);
        prim(IVZ,k,jl-j,i) = 0.0;//prim(IVZ,k,jl-1+j,i);

        dy = prim(IPR,k,jl+1,i) - prim(IPR,k,jl,i); //prim(IPR,k,jl-1+j,i) - presProfile(x1,pco->x2v(jl-1+j));
        myVal = dy/dx*delta + prim(IPR,k,jl,i);
        if (myVal < pfloor) {
          prim(IPR,k,jl-j,i) = pfloor;
        } else {
          prim(IPR,k,jl-j,i) = myVal;
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0   = pow(2*alpha,0.5);
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          //Real x2 = pco->x2v(jl-j);
          //Real x1 = pco->x1f(i);
          //Real B = B0*sqrt(presProfile(x1,x2));
          Real dx = pco->x2v(jl+1) - pco->x2v(jl);
          Real delta = pco->x2v(jl-j) - pco->x2v(jl);
          Real dy = b.x1f(k,jl+1,i) - b.x1f(k,jl,i);
          Real myVal = dy/dx*delta + b.x1f(k,jl,i) ;
          if (myVal*b.x1f(k,jl,i) <= 0) {
            myVal = 0.0;
          }
          //Real Diff = b.x1f(k,jl-1+j,i) - B0*sqrt(presProfile(x1,pco->x1v(jl-1+j)));
          b.x1f(k,(jl-j),i) =  myVal;// - Diff;
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(jl-j),i) = 0.0;//b.x2f(k,jl-1+j,i);// std::pow(alpha*pressure,0.5);  // reflect 2-field
        }
      }
    }
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(jl-j),i) = 0.0;//b.x3f(k,jl-1+j,i);// b.x3f(k,(jl+j-1),i);
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
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        //Real x1 = pco->x1v(i);
        //Real x2 = pco->x2v(ju+j);
        //Real density = densProfile(x1,x2);
        //Real pressure = presProfile(x1,x2);
        Real dx = pco->x2v(ju-1) - pco->x2v(ju);
        Real delta = pco->x2v(ju+j) - pco->x2v(ju);
        Real dy = prim(IDN,k,ju-1,i) - prim(IDN,k,ju,i);
        Real myVal = dy/dx*delta + prim(IDN,k,ju,jl);
        if (myVal < dfloor) {
          prim(IDN,k,ju+j,i) = dfloor;
        } else {
          prim(IDN,k,ju+j,i) = myVal;
        }

        //Diff = prim(IVX,k,ju-j+1,i) - vx;
        prim(IVX,k,ju+j,i) = 0.0;//vx - Diff;
        prim(IVY,k,ju+j,i) = prim(IVY,k,ju,i);//prim(IVY,k,ju-j+1,i);  // reflect 2-velocity
        prim(IVZ,k,ju+j,i) = 0.0;//prim(IVZ,k,ju-j+1,i);

        //Diff = prim(IPR,k,ju-j+1,i) - presProfile(x1,pco->x2v(ju-j+1));
        dy = prim(IPR,k,ju-1,i) - prim(IPR,k,ju,i);
        myVal = dy/dx*delta + prim(IPR,k,ju,jl);
        if (myVal < pfloor) {
          prim(IPR,k,ju+j,i) = pfloor;
        } else {
          prim(IPR,k,ju+j,i) = myVal;
        }

      }
    }
  }


  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0   = pow(2*alpha,0.5);
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          //Real x2 = pco->x2v(ju+j);
          //Real x1 = pco->x1f(i);
          //Real B = B0*sqrt(presProfile(x1,x2));
          Real dx = pco->x2v(ju-1) - pco->x2v(ju);
          Real delta = pco->x2v(ju+j) - pco->x2v(ju);
          Real dy = b.x1f(k,ju-1,i) - b.x1f(k,ju,i);
          Real myVal = dy/dx*delta + b.x1f(k,ju,i) ;
          if (myVal*b.x1f(k,ju,i) <= 0) {
            myVal = 0.0;
          }
          //Real Diff = b.x1f(k,ju+1-j,i) - B0*sqrt(presProfile(x1,pco->x1v(ju+1-j)));
          b.x1f(k,(ju+j),i) =  myVal;//- Diff;
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=2; j<=ngh+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(ju+j),i) = 0.0;//b.x2f(k,ju+1-j,i) ;// std::pow(alpha*pressure,0.5);  // reflect 2-field
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(ju+j),i) =  0.0;//b.x3f(k,ju-j+1,i);
        }
      }
    }

  }
  return;
}

//======================================================================================
// CR Boundary Conditions
//======================================================================================
// void ProfilesCRInnerX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
//     const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//     AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
//     int js, int je, int ks, int ke, int ngh)
// {
//   if(CR_ENABLED){
//     for (int k=ks; k<=ke; ++k) {
//       for (int j=js; j<=je; ++j) {
//         for (int i=1; i<=ngh; ++i) {
//           Real x1 = pco->x1v(is-i);
//           Real x2 = pco->x2v(j);
//           Real dPcdz = -1.0*beta*presProfile(x1,x2)/H*tanh(x2/(nGrav*H));
//           Real crp = beta * presProfile(x1,x2);
//           Real crfx = 4.0 * crp * vx;
//           Real crfy = -1.0*dPcdz/sigma;
//
//           Real Diff = u_cr(CRE,k,j,is-1+i) - 3.0*beta*presProfile(pco->x1v(is-1+i),x2);
//           u_cr(CRE,k,j,is-i) = 3.0*crp - Diff;
//
//           Diff = u_cr(CRF1,k,j,is-1+i) - 4.0 * vx * beta*presProfile(pco->x1v(is-1+i),x2);
//           u_cr(CRF1,k,j,is-i) = crfx - Diff;
//
//           Diff = u_cr(CRF2,k,j,is-1+i) + beta*presProfile(pco->x1v(is-1+i),x2)/H*tanh(x2/(nGrav*H)) / sigma;
//           u_cr(CRF2,k,j,is-i) = crfy - Diff;
//
//           u_cr(CRF3,k,j,is-i) = 0.0;
//
//         }
//       }
//     }
//   }
//   return;
// }
//
// void ProfilesCROuterX1(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
//     const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//     AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
//     int js, int je, int ks, int ke, int ngh)
// {
//   if(CR_ENABLED){
//     for (int k=ks; k<=ke; ++k) {
//       for (int j=js; j<=je; ++j) {
//         for (int i=1; i<=ngh; ++i) {
//           Real x1 = pco->x1v(ie+i);
//           Real x2 = pco->x2v(j);
//           Real dPcdz = -1.0*beta*presProfile(x1,x2)/H*tanh(x2/(nGrav*H));
//           Real crp = beta * presProfile(x1,x2);
//           Real crfx = 4.0 * crp * vx;
//           Real crfy = -1.0*dPcdz/sigma;
//
//           Real Diff = u_cr(CRE,k,j,ie+1-i) - 3.0*beta*presProfile(pco->x1v(ie+1-i),x2);
//           u_cr(CRE,k,j,ie+i) = 3.0*crp - Diff;
//
//           Diff = u_cr(CRF1,k,j,ie+1-i) - 4.0 * vx * beta*presProfile(pco->x1v(ie+1-i),x2);
//           u_cr(CRF1,k,j,ie+i) = crfx - Diff; //u_cr(CRF1,k,j,is);
//
//           Diff = u_cr(CRF2,k,j,ie+1-i) + beta*presProfile(pco->x1v(ie+1-i),x2)/H*tanh(x2/(nGrav*H)) / sigma;
//           u_cr(CRF2,k,j,ie+i) = crfy - Diff;
//
//           u_cr(CRF3,k,j,ie+i) = 0.0; //u_cr(CRF3,k,j,is);
//
//         }
//       }
//     }
//   }
//   return;
// }

void ProjectCRInnerX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real dx = pco->x2v(js+1) - pco->x2v(js);
          Real delta = pco->x2v(js-j) - pco->x2v(js);
          Real dy =   u_cr(CRE,k,js+1,i) -   u_cr(CRE,k,js,i);
          Real myVal = dy/dx*delta + u_cr(CRE,k,js,i);
          if (myVal < 3.0*pfloor) {
            u_cr(CRE,k,js-j,i) = 3.0*pfloor;
          } else {
            u_cr(CRE,k,js-j,i) = myVal;
          }
          u_cr(CRF1,k,js-j,i) = 0.0;
          u_cr(CRF2,k,js-j,i) = u_cr(CRF2,k,js,i);//beta*presProfile(x1,x2)/(sigma*H)*tanh(x2/(nGrav*H));
          u_cr(CRF3,k,js-j,i) = 0.0; //u_cr(CRF3,k,j,is);

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
          Real dx = pco->x2v(je-1) - pco->x2v(je);
          Real delta = pco->x2v(je+j) - pco->x2v(je);
          Real dy = u_cr(CRE,k,je-1,i) - u_cr(CRE,k,je,i);
          Real myVal = dy/dx*delta + u_cr(CRE,k,je,i);
          if (myVal < 3.0*pfloor) {
            u_cr(CRE,k,je+j,i) = 3.0*pfloor;
          } else {
            u_cr(CRE,k,je+j,i) = myVal;
          }
          u_cr(CRF1,k,je+j,i) = 0.0;//crfx - Diff;
          u_cr(CRF2,k,je+j,i) = u_cr(CRF2,k,je,i);//crfy - Diff;
          u_cr(CRF3,k,je+j,i) = 0.0; //u_cr(CRF3,k,j,is);
        }
      }
    }
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
  Coordinates *pco = pmb->pcoord;
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
