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
//NOTE x1 is (x) direction, along magnetic field
//     x2 is (z) direction, along gravitational field
//     x3 is (y) direction, perpendicular to both, radial in galactic disk coords
//     Length scale is H = P_0/(rho_0 g_0) (1+ alpha + beta)
//
//Hydrostatic Equilibrium variables
Real dens0, pres0, vx; // Initial hydro quantities
                       // dens0 in multiples of 1e-24 g/cm^3
                       // pres0 in multiples of 1e-12 erg/cm^3
                       // vx in multiples of 1e6 cm/s
Real H, nGrav; // H is scalar multiple change from
               // H_0 = Sqr(1e6 cm/s) / (4 * 1e-9 cm/s^2)
               //     = 0.25e21 cm = 81.019482 pc
               // nGrav is scale height of stars divided by scale height of gas
               //      approx 1 in Milky Way
Real g0;       // in multiples of 4e-9 cm/s^2
Real alpha;    // Ratio of magnetic pressure to gas pressure
Real beta;     // Ratio of cosmic ray pressure to gas pressure

//Floors for Diode boundary conds
Real dfloor, pfloor; // Floor values for density and rpessure

const Real float_min{std::numeric_limits<float>::min()};
//Useful minimum floating point number

//Cooling and perturbation scalar
int cooling; //Boolean - if cooling==1 do Inoue 2006 2 phase gas cooling profile


//Perturbation variables
Real crPertCenterZ; // pert is at X=0, Y=0 this determines height from disk
Real crD; // percentage of SN energy in CRs - in multiples of 10%
Real crEsn; // energy of SN in units of 1e51 erg
Real crPertRad; // radius of supernova expansion (width of gaussian profile)
                // in units of 10pc

//Profile functions
Real densProfile(Real x1, Real x2, Real x3);
Real presProfile(Real x1, Real x2, Real x3);
Real gravProfile(Real x1, Real x2, Real x3);
Real pertProfile(Real x1, Real x2, Real x3);

//For logarithmic height spacing
Real LogMeshSpacingX2(Real x, RegionSize rs);

//Gravity function tanh, and cooling
void mySource(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons);

//x2 boundaries with vacuum
void DiodeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiodeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//CR boundary conditions

// vacuum at x2 bounds
void DiodeCROuterX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

void DiodeCRInnerX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh);

//cr Diffusion variables and function
Real sigma;
void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);


//Implement functions
Real densProfile(Real x1, Real x2, Real x3)
{
  Real rho = pow(cosh(x2/(nGrav)),-1.0*nGrav);
  return rho;
}

Real presProfile(Real x1, Real x2, Real x3)
{
  Real pres = densProfile(x1,x2,x3);
  return pres;
}

Real gravProfile(Real x1, Real x2, Real x3)
{
  Real g = -1*tanh(x2/(nGrav));//*g0;
  return g;
}

Real pertProfile(Real x1, Real x2, Real x3)
{
  Real dist = pow(SQR(x1)+SQR(x2-crPertCenterZ)+SQR(x3),0.5);
  Real p = pow(crPertRad,-3.0)*exp(-32.82*SQR(dist/crPertRad)*SQR(H));
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
  Real x2rat = pin->GetOrAddReal("mesh","x2rat",0.0);
  if (x2rat< 0.0) {
    EnrollUserMeshGenerator(X2DIR,LogMeshSpacingX2);
  }
  // MHD boundary conditions
  if (pin->GetString("mesh","ix2_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x2, DiodeInnerX2);
  }
  if (pin->GetString("mesh","ox2_bc")=="user"){
    EnrollUserBoundaryFunction(outer_x2, DiodeOuterX2);
  }

  // Source Functions
  EnrollUserExplicitSourceFunction(mySource);
  cooling = pin->GetOrAddInteger("problem","cooling",0);

  if(CR_ENABLED){
    //CR Boundary conditions
    if (pin->GetString("mesh","ix2_bc")=="user")
      EnrollUserCRBoundaryFunction(inner_x2, DiodeCRInnerX2);
    if (pin->GetString("mesh","ox2_bc")=="user")
      EnrollUserCRBoundaryFunction(outer_x2, DiodeCROuterX2);
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
  H = pres0/dens0/(g0)*(1+alpha+beta);

  if(CR_ENABLED){
    //Load CR Variables
    crPertCenterZ = pin->GetReal("problem","pertZ");
    sigma = pin->GetReal("cr","sigma");
    crEsn= pin->GetReal("problem","snEner");
    crD= pin->GetReal("problem","snEnerFrac");
    crPertRad = pin->GetReal("problem","pertR");
  }
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real density = densProfile(x1,x2,x3);
        Real pressure = presProfile(x1,x2,x3);

        //set hydro variables
        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = vx*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pressure/(gamma-1) + 0.5*density*SQR(vx);

        if(CR_ENABLED){
          // get CR parameters
          Real crp = beta*pressure; //*(1+ampCR*(1-x2/centerCR));
          Real pertVal = pertProfile(x1,x2,x3);
          Real pertValNoX = pertProfile(0.0,x2,x3);
          Real dPcdz = -1.0*beta*presProfile(x1,x2,x3)*tanh(x2/(nGrav));

          // set CR variables
          pcr->u_cr(CRE,k,j,i) = 3.0*crp+pertVal * crD * crEsn * 216.1118 / pres0;
          pcr->u_cr(CRF1,k,j,i) = vx*4.0*crp;
          pcr->u_cr(CRF2,k,j,i) = -1.0*dPcdz/sigma;
          pcr->u_cr(CRF3,k,j,i) = 0.0;

          // Setup scalar tracker for perturbation
          if ((NSCALARS > 0) ) {
            pscalars->s(0,k,j,i) = pertVal;
            pscalars->s(1,k,j,i) = pertValNoX;
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
          Real x3 = pcoord->x3v(k);
          Real pressure= presProfile(x1,x2,x3);

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
        Real x3 = pmb->pcoord->x3v(k);
        Real gravity = gravProfile(x1,x2,x3);
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


//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void DiodeInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        prim(IDN,k,jl-j,i) = dfloor;
        prim(IPR,k,jl-j,i) = pfloor;
        prim(IVX,k,jl-j,i) = 0.0;
        prim(IVY,k,jl-j,i) = 0.0;
        prim(IVZ,k,jl-j,i) = 0.0;
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
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
          b.x2f(k,(jl-j),i) = 0.0;
        }
      }
    }
    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(jl-j),i) = 0.0;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureOuterX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void DiodeOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh)
{
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        prim(IDN,k,ju+j,i) = dfloor;
        prim(IPR,k,ju+j,i) = pfloor;
        prim(IVX,k,ju+j,i) = 0.0;
        prim(IVY,k,ju+j,i) = 0.0;
        prim(IVZ,k,ju+j,i) = 0.0;
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(ju+j),i) =  0.0;
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=2; j<=ngh+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(ju+j),i) = 0.0;
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
  return;
}

//======================================================================================
// CR Boundary Conditions
//======================================================================================
void DiodeCRInnerX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          u_cr(CRE,k,js-j,i) = 3.0*pfloor;
          u_cr(CRF1,k,js-j,i) = 0.0;
          u_cr(CRF2,k,js-j,i) = 0.0;
          u_cr(CRF3,k,js-j,i) = 0.0;

        }
      }
    }
  }
  return;
}

void DiodeCROuterX2(MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u_cr, Real time, Real dt, int is, int ie,
    int js, int je, int ks, int ke, int ngh)
{
  if(CR_ENABLED){
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=is; i<=ie; ++i) {
          u_cr(CRE,k,je+j,i) = 3.0*pfloor;
          u_cr(CRF1,k,je+j,i) = 0.0;
          u_cr(CRF2,k,je+j,i) = 0.0;
          u_cr(CRF3,k,je+j,i) = 0.0;
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
