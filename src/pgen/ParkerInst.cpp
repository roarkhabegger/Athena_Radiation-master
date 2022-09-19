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
Real dens0, pres0, g0, vx, myGamma; // Initial hydro quantities
                       // dens0 in multiples of 1e-24 g/cm^3
                       // pres0 in multiples of 1e-12 erg/cm^3
                       // g0 = (1+alpha + beta)*pres0/dens0
                       // vx in multiples of 1e6 cm/s
Real nGrav;    // nGrav is scale height of stars divided by scale height of gas
               //      approx 1 in Milky Way
Real alpha;    // Ratio of magnetic pressure to gas pressure
Real beta;     // Ratio of cosmic ray pressure to gas pressure

Real h;

//Floors for Diode boundary conds
Real dfloor, pfloor; // Floor values for density and rpessure

const Real float_min{std::numeric_limits<float>::min()};
//Useful minimum floating point number

//Cooling and perturbation scalar
int cooling; //Boolean - if cooling==1 do Inoue 2006 2 phase gas cooling profile


//Perturbation variables
Real crPertCenterX;
Real crPertCenterY;
Real crPertCenterZ; // this determines height from disk
Real crD; // percentage of SN energy in CRs - in multiples of 10%
Real crEsn; // energy of SN in units of 1e51 erg
Real crPertRad; // radius of supernova expansion (width of gaussian profile)
                // in units of 10pc
Real crThermal;
Real crLinear;
Real randAmplitude;
int XNmax;
int YNmax;
Real xRange;
Real yRange;
Real *randsX;
Real *randsY;


Real s1Y;
Real s1Z;
Real s1R;
Real s1dR;

Real s2Y;
Real s2Z;
Real s2R;
Real s2dR;

//Profile functions
Real densProfile(Real x1, Real x2, Real x3);
Real presProfile(Real x1, Real x2, Real x3);
Real gravProfile(Real x1, Real x2, Real x3);
Real potProfile(Real x1, Real x2, Real x3);
// Real pertProfile(Real x1, Real x2, Real x3);
Real fcProfile(Real x1, Real x2, Real x3);
Real s1Profile(Real x1, Real x2, Real x3);
Real s2Profile(Real x1, Real x2, Real x3);

//For logarithmic height spacing
Real LogMeshSpacingX2(Real x, RegionSize rs);

//Gravity function tanh, and cooling
void mySource(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons);
void CRSource(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &u_cr);

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
Real sigmaParl, sigmaPerp;

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);


//Implement functions
Real densProfile(Real x1, Real x2, Real x3)
{
  Real rho = dens0*pow(cosh(x2/(nGrav*h)),-1.0*nGrav);
  return rho;
}

Real presProfile(Real x1, Real x2, Real x3)
{
  Real pres = pres0*pow(cosh(x2/(nGrav*h)),-1.0*nGrav);
  return pres;
}

Real gravProfile(Real x1, Real x2, Real x3)
{
  Real g = -1*g0*tanh(x2/(nGrav*h));//*g0;
  return g;
}
// Real potProfile(Real x1, Real x2, Real x3)
// {
//   Real phi = g0*nGrav*h*log(abs( cosh(x2/(nGrav*h)) ));//*g0;
//   return phi;
// }

Real pertProfile(Real x1, Real x2, Real x3)
{
  Real dist = pow(SQR(x1-crPertCenterX)+SQR(x2-crPertCenterZ)+SQR(x3-crPertCenterY),0.5);
  Real p = exp(-50.0*SQR(dist/crPertRad)) / pow(2*M_PI*SQR(crPertRad),1.5) * 1000;
  return p;
  //Coefficient is (100pc)^2/200 pc^2
}
Real fcProfile(Real x1, Real x2, Real x3)
{
  Real  fcz = pow(cosh(x2/(nGrav*h)),-1.0*nGrav)*tanh(x2/(nGrav*h));
  fcz *= beta*pres0/(h*sigmaPerp) ;
  return fcz;
}


Real s1Profile(Real x1, Real x2, Real x3)
{
  Real dist = pow(SQR(x2-s1Z)+SQR(x3-s1Y),0.5);
  Real p = 0.5*(1 - tanh((dist - s1R)/s1dR));
  //Real p = pow(s1R,-3.0)*pow(2*M_PI,1.5)*exp(-0.5*SQR(dist/s1R));
  return p;
}

Real s2Profile(Real x1, Real x2, Real x3)
{
  Real dist = pow(SQR(x2-s2Z)+SQR(x3-s2Y),0.5);
  Real p = 0.5*(1 - tanh((dist - s2R)/s2dR));
  // Real p = pow(s2R,-3.0)*pow(2*M_PI,1.5)*exp(-0.5*SQR(dist/s2R));
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
   pcr->EnrollUserCRSource(CRSource);
  }
}

//Set up initial MESH data
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  myGamma = pin->GetReal("hydro","gamma");

  // Load variables
  vx=pin->GetReal("problem","xVel");
  nGrav = pin->GetReal("problem","GravNumScaleHeight");
  beta = pin->GetOrAddReal("problem","beta",0.0);
  alpha = pin->GetOrAddReal("problem","alpha",0.0);
  pres0 = pin->GetReal("problem","Pres");

  dens0 = pin->GetReal("problem","Dens");

  g0 =pin->GetReal("problem","Grav");
  h  =pin->GetReal("problem","ScaleH");

  dfloor = pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*float_min)) ;
  pfloor = pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*float_min)) ;

  if(CR_ENABLED){
    //Load CR Variables
    crPertCenterX = pin->GetReal("problem","pertX");
    crPertCenterY = pin->GetReal("problem","pertY");
    crPertCenterZ = pin->GetReal("problem","pertZ");
    sigmaPerp = pin->GetReal("cr","sigmaPerp");
    sigmaParl = pin->GetReal("cr","sigmaParl");
    crEsn = pin->GetReal("problem","snEner");
    crD = pin->GetReal("problem","snEnerFrac");
    crPertRad = pin->GetReal("problem","pertR");
    crThermal = pin->GetOrAddReal("problem","ThermalBlast",-1);
    crLinear = pin->GetOrAddReal("problem","LinearPert",-1);
    if (crLinear > 0.0) {
      randAmplitude = pin->GetReal("problem", "randAmplitude");
      XNmax = (int) pin->GetOrAddReal("problem","XNMax",floor(pin->GetReal("mesh", "nx1")/5));
      YNmax = (int) pin->GetOrAddReal("problem","YNMax",floor(pin->GetReal("mesh", "nx3")/5));
    }
    // std::cout << crEsn << "  "<< crD << std::endl;
  }
  // Setup scalar tracker for perturbation
  if ((NSCALARS > 0) ) {
    s1Y = pin->GetOrAddReal("problem","scalar1Y",0.0);
    s1Z = pin->GetOrAddReal("problem","scalar1Z",0.0);
    s1R = pin->GetOrAddReal("problem","scalar1R",1.0);
    s1dR = pin->GetOrAddReal("problem","scalar1dR",0.000001);
    s2Y = pin->GetOrAddReal("problem","scalar2Y",0.0);
    s2Z = pin->GetOrAddReal("problem","scalar2Z",0.0);
    s2R = pin->GetOrAddReal("problem","scalar2R",1.0);
    s2dR = pin->GetOrAddReal("problem","scalar2dR",0.000001);
  }

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
  //From Sherry's ParkerInst_Perturb
  if (crLinear > 0.0) {
    //setup perturbation parameters before loop
    Real A = randAmplitude; //CHANGE TO COMPUTATIONAL UNITS (10^-12:?)
    //minimum 2
    xRange = pin->GetReal("mesh", "x1max") - pin->GetReal("mesh", "x1min");
    yRange = pin->GetReal("mesh", "x3max") - pin->GetReal("mesh", "x3min");
    //srand(gid); //arbitrary seed for each meshblock
    srand(10); //consistent seed
    //setup random phases for each wavelength
    randsX = (Real*) malloc(sizeof(double[XNmax]));
    randsY = (Real*) malloc(sizeof(double[YNmax]));
    //Real randsY[YNmax];
    for (int x=0; x<XNmax; x++) randsX[x] = (rand() * 2 * M_PI) / RAND_MAX;
    for (int y=0; y<YNmax; y++) randsY[y] = (rand() * 2 * M_PI) / RAND_MAX;
  }


  // Derived variables
  //H = pres0/(dens0*g0)*(1+alpha+beta); // H ==1 always if length scale is H
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real density = densProfile(x1,x2,x3);
        Real pressure = presProfile(x1,x2,x3);
        // Real pot = potProfile(x1,x2,x3)*density;

        //set hydro variables
        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = vx*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pressure/(myGamma-1) + 0.5*density*SQR(vx) ;//+ pot;

        //FROM Sherry's ParkerInst_Perturb
        if(crLinear > 0.0){
          //set perturbations in z (vertical) velocities
          Real dv = 0; //init & reset every loop
          Real A = randAmplitude;
          //3D case; sum over these terms
          for (int x=0; x<XNmax; x++){
            if(YNmax==0){
              // std::cout << "XNMAX IS" << std::endl;
              // std::cout << XNmax << std::endl;
              // std::cout << "YNMAX IS" << std::endl;
              // std::cout << YNmax << std::endl;
              Real Lx = xRange/(x+1); // x;
              dv += (A / (XNmax)) //amplitude scaled
                    * sin(((2*M_PI*x1) / Lx) - randsX[x]);

            }else{
              for (int y=0; y<YNmax; y++){
                Real Lx = xRange/(x+1); // x;
                Real Ly = yRange/(y+1); // y;
                dv += (A / (XNmax*YNmax)) //amplitude scaled
                      * sin(((2*M_PI*x1) / Lx) - randsX[x])
                      * sin(((2*M_PI*x3) / Ly) - randsY[y]);
              }
            }
          }
          //change momentum
          phydro->u(IM2, k, j, i) += dv*density;
          phydro->u(IEN, k, j, i) += 0.5*density*SQR(dv);
        }

        if(CR_ENABLED){
          // get CR parameters
          Real crp = beta*pressure; //*(1+ampCR*(1-x2/centerCR));
          Real pertVal = pertProfile(x1,x2,x3);;
          pcr->u_cr(CRE,k,j,i) = 3.0*crp;
          pcr->u_cr(CRF1,k,j,i) = vx*4.0*crp;
          pcr->u_cr(CRF2,k,j,i) = fcProfile(x1,x2,x3);//-1.0*dPcdz/sigmaParl;
          pcr->u_cr(CRF3,k,j,i) = 0.0;
          // set CR variables
          if (crEsn > 0.0) {
            if (crThermal < 0.0) {
              pcr->u_cr(CRE,k,j,i) += pertVal * crD * crEsn;
            } else {
              phydro->u(IEN,k,j,i) += pertVal * crD * crEsn;
            }
          }
          //perturbation coefficient is 2.161118 1e-10 erg/cm^3 / (1e-12 erg/cm^3)

        }
        // Setup scalar tracker for flux tubes
        if ((NSCALARS > 0) ) {
          pscalars->s(0,k,j,i) = s1Profile(x1,x2,x3);
          pscalars->s(1,k,j,i) = s2Profile(x1,x2,x3);
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
          pcr->sigma_diff(0,k,j,i) = sigmaParl;
          pcr->sigma_diff(1,k,j,i) = sigmaPerp;
          pcr->sigma_diff(2,k,j,i) = sigmaPerp;
        }
      }
    }// end k,j,i

  }// End CR
  return;
}

//Source function with gravity and cooling
void CRSource(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &u_cr)
{
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
  #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        //GRAVITY
        Real xi = pmb->pcoord->x2v(j)/(nGrav*h);
        Real arg = (1/(nGrav*SQR(cosh(xi))) - SQR(tanh(xi)))*pow(cosh(xi),-1.0*nGrav);
        Real coeff = (beta* pres0)/ (sigmaPerp*SQR(h));
        u_cr(CRE,k,j,i) += arg*coeff*dt;

        }
      }
    }

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

        pcr->sigma_diff(0,k,j,i) = sigmaParl;
        pcr->sigma_diff(1,k,j,i) = sigmaPerp;
        pcr->sigma_diff(2,k,j,i) = sigmaPerp;

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

//----------------------------------------------------------------------------------------
void Mesh::UserWorkInLoop(void)
{
  MeshBlock *pmb = pblock;
  Real vmax = pmb->pcr->vmax;
  Real vcX = 0.0, vcY = 0.0, vcZ = 0.0;
  Real vFast  = 0.0;
  if (pmb->pcr->vmax_Edit) {
    while(pmb != nullptr){
      int kl = pmb->ks, ku = pmb->ke;
      int jl = pmb->js, ju = pmb->je;
      int il = pmb->is, iu = pmb->ie;
      for(int k=kl; k<=ku; ++k){
        for(int j=jl; j<=ju; ++j){
          for(int i=il; i<=iu; ++i){
            //vcX = (pmb->pcr->u_cr(CRE,k,j,i)-pmb->pcr->u_cr(CRE,k,j,i-1))/pmb->pcr->u_cr(CRE,k,j,i);
            //vcY = (pmb->pcr->u_cr(CRE,k,j,i)-pmb->pcr->u_cr(CRE,k,j-1,i))/pmb->pcr->u_cr(CRE,k,j,i);
            //vcZ = (pmb->pcr->u_cr(CRE,k,j,i)-pmb->pcr->u_cr(CRE,k-1,j,i))/pmb->pcr->u_cr(CRE,k,j,i);
            vcX = pmb->pfield->b.x1f(k,j,i)/std::sqrt(pmb->phydro->u(IDN,k,j,i));
            vcY = pmb->pfield->b.x2f(k,j,i)/std::sqrt(pmb->phydro->u(IDN,k,j,i));
            vcZ = pmb->pfield->b.x3f(k,j,i)/std::sqrt(pmb->phydro->u(IDN,k,j,i));
            vcX = pmb->phydro->u(IM1,k,j,i)/pmb->phydro->u(IDN,k,j,i);
            vcY = pmb->phydro->u(IM2,k,j,i)/pmb->phydro->u(IDN,k,j,i);
            vcZ = pmb->phydro->u(IM3,k,j,i)/pmb->phydro->u(IDN,k,j,i);
            vFast = std::max(vFast,std::sqrt( SQR(vcX) + SQR(vcY) + SQR(vcZ) ));
          }
        }
      }
      pmb = pmb->next;
    }
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE,&vFast,1,MPI_ATHENA_REAL,MPI_MAX,
                  MPI_COMM_WORLD);
#endif
    vFast = vFast;
    pmb = pblock;
    //std::cout<< "vMax  = " << vmax << std::endl;
    //std::cout<< "vFast = " << vFast << std::endl;

    if (vFast < vmax ) {
      while(pmb != nullptr){
        pmb->pcr->vmax = vFast;
        pmb = pmb->next;
      }
    }
  }

}
