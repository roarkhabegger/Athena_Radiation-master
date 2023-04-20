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
#include <random> 
#include <vector>
#include <chrono>

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
Real angle;    // angle off of x axis, in the y direction
Real beta;     // Ratio of cosmic ray pressure to gas pressure

Real h;

//cr Diffusion variables and function
Real sigmaParl, sigmaPerp;

//Floors for Diode boundary conds
Real dfloor, pfloor; // Floor values for density and rpessure

const Real float_min{std::numeric_limits<float>::min()};
//Useful minimum floating point number

//Cooling and perturbation scalar
int cooling; //Boolean - if cooling==1 do Inoue 2006 2 phase gas cooling profile
std::vector<double> X1Inj = {};
std::vector<double> X2Inj = {};
std::vector<double> X3Inj = {};
int NInjs = 0;
int TotalInjs = 0;
double lastInjT = 0.0;
double SNRate = 0.0;
double Esn = 0.0;
// std::vector<float> injE  = { };
// std::vector<float> injTr = { };
// std::vector<float> injTf = { };
// std::vector<float> injT  = { };
// std::vector<float> injX1 = { };
// std::vector<float> injX2 = { };
// std::vector<float> injX3 = { };


//Profile functions
Real densProfile(Real x1, Real x2, Real x3);
Real presProfile(Real x1, Real x2, Real x3);
Real fcProfile(Real x1, Real x2, Real x3);


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

Real fcProfile(Real x1, Real x2, Real x3)
{
  Real  fcz = pow(cosh(x2/(nGrav*h)),-1.0*nGrav)*tanh(x2/(nGrav*h));
  fcz *= beta*pres0/(h*sigmaPerp) ;
  return fcz;
}

// Set local (each cell independently) opacity funciton
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  if(CR_ENABLED){
   pcr->EnrollOpacityFunction(Diffusion);
   pcr->EnrollUserCRSource(CRSource);
  }
}



//----------------------------------------------------------------------------------------
void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank==0) {
    std::cout << "Total Number of Injections = " << TotalInjs << std::endl;
  }
}

//----------------------------------------------------------------------------------------
void Mesh::UserWorkInLoop(void)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  X1Inj.clear();
  X2Inj.clear();
  X3Inj.clear();
  NInjs = 0;
  if (rank == 0) {
    Real x1d = (mesh_size.x1max - mesh_size.x1min)/float(mesh_size.nx1);
    Real x2d = (mesh_size.x2max - mesh_size.x2min)/float(mesh_size.nx2);;
    Real x3d = (mesh_size.x3max - mesh_size.x3min)/float(mesh_size.nx3);;
    //std::cout << mesh_size.x1min << "," << mesh_size.x1max << std::endl;

    std::exponential_distribution<double> distDt(SNRate);
    std::uniform_real_distribution<double> distx1(mesh_size.x1min,mesh_size.x1max);
    std::uniform_real_distribution<double> distx2(-0.1,0.1);
    std::uniform_real_distribution<double> distx3(mesh_size.x3min,mesh_size.x3max);

    unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed1);

    while ((lastInjT < time+dt) && (SNRate> 0.0)){
      X1Inj.insert(X1Inj.end(), round((distx1(gen)-mesh_size.x1min)/x1d)*x1d + mesh_size.x1min + 0.5*x1d);
      X2Inj.insert(X2Inj.end(), round((distx2(gen)-mesh_size.x2min)/x2d)*x2d + mesh_size.x2min + 0.5*x2d);
      X3Inj.insert(X3Inj.end(), round((distx3(gen)-mesh_size.x3min)/x3d)*x3d + mesh_size.x3min + 0.5*x3d);
      lastInjT += distDt(gen);
      NInjs +=1;
    }
    if  (SNRate <= 0.0) {
        lastInjT = time + dt;
    }
    //std::cout << time << " " << NInjs <<std::endl;
  } 
  
  MPI_Bcast(&lastInjT,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&NInjs,1,MPI_INT,0,MPI_COMM_WORLD);

  if ((NInjs > 0) && (rank != 0)){
    X1Inj.insert(X1Inj.end(),NInjs,-15.0);
    X2Inj.insert(X2Inj.end(),NInjs,-15.0);
    X3Inj.insert(X3Inj.end(),NInjs,-15.0);
  }

  MPI_Bcast(&X1Inj[0],NInjs,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&X2Inj[0],NInjs,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&X3Inj[0],NInjs,MPI_DOUBLE,0,MPI_COMM_WORLD);
  TotalInjs += NInjs;
}


//Set up initial MESH data
void Mesh::InitUserMeshData(ParameterInput *pin) {
  myGamma = pin->GetReal("hydro","gamma");

  // Load variables
  vx=pin->GetReal("problem","xVel");
  nGrav = pin->GetReal("problem","GravNumScaleHeight");
  beta = pin->GetOrAddReal("problem","beta",0.0);
  alpha = pin->GetOrAddReal("problem","alpha",0.0);
  angle = pin->GetOrAddReal("problem","angle",0.0);
  pres0 = pin->GetReal("problem","Pres");

  dens0 = pin->GetReal("problem","Dens");

  g0 =pin->GetReal("problem","Grav");
  h  =pin->GetReal("problem","ScaleH");
  
  dfloor = pin->GetOrAddReal("hydro", "dfloor", std::sqrt(1024*float_min)) ;
  pfloor = pin->GetOrAddReal("hydro", "pfloor", std::sqrt(1024*float_min)) ;
  
  if(CR_ENABLED){
    //Load CR Variables 
    sigmaPerp = pin->GetReal("cr","sigmaPerp");
    sigmaParl = pin->GetReal("cr","sigmaParl");
    SNRate = pin->GetReal("problem","SNRate");
    Esn = pin->GetOrAddReal("problem","Esn",2.12842045e-7);
    //.        Esn, Rad, tRise, tFall, t0, x10, x20, x30 
  }  
  // MHD boundary conditions
  if (pin->GetString("mesh","ix2_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x2, DiodeInnerX2);
  }
  if (pin->GetString("mesh","ox2_bc")=="user"){
    EnrollUserBoundaryFunction(outer_x2, DiodeOuterX2);
  }

  // Source Functions
  cooling = pin->GetOrAddInteger("problem","cooling",0);
  EnrollUserExplicitSourceFunction(mySource);
  

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
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real density = densProfile(x1,x2,x3);
        Real pressure = presProfile(x1,x2,x3);
        
        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = vx*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pressure/(myGamma-1) + 0.5*density*SQR(vx) ;//+ pot;

        if(CR_ENABLED){
          // get CR parameters
          Real crp = beta*pressure;
          pcr->u_cr(CRE,k,j,i) = 3.0*crp;
          pcr->u_cr(CRF1,k,j,i) = 0.0;
          pcr->u_cr(CRF2,k,j,i) = fcProfile(x1,x2,x3);//-1.0*dPcdz/sigmaParl;
          pcr->u_cr(CRF3,k,j,i) = 0.0;
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

          Real B = B0*sqrt(pressure)*cos(angle);

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
          Real x2 = pcoord->x2v(j);
          Real x1 = pcoord->x1v(i);
          Real x3 = pcoord->x3f(k);
          Real pressure= presProfile(x1,x2,x3);

          Real B = B0*sqrt(pressure)*sin(angle);
          pfield->b.x3f(k,j,i) = B;
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

void CRSource(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &u_cr)
{   // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // std::cout << time << " " << NInjs <<std::endl;
    // if (NInjs == 1) {std::cout << time << " Node " << rank << " has X1=" << X1Inj.at(0) << std::endl;}
  
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
  #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        //HSE forcing
        Real xi = pmb->pcoord->x2v(j)/(nGrav*h);
        Real arg = (1/(nGrav*SQR(cosh(xi))) - SQR(tanh(xi)))*pow(cosh(xi),-1.0*nGrav);
        Real coeff = (beta* pres0)/ (sigmaPerp*SQR(h));
        u_cr(CRE,k,j,i) += arg*coeff*dt;
        for (int m = 0 ; m < NInjs; ++m) {
          Real x1fl = pmb->pcoord->x1f(i);
          Real x2fl = pmb->pcoord->x2f(j);
          Real x3fl = pmb->pcoord->x3f(k);
          Real x1fr = pmb->pcoord->x1f(i+1);
          Real x2fr = pmb->pcoord->x2f(j+1);
          Real x3fr = pmb->pcoord->x3f(k+1);
          Real vol = (x1fr-x1fl)*(x2fr-x2fl)*(x3fr-x3fl);
          Real x10   = X1Inj.at(m);
          Real x20   = X2Inj.at(m);
          Real x30   = X3Inj.at(m);
          if ((x10<x1fr) && (x10>x1fl) && (x20<x2fr) && (x20>x2fl) && (x30<x3fr) && (x30>x3fl) ) {
          // std::cout << "  injection " << m << " at " <<  x10 << "," << x20 << "," << x30 <<std::endl;
            u_cr(CRE,k,j,i) += Esn/vol;
          }
        }
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
  Real T0  =  2.59008249e-6 ;
  Real T1  =  1.72672166e-5 ;
  Real T2  =  6.90688664e-5 ;
  Real T3  =  8.63360831e-4 ;
  Real T4  =  3.45344332e-1 ;
  Real A1  =  4.41644786e-11 ;
  Real A2  =  1.97575854e-9 ;
  Real A3  =  9.12495752e-15 ;
  Real A4  =  3.51263503e3 ;
  Real A5  =  6.35767207e-6 ;

  Real Heat = 0.0;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        //GRAVITY
        Real x2 = pmb->pcoord->x2v(j);
        // Real x1 = pmb->pcoord->x1v(i);
        // Real x3 = pmb->pcoord->x3v(k);
        Real gravity = -1*g0*tanh(x2/(nGrav*h));
        Real src = dt*prim(IDN,k,j,i)*gravity;

        cons(IM2,k,j,i) += src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);

        //COOLING
        if (cooling == 1) {
          Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i);
          Real n = prim(IDN,k,j,i);
          Real Lamb = 0.0;
          if ((temp >= T0) && (temp < T1)) { 
            Lamb = A1 * pow(temp,2.0);
          } else if ((temp >= T1) && (temp < T2)){
            Lamb = A2 * pow(temp,1.5);
          } else if ((temp >= T2) && (temp < T3)){
            Lamb = A3 * pow(temp,2.867);
          } else if ((temp >= T3) && (temp < T4)){
            Lamb = A4 * pow(temp,-0.65);
          } else if ((temp >= T4) ){
            Lamb = A5 * pow(temp,0.5);
          }
          // Lamb = 1e-2;
          cons(IEN,k,j,i) -= dt*(pow(n,2.0)*Lamb-n*Heat);
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
