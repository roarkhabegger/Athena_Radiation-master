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
Real d0, d1, p0, p1, vx0, vx1;
Real w;
int cooling;
Real T0, T1, T2, A, B, C;


//Profile functions
Real densProfile(Real x1, Real x2);
Real presProfile(Real x1, Real x2);

//For logarithmic height spacing
Real LogMeshSpacingX2(Real x, RegionSize rs);

void CoolingFunction(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &cons);


Real densProfile(Real x){
  Real dens = (d0-d1)*0.5*(1- tanh(x/w)) + d1;
  return dens;
}

Real presProfile(Real x){
  Real pres = (p0-p1)*0.5*(1- tanh(x/w)) + p1;
  return pres;
}

Real velProfile(Real x){
  Real vel = (vx0-vx1)*0.5*(1- tanh(x/w)) + vx1;
  return vel;
}


Real LogMeshSpacingX2(Real x, RegionSize rs){
  Real xf, xrat;
  xrat = pow(rs.x2max/rs.x2min,1.0/((Real) rs.nx2));
  xf = rs.x2min*pow(xrat,x*rs.nx2);
  return xf;
}

//X1 boundaries to match HSE
void ProfilesOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void ProfilesInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
      Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);



void Mesh::InitUserMeshData(ParameterInput *pin){

  Real x2rat = pin->GetOrAddReal("mesh","x2rat",0.0);
  if (x2rat< 0.0) {
    EnrollUserMeshGenerator(X2DIR,LogMeshSpacingX2);
  }

  cooling = pin->GetOrAddInteger("problem","cooling",0);

  if (cooling == 1)   EnrollUserExplicitSourceFunction(CoolingFunction);
  if (pin->GetString("mesh","ix1_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x1, ProfilesInnerX1);
  }
  if (pin->GetString("mesh","ox1_bc")=="user"){
    EnrollUserBoundaryFunction(outer_x1, ProfilesOuterX1);
  }

}


void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Real gamma = peos->GetGamma();

  d0 = pin->GetReal("problem","LeftDens");
  p0 = pin->GetReal("problem","LeftPres");
  vx0 = pin->GetReal("problem","LeftXVel");
  d1 = pin->GetReal("problem","RightDens");
  p1 = pin->GetReal("problem","RightPres");
  vx1 = pin->GetReal("problem","RightXVel");
  w = pin->GetReal("problem","LayerWidth");

  if (cooling == 1){
    //Computational scale for T_s = 12441
    T0 = pin->GetOrAddReal("problem","CoolingT0",9.773);
    T1 = pin->GetOrAddReal("problem","CoolingT1",0.1238);
    T2 = pin->GetOrAddReal("problem","CoolingT2",0.007594);

    //Computational scale for P_s*n_s/t_s
    A = pin->GetOrAddReal("problem","CoolingA",651790.0);
    B = pin->GetOrAddReal("problem","CoolingB",0.705);
    C = pin->GetOrAddReal("problem","CoolingC",3.0);
  }
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real density = densProfile(x1);
        Real pressure = presProfile(x1);
        Real velocity = velProfile(x1);

        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = velocity*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pressure/(gamma-1) + 0.5*density*pow(velocity,2.0);

      }// end i
    }
  }

  return;
}

void CoolingFunction(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &cons)
{
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real temp = prim(IPR,k,j,i)/prim(IDN,k,j,i);
        Real n = prim(IDN,k,j,i);
        cons(IEN,k,j,i) -= pow(n,2.0)*(A*exp(-1*T0/(temp-T1))+B*exp(-1*T2/temp))-n*C;
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
        Real x1 = pco->x1v(is-i);
        Real density = densProfile(x1);
        Real pressure = presProfile(x1);
        Real velocity = velProfile(x1);
        prim(IDN,k,j,is-i) = density+(prim(IDN,k,j,is) - densProfile(edge));
        prim(IVX,k,j,is-i) = velocity+(prim(IVX,k,j,is) - velProfile(edge));
        prim(IVY,k,j,is-i) = prim(IVY,k,j,is);
        prim(IVZ,k,j,is-i) = prim(IVZ,k,j,is);
        prim(IEN,k,j,is-i) = pressure+(prim(IEN,k,j,is) - presProfile(edge));

      }
    }
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
        Real x1 = pco->x1v(ie+i);
        Real density = densProfile(x1);
        Real pressure = presProfile(x1);
        Real velocity = velProfile(x1);
        prim(IDN,k,j,ie+i) = density+(prim(IDN,k,j,ie) - densProfile(edge));
        prim(IVX,k,j,ie+i) = velocity+(prim(IVX,k,j,ie) - velProfile(edge));
        prim(IVY,k,j,ie+i) = prim(IVY,k,j,ie);
        prim(IVZ,k,j,ie+i) = prim(IVZ,k,j,ie);
        prim(IEN,k,j,ie+i) = pressure+(prim(IEN,k,j,ie) - presProfile(edge));

      }
    }
  }


  return;
}
