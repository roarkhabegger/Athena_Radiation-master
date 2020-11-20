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

Real dens0;
Real H, gH;
Real g0;

Real LogMeshSpacingX2(Real x, RegionSize rs);

void myGravity(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons);

void FixLowerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh);

Real LogMeshSpacingX2(Real x, RegionSize rs){
  Real xf, xrat;
  xrat = pow(rs.x2max/rs.x2min,1.0/((Real) rs.nx2));
  xf = rs.x2min*pow(xrat,x*rs.nx2); 
  return xf;
}

void Mesh::InitUserMeshData(ParameterInput *pin){

  Real x2rat = pin->GetOrAddReal("mesh","x2rat",0.0);
  if (x2rat< 0.0) {
    EnrollUserMeshGenerator(X2DIR,LogMeshSpacingX2);
  }
  if (pin->GetString("mesh","ix2_bc")=="user"){
    EnrollUserBoundaryFunction(inner_x2, FixLowerX2);
  }
  EnrollUserExplicitSourceFunction(myGravity);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
  Real gamma = peos->GetGamma();

  Real vx=pin->GetReal("problem","xVel");
  H = pin->GetReal("problem","ScaleHeight");
  gH = pin->GetReal("problem","gScaleHeight");
  
  
  dens0 = pin->GetReal("problem","Dens");
  g0 = pin->GetReal("problem","grav");
 
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x2 = pcoord->x2v(j);
         
        Real prof = pow(cosh(x2/gH),-1.0*gH/H);
        Real density = dens0*prof;  
        Real vel = vx;

        Real gProf = g0*tanh(x2/(gH));

        phydro->u(IDN,k,j,i) = density;

        phydro->u(IM1,k,j,i) = vx*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        
      }// end i
    }
  }
    
  return;
}

void myGravity(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &cons) 
{ for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x2 = pmb->pcoord->x2v(j);
        Real gProf = g0*tanh(x2/(gH));
        Real src = dt*prim(IDN,k,j,i)*gProf;
        cons(IM2,k,j,i) += src;
      }
    }
  }
  return;
}

void FixLowerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, 
     int ks, int ke, int ngh){
  // copy hydro variables into ghost zones, reflecting v1
  //std::cout << "j=" << 0 << " flux=" << pmb->phydro->flux[X2DIR](IDN,0,3,30) << std::endl;
  for (int k=ks; k<=ke; ++k) {
    for (int i=is; i<=ie; ++i) {
      for (int j=1; j<=ngh; ++j) {
        Real x2 = pco->x2v(js-j);
        Real prof = pow(cosh(x2/gH),-1.0*gH/H);
        //Real prof = exp(-1.0*x2/H);
        //std::cout << "j="<< js-j << " dens0=" << dens0 <<std::endl;
        prim(IDN,k,js-j,i) = dens0*prof;
        prim(IVX,k,js-j,i) = 0.0; // reflect 1-velocity
        prim(IVY,k,js-j,i) = 0.0;
        prim(IVZ,k,js-j,i) = 0.0;
      }
    }
  }
}


