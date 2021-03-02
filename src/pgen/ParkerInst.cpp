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
Real H, gH;
Real g0;
Real alpha;

Real rhoProfile(Real x, Real Hg, Real H0, Real rho0);
Real presProfile(Real x, Real Hg, Real H0, Real rho0, Real pres0);
Real gravProfile(Real x, Real Hg, Real g0);

Real LogMeshSpacingX2(Real x, RegionSize rs);

void myGravity(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons);


void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh);

Real rhoProfile(Real x, Real Hg, Real H0, Real rho0){
  Real rho = rho0;
  Real prof = pow(cosh(x/Hg),-1.0*Hg/H0);
  rho *= prof;

  return rho;
}
Real presProfile(Real x, Real Hg, Real H0, Real rho0, Real pres0){
  Real pres = pres0;
  Real rho = rhoProfile(x,Hg,H0,rho0);
  pres *= rho/rho0;

  return pres;
}

Real gravProfile(Real x, Real Hg, Real g0){
  Real g = g0 * tanh(x/Hg);
  return g;
}

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
    EnrollUserBoundaryFunction(inner_x2, ProjectPressureInnerX2);
  }
  if (pin->GetString("mesh","ox2_bc")=="user"){
    EnrollUserBoundaryFunction(outer_x2, ProjectPressureOuterX2);
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
  pres0 = pin->GetReal("problem","Pres");
  
 
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x2 = pcoord->x2v(j);
        Real density = rhoProfile(x2,gH,H,dens0);
        Real pressure = presProfile(x2,gH,H,dens0,pres0);
        phydro->u(IDN,k,j,i) = density;
        phydro->u(IM1,k,j,i) = vx*density;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pressure/(gamma-1) + 0.5*density*pow(vx,2.0);
      }// end i
    }
  }    
    // initialize uniform interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    alpha = pin->GetOrAddReal("problem","alpha",0.0);
    Real B0   = pow(alpha*pres0,0.5);
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          Real x2 = pcoord->x2v(j);
          Real pressure= presProfile(x2,gH,H,dens0,pres0);

          Real B = B0*sqrt(pressure/pres0);

          pfield->b.x1f(k,j,i) = B;
          
          
        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          Real x2 = pcoord->x2v(j);
          Real x1 = pcoord->x1v(i);
          Real B = 0.2*B0;
          B     *= exp(-1.0*pow(x2-3.0,2.0)/(2.0*pow(1.0,2.0)));
          B     *= exp(-1.0*pow(x1-2.0,2.0)/(2.0*pow(0.5,2.0)));
          pfield->b.x2f(k,j,i) = B;
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
  return;
}

void myGravity(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &cons) 
{ for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js-NGHOST; j<=pmb->je+NGHOST; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real x2 = pmb->pcoord->x2v(j);
        Real gravity = gravProfile(x2,gH,g0);
        Real src = dt*prim(IDN,k,j,i)*gravity;
        cons(IM2,k,j,i) += src;
        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) += src*prim(IVY,k,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                            FaceField &b, Real time, Real dt,
                            int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        Real x2 = pmb->pcoord->x2v(jl-j);
        Real pressure = presProfile(x2,gH,H,dens0,pres0);
        Real gravity = gravProfile(x2,gH,g0);
        if (n==(IVY)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IVY,k,jl-j,i) = -prim(IVY,k,jl+j-1,i);  // reflect 2-velocity
          }
        } else if (n==(IPR)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IPR,k,jl-j,i) = prim(IPR,k,jl+j-1,i)
                                 + prim(IDN,k,jl+j-1,i)*gravity*(2*j-1)*pco->dx2f(jl-j)/(1+alpha);
          }
        } else {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(n,k,jl-j,i) = prim(n,k,jl+j-1,i);
          }
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          Real x2 = pmb->pcoord->x2v(jl-j);
          Real pressure= presProfile(x2,gH,H,dens0,pres0);
          Real gravity = gravProfile(x2,gH,g0);
          pressure = prim(IPR,k,jl+j-1,i)
                                 + prim(IDN,k,jl+j-1,i)*gravity*(2*j-1)*pco->dx2f(jl-j)/(1+alpha);
          Real B = sqrt(alpha*pressure);
          b.x1f(k,(jl-j),i) =  B;
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(jl-j),i) = -b.x2f(k,(jl+j  ),i);  // reflect 2-field
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(jl-j),i) =  b.x3f(k,(jl+j-1),i);
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
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        Real x2 = pmb->pcoord->x2v(ju+j);
        Real gravity = gravProfile(x2,gH,g0);
        Real pressure = presProfile(x2,gH,H,dens0,pres0);
        if (n==(IVY)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IVY,k,ju+j,i) = -prim(IVY,k,ju-j+1,i);  // reflect 2-velocity
          }
        } else if (n==(IPR)) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(IPR,k,ju+j,i) = prim(IPR,k,ju-j+1,i)
                                 - prim(IDN,k,ju-j+1,i)*gravity*(2*j-1)*pco->dx2f(ju+j)/(1+alpha);
          }
        } else {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            prim(n,k,ju+j,i) = prim(n,k,ju-j+1,i);
          }
        }
      }
    }
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          Real x2 = pmb->pcoord->x2v(ju+j);
          Real pressure= presProfile(x2,gH,H,dens0,pres0);
          Real gravity = gravProfile(x2,gH,g0);
          pressure = prim(IPR,k,ju-j+1,i)
                                 - prim(IDN,k,ju-j+1,i)*gravity*(2*j-1)*pco->dx2f(ju+j)/(1+alpha);
          Real B = sqrt(alpha*pressure);
          b.x1f(k,(ju+j),i) =  B;
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(ju+j),i) = -b.x2f(k,(ju-j  ),i);  // reflect 2-field
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(ju+j),i) =  b.x3f(k,(ju-j+1),i);
        }
      }
    }
  }
  return;
}


