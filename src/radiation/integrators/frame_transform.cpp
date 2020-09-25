//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file absorption.cpp
//  \brief  Add absorption source terms
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../radiation.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../utils/utils.hpp"

// this class header
#include "./rad_integrators.hpp"

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::LabToCom(const Real vx, const Real vy, const Real vz,
//                          AthenaArray<Real> &ir, AthenaArray<Real> &ir_cm)
//  \brief Transform specific intensity from lab frame to co-moving frame
// with flow velocity vx, vy, vz



void RadIntegrator::LabToCom(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          Real *ir_lab, AthenaArray<Real> &ir_cm)
{

  Real& prat = pmy_rad->prat;
  Real invcrat = 1.0/pmy_rad->crat;
  int& nang=pmy_rad->nang;
  int& nfreq=pmy_rad->nfreq;
  
  
  
  // square of Lorentz factor
  Real lorzsq = 1.0/(1.0 - (vx * vx + vy * vy + vz * vz) * invcrat * invcrat);
  


  for(int ifr=0; ifr<nfreq; ++ifr){
#pragma omp simd
    for(int n=0; n<nang; n++){
       Real vnc = vx * mux[n] + vy * muy[n] + vz * muz[n];
      
       vnc = 1.0 - vnc * invcrat;
       Real coef = vnc * vnc * vnc * vnc * lorzsq * lorzsq;
      
       ir_cm(n+ifr*nang) = ir_lab[n+ifr*nang] * coef;
    }
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::ComToLab(const Real vx, const Real vy, const Real vz,
//                          AthenaArray<Real> &ir, AthenaArray<Real> &ir_cm)
//  \brief Transform specific intensity from co-moving frame to lab frame
// with flow velocity vx, vy, vz



void RadIntegrator::ComToLab(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          AthenaArray<Real> &ir_cm, Real *ir_lab)
{

  Real& prat = pmy_rad->prat;
  Real invcrat = 1.0/pmy_rad->crat;
  int& nang=pmy_rad->nang;
  int& nfreq=pmy_rad->nfreq;
  
  
  
  // square of Lorentz factor
  Real lorzsq = 1.0/(1.0 - (vx * vx + vy * vy + vz * vz) * invcrat * invcrat);
  


  for(int ifr=0; ifr<nfreq; ++ifr){
#pragma omp simd
    for(int n=0; n<nang; n++){
       Real vnc = vx * mux[n] + vy * muy[n] + vz * muz[n];
      
       vnc = 1.0 - vnc * invcrat;
       Real coef = vnc * vnc * vnc * vnc * lorzsq * lorzsq;
      
       ir_lab[n+ifr*nang] = ir_cm(n+ifr*nang) / coef;
    }
  }

  return;
}


//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::ComAngle(const Real vx, const Real vy, const Real vz,
//                  Real mux, Real muy, Real muz, Real mux0, Real muy0, Real muz0)
//  \brief Transform angles from lab frame to co-moving frame


void RadIntegrator::ComAngle(const Real vx, const Real vy, const Real vz,
                  Real mux, Real muy, Real muz, Real *mux0, Real *muy0, Real *muz0)
{

  Real invcrat = 1.0/pmy_rad->crat;
  
  
  // square of Lorentz factor
  Real lorz = 1.0/(1.0 - (vx * vx + vy * vy + vz * vz) * invcrat * invcrat);
  lorz = sqrt(lorz);
  Real vdotn = vx * mux + vy * muy + vz * muz;
  
  Real angcoef = lorz * invcrat * (1.0 - lorz * vdotn * invcrat/(1.0 + lorz));
  Real incoef = 1.0/(lorz*(1.0-vdotn * invcrat));
  
  (*mux0) = (mux - vx * angcoef) * incoef;
  (*muy0) = (muy - vy * angcoef) * incoef;
  (*muz0) = (muz - vz * angcoef) * incoef;
  

  return;
}
