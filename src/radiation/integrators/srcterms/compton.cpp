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
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../radiation.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../eos/eos.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../utils/utils.hpp"

// this class header
#include "../rad_integrators.hpp"

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::Absorption()
//  \brief 

// wmu_cm is the weight in the co-moving frame
// wmu_cm=wmu * 1/(1-vdotn/Crat)^2 / Lorz^2
// tran_coef is (1-vdotn/Crat)*Lorz
// rho is gas density
// tgas is gas temperature
// This function only update the absorption term in each cell

void RadIntegrator::Compton(AthenaArray<Real> &wmu_cm,
          AthenaArray<Real> &tran_coef, Real *sigma_s,
          Real dt, Real rho, Real &tgas, AthenaArray<Real> &ir_cm)
{

  Real& prat = pmy_rad->prat;
  Real ct = dt * pmy_rad->crat;
  Real redfactor=pmy_rad->reduced_c/pmy_rad->crat;

  int& nang=pmy_rad->nang;
  int& nfreq=pmy_rad->nfreq;
  Real gamma = pmy_rad->pmy_block->peos->GetGamma();
  Real telectron = 1.0/pmy_rad->telectron;
  
  
  // Polynomial coefficients for Compton
  Real coef[2];
  for (int i=0; i<2; ++i)
    coef[i] = 0.0;

  
  Real tgasnew = tgas;
  
  for(int ifr=0; ifr<nfreq; ++ifr){
  
    Real dtcsigma = ct * sigma_s[ifr];
    Real rdtcsigma = redfactor * dtcsigma;


    // Add Simple Compton scattering using the partically updated jr and ir

    if(rdtcsigma > TINY_NUMBER){
    
     // Calculate the sum \int gamma (1-vdotn/c) dw_0 4 dt csigma_s /T_e
      Real suma1 = 0.0, suma2 = 0.0, jr_cm=0.0, source=0.0;
      Real *irn = &(ir_cm(nang*ifr));
      Real *wmun = &(wmu_cm(0));
      Real *tcoef = &(tran_coef(0));
#pragma omp simd reduction(+:jr_cm,suma1) aligned(irn,wmun,tcoef:ALI_LEN)
      for(int n=0; n<nang; n++){
         jr_cm += irn[n] * wmun[n];
         suma1 += tcoef[n] * wmun[n] * 4.0 * rdtcsigma * telectron;
      }
      suma2 = 4.0 * prat * dtcsigma*(gamma-1.0)*telectron/rho;
      
      Real tr = sqrt(sqrt(jr_cm));
      Real trnew;
      
      if(fabs(tr - tgas) > 1.e-12){
        coef[1] = (1.0 + suma2* jr_cm)/(suma1 * jr_cm);
        coef[0] = -(1.0+suma2*jr_cm)/suma1-tgas;
        
        int flag = FouthPolyRoot(coef[1], coef[0], trnew);
        if(flag == -1){
          trnew = tr;
          tgasnew = tgas;
          source = 0.0;
        }else{
        
          Real jrnew = trnew * trnew * trnew * trnew;
    
          tgasnew = (jrnew - jr_cm)/(suma1*jr_cm) + trnew;
          source = rdtcsigma * 4.0 * jr_cm * telectron * (tgasnew - trnew);
       }
      }
      // Update the co-moving frame specific intensity

#pragma omp simd aligned(irn,tcoef:ALI_LEN)
      for(int n=0; n<nang; n++){
        irn[n] += source * tcoef[n];
      }
      
    }// End Compton
    
  }// End Frequency
  
  // Update gas temperature
  tgas = tgasnew;

  return;
}



