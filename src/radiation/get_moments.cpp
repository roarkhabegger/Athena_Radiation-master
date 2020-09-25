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
//! \file get_moments.cpp
//  \brief calculate the moments of the radiation field
//======================================================================================


// Athena++ headers
#include "./radiation.hpp"
#include "./integrators/rad_integrators.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../defs.hpp"

//--------------------------------------------------------------------------------------
// \!fn void CalculateMoment()

// \brief function to create the radiation moments


// calculate the frequency integrated moments of the radiation field
// including the ghost zones
void Radiation::CalculateMoment(AthenaArray<Real> &ir_in)
{
  Real er, frx, fry, frz, prxx, pryy, przz, prxy, prxz, pryz;
  int n1z = pmy_block->block_size.nx1 + 2*(NGHOST);
  int n2z = pmy_block->block_size.nx2;
  int n3z = pmy_block->block_size.nx3;
  if(n2z > 1) n2z += (2*(NGHOST));
  if(n3z > 1) n3z += (2*(NGHOST));
  
  Real *weight = &(wmu(0));
  

  
  // reset the moment arrays to be zero
  // There are 13 3D arrays
  for(int n=0; n<13; ++n)
    for(int k=0; k<n3z; ++k)
      for(int j=0; j<n2z; ++j){
          Real *i_mom = &(rad_mom(n,k,j,0));
#pragma omp simd aligned(i_mom:ALI_LEN)
        for(int i=0; i<n1z; ++i){
          i_mom[i] = 0.0;
        }}

  
  for(int k=0; k<n3z; ++k){
    for(int j=0; j<n2z; ++j){
      for(int i=0; i<n1z; ++i){
        for(int ifr=0; ifr<nfreq; ++ifr){
          er=0.0; frx=0.0; fry=0.0; frz=0.0;
          prxx=0.0; pryy=0.0; przz=0.0; prxy=0.0;
          prxz=0.0; pryz=0.0;
          Real *intensity = &(ir_in(k,j,i,ifr*nang));
          Real *cosx = &(mu(0,k,j,i,0));
          Real *cosy = &(mu(1,k,j,i,0));
          Real *cosz = &(mu(2,k,j,i,0));
#pragma omp simd aligned(cosx,weight,intensity,cosy,cosz:ALI_LEN) reduction(+:er,frx,fry,frz,prxx,pryy,przz,prxy,prxz,pryz)
          for(int n=0; n<nang; ++n){
            Real irweight = weight[n] * intensity[n];
            er   += irweight;
            frx  += irweight * cosx[n];
            fry  += irweight * cosy[n];
            frz  += irweight * cosz[n];
            prxx += irweight * cosx[n] * cosx[n];
            pryy += irweight * cosy[n] * cosy[n];
            przz += irweight * cosz[n] * cosz[n];
            prxy += irweight * cosx[n] * cosy[n];
            prxz += irweight * cosx[n] * cosz[n];
            pryz += irweight * cosy[n] * cosz[n];
          }
          //multiply the frequency weight
          er *= wfreq(ifr);
          frx *= wfreq(ifr);
          fry *= wfreq(ifr);
          frz *= wfreq(ifr);
          prxx *= wfreq(ifr);
          pryy *= wfreq(ifr);
          przz *= wfreq(ifr);
          prxy *= wfreq(ifr);
          prxz *= wfreq(ifr);
          pryz *= wfreq(ifr);
          

          
          //assign the moments
          rad_mom(IER,k,j,i) += er;
          rad_mom(IFR1,k,j,i) += frx;
          rad_mom(IFR2,k,j,i) += fry;
          rad_mom(IFR3,k,j,i) += frz;
          rad_mom(IPR11,k,j,i) += prxx;
          rad_mom(IPR22,k,j,i) += pryy;
          rad_mom(IPR33,k,j,i) += przz;
          rad_mom(IPR12,k,j,i) += prxy;
          rad_mom(IPR13,k,j,i) += prxz;
          rad_mom(IPR23,k,j,i) += pryz;
          rad_mom(IPR21,k,j,i) += prxy;
          rad_mom(IPR31,k,j,i) += prxz;
          rad_mom(IPR32,k,j,i) += pryz;
          
          
        }// End frequency loop

      }
    }
    
  }
  
  
}


//--------------------------------------------------------------------------------------
// \!fn void CalculateComMoment()

// \brief Calculate the radiation moments in the co-moving frame
// Also load specific intensity for dump

void Radiation::CalculateComMoment()
{

  Hydro *phydro=pmy_block->phydro;
  Real invcrat = 1.0/crat;
  
  Real er, frx, fry, frz;
  int n1z = pmy_block->block_size.nx1 + 2*(NGHOST);
  int n2z = pmy_block->block_size.nx2;
  int n3z = pmy_block->block_size.nx3;
  if(n2z > 1) n2z += (2*(NGHOST));
  if(n3z > 1) n3z += (2*(NGHOST));
  
  Real *weight = &(wmu(0));
  
  AthenaArray<Real> &i_mom = rad_mom_cm;
  
  Real *ir_lab;

  
  // reset the moment arrays to be zero
  // There are 4 3D arrays
  for(int n=0; n<4; ++n)
    for(int k=0; k<n3z; ++k)
      for(int j=0; j<n2z; ++j)
#pragma omp simd
        for(int i=0; i<n1z; ++i){
          i_mom(n,k,j,i) = 0.0;
        }

  
  for(int k=0; k<n3z; ++k){
    for(int j=0; j<n2z; ++j){
      for(int i=0; i<n1z; ++i){

          Real *ir_lab = &(ir(k,j,i,0));
          Real *cosx = &(mu(0,k,j,i,0));
          Real *cosy = &(mu(1,k,j,i,0));
          Real *cosz = &(mu(2,k,j,i,0));
        
          Real vx = phydro->u(IM1,k,j,i)/phydro->u(IDN,k,j,i);
          Real vy = phydro->u(IM2,k,j,i)/phydro->u(IDN,k,j,i);
          Real vz = phydro->u(IM3,k,j,i)/phydro->u(IDN,k,j,i);
          Real vel = vx * vx + vy * vy + vz * vz;
        
          Real ratio = sqrt(vel) * invcrat;
         // Limit the velocity to be smaller than the speed of light
         if(ratio > vmax){
           Real factor = vmax/ratio;
           vx *= factor;
           vy *= factor;
           vz *= factor;
           
           vel *= (factor*factor);
         }

            // square of Lorentz factor
          Real lorzsq = 1.0/(1.0 - vel * invcrat * invcrat);
          Real lorz = sqrt(lorzsq);
 
        for(int ifr=0; ifr<nfreq; ++ifr){
          er=0.0; frx=0.0; fry=0.0; frz=0.0;
          Real numsum = 0.0;
#pragma omp simd aligned(cosx,cosy,cosz:ALI_LEN) reduction(+:numsum,er,frx,fry,frz)
          for(int n=0; n<nang; ++n){
            Real vdotn = vx * cosx[n] + vy * cosy[n] + vz * cosz[n];
            Real vnc = 1.0 - vdotn * invcrat;
            Real coef = vnc * vnc * vnc * vnc * lorzsq * lorzsq;
            Real wmu_cm = weight[n]/(lorzsq * vnc * vnc);
            numsum += wmu_cm;
            Real angcoef = lorz * invcrat * (1.0
                         - lorz * vdotn * invcrat/(1.0 + lorz));
            Real incoef = 1.0/(lorz * vnc);
            Real cosx_cm = (cosx[n] - vx * angcoef) * incoef;
            Real cosy_cm = (cosy[n] - vy * angcoef) * incoef;
            Real cosz_cm = (cosz[n] - vz * angcoef) * incoef;
            
            Real ir_weight = ir_lab[n] * coef * wmu_cm;
    
            er   += ir_weight;
            frx  += ir_weight * cosx_cm;
            fry  += ir_weight * cosy_cm;
            frz  += ir_weight * cosz_cm;
            
          }
          // normalize weight
          numsum = 1.0/numsum;
          //multiply the frequency weight
          er *= wfreq(ifr) * numsum;
          frx *= wfreq(ifr) * numsum;
          fry *= wfreq(ifr) * numsum;
          frz *= wfreq(ifr) * numsum;
          
          i_mom(IER,k,j,i) += er;
          i_mom(IFR1,k,j,i) += frx;
          i_mom(IFR2,k,j,i) += fry;
          i_mom(IFR3,k,j,i) += frz;
          
        }// End frequency loop
        
        // Now calculate frequency inetgrated opacity
        Real sum_sigma_s=0.0, sum_sigma_a = 0.0, sum_planck = 0.0;
        Real *sigmas=&(sigma_s(k,j,i,0));
        Real *sigmaa=&(sigma_a(k,j,i,0));
        Real *sigmap=&(sigma_planck(k,j,i,0));
        for(int ifr=0; ifr<nfreq; ++ifr){
          sum_sigma_s += sigmas[ifr] * wfreq(ifr);
          sum_sigma_a += sigmaa[ifr] * wfreq(ifr);
          sum_planck  += sigmap[ifr] * wfreq(ifr);
        }
        grey_sigma(OPAS,k,j,i) = sum_sigma_s;
        grey_sigma(OPAA,k,j,i) = sum_sigma_a;
        grey_sigma(OPAP,k,j,i) = sum_planck;

      }
    }
  }

  // Load the specific intensity for dump
//  if(ir_output > 0){
//    for(int n=0; n<ir_output; ++n){
//      for(int k=0; k<n3z; ++k){
//        for(int j=0; j<n2z; ++j){
//          for(int i=0; i<n1z; ++i){
//            dump_ir(n,k,j,i)=ir(k,j,i,ir_index(n));
//          }
//        }
//      }
//    }


//  }
  

  
}

