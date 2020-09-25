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
//! \file rad_source.cpp
//  \brief Add radiation source terms to both radiation and gas
//======================================================================================




#include <algorithm>
#include <string>
#include <vector>

//Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../radiation.hpp"
#include "../../coordinates/coordinates.hpp" //
#include "../../hydro/hydro.hpp"
#include "../../field/field.hpp"
#include "../../eos/eos.hpp"

// class header
#include "rad_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void RadIntegrator::AddSourceTerms(MeshBlock *pmb, const Real dt, AthenaArray<Real> &u,
        AthenaArray<Real> &w, AthenaArray<Real> &bcc, AthenaArray<Real> &ir)
{
  Radiation *prad=pmb->prad;
  Coordinates *pco = pmb->pcoord;
  
  Real& prat = prad->prat;
  Real invcrat = 1.0/prad->crat;
  Real invredc = 1.0/prad->reduced_c;
  Real invredfactor = invredc/invcrat;

  Real gm1 = pmb->peos->GetGamma() - 1.0;

  Real rho_floor = pmb->peos->GetDensityFloor();
  
  
  Real *sigma_at, *sigma_aer, *sigma_s, *sigma_p;
  Real *lab_ir;
  
  int &nang =prad->nang;
  int &nfreq=prad->nfreq;
  
  
  // Get the temporary array
  AthenaArray<Real> &wmu_cm = wmu_cm_;
  AthenaArray<Real> &tran_coef = tran_coef_;
  AthenaArray<Real> &ir_cm = ir_cm_;
  AthenaArray<Real> &cm_to_lab = cm_to_lab_;


  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
 
  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){

        // for implicit update, using the quantities from the partially 
        // updated u, not from w 

         Real rho = u(IDN,k,j,i);
         rho = std::max(rho,rho_floor);
         Real vx = u(IM1,k,j,i)/rho;
         Real vy = u(IM2,k,j,i)/rho;
         Real vz = u(IM3,k,j,i)/rho;
         Real pb = 0.0;
         if(MAGNETIC_FIELDS_ENABLED)
           pb = 0.5*(SQR(bcc(IB1,k,j,i))+SQR(bcc(IB2,k,j,i))
                +SQR(bcc(IB3,k,j,i)));

         Real vel = vx * vx + vy * vy + vz * vz;
         Real tgas = u(IEN,k,j,i) - pb - 0.5*rho*vel;
         tgas = gm1*tgas/rho;
         tgas = std::max(tgas,prad->t_floor_);

         // calculate radiation energy density
         Real er = 0.0;
         for(int ifr=0; ifr<nfreq; ++ifr){ 
           Real *irn = &(ir(k,j,i,ifr*nang));
           Real *weight = &(prad->wmu(0));
           Real er_freq = 0.0;
#pragma omp simd reduction(+:er_freq)
           for(int n=0; n<nang; ++n){
             er_freq += weight[n] * irn[n];
           }
           er_freq *= prad->wfreq(ifr);
           er += er_freq;
         }       

         // Do not use the velocity directly in strongly radiation pressure
         // dominated regime
         // use the predicted velocity based on moment equation
         if(prat * er * invcrat * invcrat > rho){

            PredictVel(ir,k,j,i, 0.5 * dt, rho, &vx, &vy, &vz);
            vel = vx * vx + vy * vy + vz * vz;

         }
        
         Real ratio = sqrt(vel) * invcrat;
         // Limit the velocity to be smaller than the speed of light
         if(ratio > prad->vmax){
           Real factor = prad->vmax/ratio;
           vx *= factor;
           vy *= factor;
           vz *= factor;
           
           vel *= (factor*factor);
         }


        
        
         Real lorzsq = 1.0/(1.0 - vel  * invcrat * invcrat);
        

        
         sigma_at=&(prad->sigma_a(k,j,i,0));
         sigma_aer=&(prad->sigma_ae(k,j,i,0));
         sigma_s=&(prad->sigma_s(k,j,i,0));
         sigma_p=&(prad->sigma_planck(k,j,i,0));

         // Prepare the transformation coefficients
         Real numsum = 0.0;

#pragma omp simd reduction(+:numsum)
         for(int n=0; n<nang; ++n){
            Real vdotn = vx * prad->mu(0,k,j,i,n) + vy * prad->mu(1,k,j,i,n)
                        + vz * prad->mu(2,k,j,i,n);
            Real vnc = 1.0 - vdotn * invcrat;
            tran_coef(n) = sqrt(lorzsq) * vnc;
            wmu_cm(n) = prad->wmu(n)/(tran_coef(n) * tran_coef(n));
            numsum += wmu_cm(n);
            cm_to_lab(n) = tran_coef(n)*tran_coef(n)*tran_coef(n)*tran_coef(n);
           
         }
           // Normalize weight in co-moving frame to make sure the sum is one
         numsum = 1.0/numsum;
#pragma omp simd
         for(int n=0; n<nang; ++n){
            wmu_cm(n) *= numsum;
         }
        
         lab_ir=&(ir(k,j,i,0));
        
         for(int ifr=0; ifr<nfreq; ++ifr){
#pragma omp simd
           for(int n=0; n<nang; ++n){
             ir_cm(n+ifr*nang) = lab_ir[n+ifr*nang] * cm_to_lab(n);
           }
         }// End frequency
        
         // Add absorption and scattering opacity source
         AbsorptionScattering(wmu_cm,tran_coef, sigma_at, sigma_p, sigma_aer,
                              sigma_s, dt, rho, tgas, ir_cm);
        
         // Add compton scattering
         if(compton_flag_ > 0)
           Compton(wmu_cm,tran_coef, sigma_s, dt, rho, tgas, ir_cm);
       
        //After all the source terms are applied
        // Calculate energy and momentum source terms
        
        // first, calculate Er and Fr in lab frame before the step
        Real er0 = 0.0;
        Real frx0 = 0.0;
        Real fry0 = 0.0;
        Real frz0 = 0.0;
      if(prad->set_source_flag == 0){
        
        for(int ifr=0; ifr<nfreq; ++ifr){
#pragma omp simd reduction (+:er0,frx0,fry0,frz0)
           for(int n=0; n<nang; ++n){
               Real ir_weight = lab_ir[n+ifr*prad->nang] * prad->wmu(n);
               er0 += ir_weight;
               frx0 += ir_weight * prad->mu(0,k,j,i,n);
               fry0 += ir_weight * prad->mu(1,k,j,i,n);
               frz0 += ir_weight * prad->mu(2,k,j,i,n);
           }
           er0 *= prad->wfreq(ifr);
           frx0 *= prad->wfreq(ifr);
           fry0 *= prad->wfreq(ifr);
           frz0 *= prad->wfreq(ifr);
        }
      }
       
       // now update the lab frame intensity
        for(int ifr=0; ifr<nfreq; ++ifr){
          for(int n=0; n<nang; ++n){
              lab_ir[n+ifr*nang] = std::max(ir_cm(n+ifr*nang)/cm_to_lab(n), TINY_NUMBER);
          }
        }

       // now calculate the new moments
               // first, calculate Er and Fr in lab frame before the step
      if(prad->set_source_flag == 0){
        Real er = 0.0;
        Real frx = 0.0;
        Real fry = 0.0;
        Real frz = 0.0;
        
        for(int ifr=0; ifr<nfreq; ++ifr){
#pragma omp simd reduction (+:er,frx,fry,frz)
           for(int n=0; n<nang; ++n){
               Real ir_weight = lab_ir[n+ifr*prad->nang] * prad->wmu(n);
               er += ir_weight;
               frx += ir_weight * prad->mu(0,k,j,i,n);
               fry += ir_weight * prad->mu(1,k,j,i,n);
               frz += ir_weight * prad->mu(2,k,j,i,n);
           }
           er *= prad->wfreq(ifr);
           frx *= prad->wfreq(ifr);
           fry *= prad->wfreq(ifr);
           frz *= prad->wfreq(ifr);
        }
        
        
        // Now apply the radiation source terms to gas with energy and
        // momentum conservation
        u(IM1,k,j,i) += (-prat*(frx- frx0) * invredc);
        u(IM2,k,j,i) += (-prat*(fry- fry0) * invredc);
        u(IM3,k,j,i) += (-prat*(frz- frz0) * invredc);
        u(IEN,k,j,i) += (-prat*(er - er0 ) * invredfactor);
        
        //limit the velocity by speed of light
        vx = u(IM1,k,j,i)/u(IDN,k,j,i);
        vy = u(IM2,k,j,i)/u(IDN,k,j,i);
        vz = u(IM3,k,j,i)/u(IDN,k,j,i);
        vel = vx*vx+vy*vy+vz*vz;
        vel = sqrt(vel);
        ratio = vel * invcrat;
        if(ratio > prad->vmax){
          Real factor = prad->vmax/ratio;
          u(IM1,k,j,i) *= factor;
          u(IM2,k,j,i) *= factor;
          u(IM3,k,j,i) *= factor;
        
        }
        
        // check internal energy
//        Real ekin=0.5 * (u(IM1,k,j,i) * u(IM1,k,j,i)
//                       + u(IM2,k,j,i) * u(IM2,k,j,i)
//                       + u(IM3,k,j,i) * u(IM3,k,j,i))/u(IDN,k,j,i);
//        Real pb=0.0;
//        if(MAGNETIC_FIELDS_ENABLED){
//          if(step==1){
//            const Real& bcc1 = pmb->pfield->bcc1(IB1,k,j,i);
//            const Real& bcc2 = pmb->pfield->bcc1(IB2,k,j,i);
//            const Real& bcc3 = pmb->pfield->bcc1(IB3,k,j,i);
//            pb=0.5*(SQR(bcc1) + SQR(bcc2) + SQR(bcc3));
//          }else{
//            const Real& bcc1 = pmb->pfield->bcc(IB1,k,j,i);
//            const Real& bcc2 = pmb->pfield->bcc(IB2,k,j,i);
//            const Real& bcc3 = pmb->pfield->bcc(IB3,k,j,i);
//            pb=0.5*(SQR(bcc1) + SQR(bcc2) + SQR(bcc3));
//          }
//        }
//        Real eint=u(IEN,k,j,i) - pb - ekin;
//        if(eint < 0.0){
//          Real gm1 = pmb->phydro->peos->GetGamma() - 1.0;
//          eint = u(IDN,k,j,i) * tgas /gm1;
//          u(IEN,k,j,i) = eint  + pb + ekin;
//        }
        
       }//source flag
        
      }// end i
    }// end j
  }// end k
    


  
}


void RadIntegrator::PredictVel(AthenaArray<Real> &ir, int k, int j, int i, 
      Real dt, Real rho, Real *vx, Real *vy, Real *vz)
{
    Radiation *prad = pmy_rad;
  
    Real &prat = prad->prat;
    Real invcrat = 1.0/prad->crat;
    Real ct = dt * prad->reduced_c;
    int& nang =prad->nang;
    int& nfreq=prad->nfreq;
    // first, calculate the moments
    Real er =0.0, fr1=0.0, fr2=0.0, fr3=0.0,
         pr11=0.0,pr12=0.0,pr13=0.0,pr22=0.0,
         pr23=0.0,pr33=0.0;
    Real *weight = &(prad->wmu(0));
    for(int ifr=0; ifr<nfreq; ++ifr){
      Real er_f = 0.0,fr1_f=0.0,fr2_f=0.0,fr3_f=0.0,
           pr11_f=0.0,pr12_f=0.0,pr13_f=0.0,pr22_f=0.0,
           pr23_f=0.0,pr33_f=0.0;

      Real *irn = &(ir(k,j,i,ifr*nang));
      Real *cosx = &(prad->mu(0,k,j,i,0));
      Real *cosy = &(prad->mu(1,k,j,i,0));
      Real *cosz = &(prad->mu(2,k,j,i,0));
#pragma omp simd aligned(cosx,weight,irn,cosy,cosz:ALI_LEN) reduction(+:er_f,fr1_f,fr2_f,fr3_f,pr11_f,pr12_f,pr13_f,pr22_f,pr23_f,pr33_f)
      for(int n=0; n<nang; ++n){
        Real irweight = weight[n] * irn[n];
        er_f   += irweight;
        fr1_f  += irweight * cosx[n];
        fr2_f  += irweight * cosy[n];
        fr3_f  += irweight * cosz[n];
        pr11_f += irweight * cosx[n] * cosx[n];
        pr12_f += irweight * cosx[n] * cosy[n];
        pr13_f += irweight * cosx[n] * cosz[n];
        pr22_f += irweight * cosy[n] * cosy[n];
        pr23_f += irweight * cosy[n] * cosz[n];
        pr33_f += irweight * cosz[n] * cosz[n];
      }
      er_f *= prad->wfreq(ifr);
      fr1_f *= prad->wfreq(ifr);
      fr2_f *= prad->wfreq(ifr);
      fr3_f *= prad->wfreq(ifr);
      pr11_f *= prad->wfreq(ifr);
      pr12_f *= prad->wfreq(ifr);
      pr13_f *= prad->wfreq(ifr);
      pr22_f *= prad->wfreq(ifr);
      pr23_f *= prad->wfreq(ifr);
      pr33_f *= prad->wfreq(ifr);

      er   += er_f;
      fr1  += fr1_f;
      fr2  += fr2_f;
      fr3  += fr3_f;
      pr11 += pr11_f;
      pr12 += pr12_f;
      pr13 += pr13_f;
      pr22 += pr22_f;
      pr23 += pr23_f;
      pr33 += pr33_f;

    }
  
    Real dtcsigma = ct * (prad->grey_sigma(OPAS,k,j,i) + prad->grey_sigma(OPAA,k,j,i));
  
    Real vx0 = (*vx);
    Real vy0 = (*vy);
    Real vz0 = (*vz);
  
    Real m0x = prat * fr1 * invcrat + rho * vx0;
    Real m0y = prat * fr2 * invcrat + rho * vy0;
    Real m0z = prat * fr3 * invcrat + rho * vz0;
  
  
    Real vx11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr11)
                                       * invcrat * invcrat;
    Real vy11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr22)
                                       * invcrat * invcrat;
    Real vz11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr33)
                                       * invcrat * invcrat;
    Real vx12 = dtcsigma * prat * pr12 * invcrat * invcrat;
    Real vx13 = dtcsigma * prat * pr13 * invcrat * invcrat;
    Real vy12 = dtcsigma * prat * pr23 * invcrat * invcrat;
    Real rhs1 = rho * vx0 + dtcsigma * m0x;
    Real rhs2 = rho * vy0 + dtcsigma * m0y;
    Real rhs3 = rho * vz0 + dtcsigma * m0z;
  
    Real factor = vx11 * vy11 * vz11 - vy11 * vx13 * vx13 + 2.0 * vx12 * vx13 * vy12
              - vx11 * vy12 * vy12 - vx12 * vx12 * vz11;
    factor = 1.0/factor;

    (*vx) = factor*(rhs3*(vx12*vy12 - vx13*vy11) + rhs2*(vy12*vx13
                - vx12*vz11) + rhs1*(vy11*vz11 - vy12*vy12));
          
    (*vy) = factor*(rhs3*(vx12*vx13 - vx11*vy12) + rhs2*(vx11*vz11
                - vx13*vx13) + rhs1*(vx13*vy12 - vx12*vz11));
          
    (*vz) = factor*(rhs3*(vx11*vy11 - vx12*vx12) + rhs2*(vx12*vx13
                - vx11*vy12) + rhs1*(vx12*vy12 - vx13*vy11));
  
}








