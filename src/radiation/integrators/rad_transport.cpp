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
//! \file rad_transport.cpp
//  \brief implementation of radiation integrators
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../radiation.hpp"
#include "../../coordinates/coordinates.hpp" //
#include "../../reconstruct/reconstruction.hpp"


// class header
#include "rad_integrators.hpp"



// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void RadIntegrator::CalculateFluxes(AthenaArray<Real> &w,
                     AthenaArray<Real> &ir, const int order)
{
  Radiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;
  Coordinates *pco=pmb->pcoord;
  
  int nang=prad->nang;
  int nfreq=prad->nfreq;
  Real invcrat=1.0/pmy_rad->crat;

  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2, 
  ncells3 = pmb->ncells3; 

  Real tau_fact;
  
  AthenaArray<Real> &x1flux=prad->flux[X1DIR];

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;
  jl = js, ju=je, kl=ks, ku=ke;

  if(ncells2 > 1)
  {
    if(ncells3 == 1){
      jl=js-1, ju=je+1, kl=ks, ku=ke;
    }else{
      jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
 
  }


    
    
//--------------------------------------------------------------------------------------
// i-direction



  // This part is used for all three directions when adv_flag is turned on
  if(adv_flag_ > 0){

    for(int k=0; k<ncells3; ++k){
      for(int j=0; j<ncells2; ++j){
          if (ncells2 > 1) pco->CenterWidth2(k,j,0,ncells1-1,cwidth2_);
          if (ncells3 > 1) pco->CenterWidth3(k,j,0,ncells1-1,cwidth3_);
        for(int i=0; i<ncells1; ++i){
           Real vx = w(IVX,k,j,i);
           Real vy = w(IVY,k,j,i);
           Real vz = w(IVZ,k,j,i);
           for(int ifr=0; ifr<nfreq; ++ifr){
             Real ds = pco->dx1v(i);
             if(ncells2 > 1) ds += cwidth2_(i);
             if(ncells3 > 1) ds += cwidth3_(i);
               
             GetTaufactor(vx, vy, vz, ds,
                      prad->sigma_a(k,j,i,ifr)+prad->sigma_s(k,j,i,ifr), &tau_fact);
             Real *cosx = &(prad->mu(0,k,j,i,0));
             Real *cosy = &(prad->mu(1,k,j,i,0));
             Real *cosz = &(prad->mu(2,k,j,i,0));
             Real *irn = &(ir(k,j,i,ifr*nang));
             Real *tempi1n = &(temp_i1_(k,j,i,ifr*nang));
             Real *tempi2n = &(temp_i2_(k,j,i,ifr*nang));
#pragma omp simd aligned(cosx,cosy,cosz,irn,tempi1n,tempi2n:ALI_LEN)
             for(int n=0; n<nang; ++n){
               
               Real vdotn = vx*cosx[n]+vy*cosy[n]+vz*cosz[n];
                 
               vdotn *= invcrat;
                 
               Real adv_coef = tau_fact * vdotn * (3.0 + vdotn * vdotn);
               Real q1 = irn[n] * (1.0 - adv_coef);
               tempi1n[n] = q1;
               tempi2n[n] = adv_coef;
               
             }
          }
        }// end i
      }// end j
    }// end k

    for (int k=kl; k<=ku; ++k){
      for (int j=jl; j<=ju; ++j){
        // get the velocity at the interface
        for(int i=is-1; i<=ie+1; ++i){
          Real dxl = pco->x1f(i)-pco->x1v(i-1);
          Real dxr = pco->x1v(i) - pco->x1f(i);
          Real factl = dxr/(dxl+dxr);
          Real factr = dxl/(dxl+dxr);
          for(int ifr=0; ifr<nfreq; ++ifr){
            Real *cosx = &(prad->mu(0,k,j,i-1,0));
            Real *cosx1 = &(prad->mu(0,k,j,i,0));
            Real *tempi2 = &(temp_i2_(k,j,i-1,ifr*nang));
            Real *tempi2_1 = &(temp_i2_(k,j,i,ifr*nang));
            Real *veln = &(vel_(k,j,i,ifr*nang));
#pragma omp simd aligned(cosx,cosx1,tempi2,tempi2_1,veln:ALI_LEN)
            for(int n=0; n<nang; ++n){
              // linear intepolation between x1v(i-1), x1f(i), x1v(i)
              veln[n] = prad->reduced_c *
                    (factl * cosx[n] * tempi2[n]
                     + factr * cosx1[n] * tempi2_1[n]);
            }


          }// end ifr
        }// end i
      }// end j
    }// endk

  } // End adv_flag


    // First, do reconstruction of the specific intensities

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      // reconstruct L/R states
      if(adv_flag_ > 0)
      {
        if (order == 1) {
          pmb->precon->DonorCellX1(k, j, is-1, ie+1, temp_i1_, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, temp_i1_, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, temp_i1_, -1, il_, ir_);
        }
      }else{
        if (order == 1) {
          pmb->precon->DonorCellX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
        }        
      }

      // calculate flux with velocity times the interface state
      for(int i=is; i<=ie+1; ++i){
        Real *vel = &(velx_(k,j,i,0));
        for(int n=0; n<prad->n_fre_ang; ++n){
          if(vel[n] > 0.0)
            x1flux(k,j,i,n) = vel[n] * il_(i,n);
          else if(vel[n] < 0.0)
            x1flux(k,j,i,n) = vel[n] * ir_(i,n);
          else
            x1flux(k,j,i,n) = 0.0;
        }// end n
      }// end i

    }// end j
  }// end k

  if(adv_flag_ > 0){

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        // reconstruct L/R states
        if (order == 1) {
          pmb->precon->DonorCellX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX1(k, j, is-1, ie+1, ir, -1, il_, ir_);
        }

        // calculate flux with velocity times the interface state
        for(int i=is; i<=ie+1; ++i){
          Real *vel = &(vel_(k,j,i,0));
          for(int n=0; n<prad->n_fre_ang; ++n){
            if(vel[n] > 0.0)
              x1flux(k,j,i,n) += vel[n] * il_(i,n);
            else if(vel[n] < 0.0)
              x1flux(k,j,i,n) += vel[n] * ir_(i,n);

          }// end n
        }// end i

      }// end j
    }// end k

  }// end adv_flag

// j-direction
  if(pmb->pmy_mesh->f2){
    AthenaArray<Real> &x2flux=prad->flux[X2DIR];

    il = is-1, iu = ie+1, kl = ks, ku = ke;
    if (ncells3 ==  1) // 2D
      kl = ks, ku = ke;
    else // 3D
      kl = ks-1, ku = ke+1;


    if(adv_flag_ > 0){

      for (int k=kl; k<=ku; ++k){
#pragma omp for schedule(static)
        for (int j=js; j<=je+1; ++j){
        // get the velocity
          for(int i=il; i<=iu; ++i){
            Real dxl = pco->x2f(j)-pco->x2v(j-1);
            Real dxr = pco->x2v(j) - pco->x2f(j);
            Real factl = dxr/(dxl+dxr);
            Real factr = dxl/(dxl+dxr);
            for(int ifr=0; ifr<nfreq; ++ifr){
              Real *cosy = &(prad->mu(1,k,j-1,i,0));
              Real *cosy1 = &(prad->mu(1,k,j,i,0));
              Real *tempi2 = &(temp_i2_(k,j-1,i,ifr*nang));
              Real *tempi2_1 = &(temp_i2_(k,j,i,ifr*nang));
              Real *veln = &(vel_(k,j,i,ifr*nang));
#pragma omp simd aligned(cosy,cosy1,tempi2,tempi2_1,veln:ALI_LEN)
              for(int n=0; n<nang; ++n){
              // linear intepolation between x2v(j-1), x2f(j), x2v(j)                 
                veln[n] = prad->reduced_c *
                    (factl * cosy[n] * tempi2[n]
                     + factr * cosy1[n] * tempi2_1[n]);
              }
            }// end ifr
          }// end i
        }// end j
      }// end k
    }// end adv_flag

    for (int k=kl; k<=ku; ++k) {
      //reconstruct the first row
      if(adv_flag_ > 0){
        if (order == 1) {
          pmb->precon->DonorCellX2(k, js-1, il, iu, temp_i1_, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, temp_i1_, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, temp_i1_, -1, il_, ir_);
        }
      }else{// end adv_flag
        if (order == 1) {
          pmb->precon->DonorCellX2(k, js-1, il, iu, ir, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, ir, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, ir, -1, il_, ir_);
        }        
      }

      for(int j=js; j<=je+1; ++j){

        if(adv_flag_ > 0){

          if (order == 1) {
            pmb->precon->DonorCellX2(k, j, il, iu, temp_i1_, -1, ilb_, ir_);
          } else if (order == 2) {
            pmb->precon->PiecewiseLinearX2(k, j, il, iu, temp_i1_, -1, ilb_, ir_);
          } else {
            pmb->precon->PiecewiseParabolicX2(k, j, il, iu, temp_i1_, -1, ilb_, ir_);
          }

        }else{
          if (order == 1) {
            pmb->precon->DonorCellX2(k, j, il, iu, ir, -1, ilb_, ir_);
          } else if (order == 2) {
            pmb->precon->PiecewiseLinearX2(k, j, il, iu, ir, -1, ilb_, ir_);
          } else {
            pmb->precon->PiecewiseParabolicX2(k, j, il, iu, ir, -1, ilb_, ir_);
          }       
        }

        //calculate the flux
      // calculate flux with velocity times the interface state
        for(int i=il; i<=iu; ++i){
          Real *vel = &(vely_(k,j,i,0));
          for(int n=0; n<prad->n_fre_ang; ++n){
            if(vel[n] > 0.0)
              x2flux(k,j,i,n) = vel[n] * il_(i,n);
            else if(vel[n] < 0.0)
              x2flux(k,j,i,n) = vel[n] * ir_(i,n);
            else
              x2flux(k,j,i,n) = 0.0;
          }// end n
        }// end i

        //swap the array for the next step

        il_.SwapAthenaArray(ilb_);

      }// end j loop from js to je+1

      if(adv_flag_ > 0){

        if (order == 1) {
          pmb->precon->DonorCellX2(k, js-1, il, iu, ir, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX2(k, js-1, il, iu, ir, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX2(k, js-1, il, iu, ir, -1, il_, ir_);
        }

        for(int j=js; j<=je+1; ++j){
          if (order == 1) {
            pmb->precon->DonorCellX2(k, j, il, iu, ir, -1, ilb_, ir_);
          } else if (order == 2) {
            pmb->precon->PiecewiseLinearX2(k, j, il, iu, ir, -1, ilb_, ir_);
          } else {
            pmb->precon->PiecewiseParabolicX2(k, j, il, iu, ir, -1, ilb_, ir_);
          }

        //calculate the flux
      // calculate flux with velocity times the interface state
          for(int i=il; i<=iu; ++i){
            Real *vel = &(vel_(k,j,i,0));
            for(int n=0; n<prad->n_fre_ang; ++n){
              if(vel[n] > 0.0)
                x2flux(k,j,i,n) += vel[n] * il_(i,n);
              else if(vel[n] < 0.0)
                x2flux(k,j,i,n) += vel[n] * ir_(i,n);
            }// end n
          }// end i

        //swap the array for the next step

          il_.SwapAthenaArray(ilb_);

        }// end j loop from js to je+1

      }// end adv_flag_ > 0


    }// end k



  }//end pmy_mesh->f2 
    
  if(pmb->pmy_mesh->f3){
    AthenaArray<Real> &x3flux=prad->flux[X3DIR];

    il =is-1, iu=ie+1, jl=js-1, ju=je+1;

    // First, calculate the transport velocity

    if(adv_flag_ > 0){

      for (int k=ks; k<=ke+1; ++k){
#pragma omp for schedule(static)
        for (int j=jl; j<=ju; ++j){
        // get the velocity
          for(int i=il; i<=iu; ++i){
            Real dxl = pco->x3f(k) - pco->x3v(k-1);
            Real dxr = pco->x3v(k) - pco->x3f(k);
            Real factl = dxr/(dxl+dxr);
            Real factr = dxl/(dxl+dxr);
            for(int ifr=0; ifr<nfreq; ++ifr){
              Real *cosz = &(prad->mu(2,k-1,j,i,0));
              Real *cosz1 = &(prad->mu(2,k,j,i,0));
              Real *tempi2 = &(temp_i2_(k-1,j,i,ifr*nang));
              Real *tempi2_1 = &(temp_i2_(k,j,i,ifr*nang));
              Real *veln = &(vel_(k,j,i,ifr*nang));
#pragma omp simd aligned(cosz,cosz1,tempi2,tempi2_1,veln:ALI_LEN)
              for(int n=0; n<nang; ++n){
              // linear intepolation between x2v(j-1), x2f(j), x2v(j)                 
                veln[n] = prad->reduced_c *
                  (factl * cosz[n] * tempi2[n] + factr * cosz1[n] * tempi2_1[n]);
              }
            }// end ifr
          }
        }// end j
      }// end k
    }// end adv_flag

    for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
      if(adv_flag_ > 0){
        if (order == 1) {
          pmb->precon->DonorCellX3(ks-1, j, il, iu, temp_i1_, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, temp_i1_, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, temp_i1_, -1, il_, ir_);
        }
      }else{
        if (order == 1) {
          pmb->precon->DonorCellX3(ks-1, j, il, iu, ir, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, ir, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, ir, -1, il_, ir_);
        }        
      }

      for (int k=ks; k<=ke+1; ++k) {
        // reconstruct L/R states at k
        if(adv_flag_ > 0){
          if (order == 1) {
            pmb->precon->DonorCellX3(k, j, il, iu, temp_i1_, -1, ilb_, ir_);
          } else if (order == 2) {
            pmb->precon->PiecewiseLinearX3(k, j, il, iu, temp_i1_, -1, ilb_, ir_);
          } else {
            pmb->precon->PiecewiseParabolicX3(k, j, il, iu, temp_i1_, -1, ilb_, ir_);
          }
        }else{
          if (order == 1) {
            pmb->precon->DonorCellX3(k, j, il, iu, ir, -1, ilb_, ir_);
          } else if (order == 2) {
            pmb->precon->PiecewiseLinearX3(k, j, il, iu, ir, -1, ilb_, ir_);
          } else {
            pmb->precon->PiecewiseParabolicX3(k, j, il, iu, ir, -1, ilb_, ir_);
          }

        }

      // calculate flux with velocity times the interface state
        for(int i=il; i<=iu; ++i){
            Real *vel = &(velz_(k,j,i,0));
          for(int n=0; n<prad->n_fre_ang; ++n){
            if(vel[n] > 0.0)
                x3flux(k,j,i,n) = vel[n] * il_(i,n);
            else if(vel[n] < 0.0)
                x3flux(k,j,i,n) = vel[n] * ir_(i,n);
            else
                x3flux(k,j,i,n) = 0.0;
          }// end n
        }// end i        

        // swap the arrays for the next step
        il_.SwapAthenaArray(ilb_);
      }// end k loop
    }// end j from jl to ju

    if(adv_flag_ > 0){

      for (int j=jl; j<=ju; ++j) { // this loop ordering is intentional
      // reconstruct the first row
        if (order == 1) {
          pmb->precon->DonorCellX3(ks-1, j, il, iu, ir, -1, il_, ir_);
        } else if (order == 2) {
          pmb->precon->PiecewiseLinearX3(ks-1, j, il, iu, ir, -1, il_, ir_);
        } else {
          pmb->precon->PiecewiseParabolicX3(ks-1, j, il, iu, ir, -1, il_, ir_);
        }
        for (int k=ks; k<=ke+1; ++k) {
          // reconstruct L/R states at k
          if (order == 1) {
            pmb->precon->DonorCellX3(k, j, il, iu, ir, -1, ilb_, ir_);
          } else if (order == 2) {
            pmb->precon->PiecewiseLinearX3(k, j, il, iu, ir, -1, ilb_, ir_);
          } else {
            pmb->precon->PiecewiseParabolicX3(k, j, il, iu, ir, -1, ilb_, ir_);
          }

      // calculate flux with velocity times the interface state
          for(int i=il; i<=iu; ++i){
              Real *vel = &(vel_(k,j,i,0));
            for(int n=0; n<prad->n_fre_ang; ++n){
              if(vel[n] > 0.0)
                  x3flux(k,j,i,n) += vel[n] * il_(i,n);
              else if(vel[n] < 0.0)
                  x3flux(k,j,i,n) += vel[n] * ir_(i,n);
            }// end n
          }// end i        

        // swap the arrays for the next step
          il_.SwapAthenaArray(ilb_);
        }// end k loop
      }// end j from jl to ju
    }// end adv_flag > 0


  }// end k direction

  // calculate flux along angular direction
  if(prad->angle_flag == 1 && prad->nzeta > 0){
    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){
          int &nzeta = prad->nzeta;
          int &npsi = prad->npsi;
          int psi_limit=2*npsi;
          if(npsi == 0) psi_limit=1;

          for(int m=0; m<psi_limit; ++m){
#pragma omp simd
            for(int n=0; n<nzeta*2; ++n){
              int ang_num = n*psi_limit+m;
              q_zeta_(n+NGHOST) = ir(k,j,i,ang_num);
            }// end nzeta
            // Add ghost zones
#pragma omp simd
            for(int n=1; n<=NGHOST; ++n){
              q_zeta_(NGHOST-n) = q_zeta_(NGHOST+n-1);
            }      
#pragma omp simd
            for(int n=1; n<=NGHOST; ++n){
              q_zeta_(2*nzeta+NGHOST+n-1) = q_zeta_(2*nzeta+NGHOST-n);
            }  
            int zl = NGHOST-1;
            int zu = 2*nzeta+NGHOST;
            if (order == 1) {
              pmb->precon->DonorCellZeta(prad,zl,zu,q_zeta_,
                                          ql_zeta_,qr_zeta_);
            } else {
              pmb->precon->PiecewiseLinearZeta(prad,zl,zu,q_zeta_,
                                          ql_zeta_,qr_zeta_);
            } 
  
            // zeta flux
            pco->GetGeometryZeta(prad,k,j,i,g_zeta_);
            for(int n=0; n<nzeta*2+1; ++n){
              int ang_num = n*psi_limit+m;
              Real g_coef = g_zeta_(n);
              if(g_coef > 0)
                zeta_flux_(k,j,i,ang_num) =  -prad->reduced_c * qr_zeta_(n+NGHOST) 
                                        * g_zeta_(n); 
              else if(g_coef < 0)
                zeta_flux_(k,j,i,ang_num) =  -prad->reduced_c * ql_zeta_(n+NGHOST) 
                                       * g_zeta_(n); 

            }
          }// end npsi
        }
      }
    }


  }//end angule_flag==1 and nzeta > 0


  // Now calculate psi flux
  if(prad->angle_flag == 1 && prad->npsi > 0){
    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){
          int &nzeta = prad->nzeta;
          int &npsi = prad->npsi;
          int zeta_limit=2*nzeta;
          if(nzeta == 0) zeta_limit=1;

          for(int n=0; n<zeta_limit; ++n){
            Real *qpsi = &(q_psi_(NGHOST));
            Real *irm = &(ir(k,j,i,n*2*npsi));
#pragma omp simd aligned(qpsi,irm:ALI_LEN)
            for(int m=0; m<npsi*2; ++m){
              qpsi[m] = irm[m];
            }// end nzeta
            // Add ghost zones
            // phi is periodic
            // the first one qpsi[NGHOST]
            // The last one is qpsi[NGHOST+2*npsi-1]
            for(int m=1; m<=NGHOST; ++m){
              q_psi_(NGHOST-m) = q_psi_(2*npsi+NGHOST-m);
            }      
            for(int m=1; m<=NGHOST; ++m){
              q_psi_(2*npsi+NGHOST+m-1) = q_psi_(NGHOST+m-1);
            }  
            int pl = NGHOST-1;
            int pu = 2*npsi+NGHOST;
            if (order == 1) {
              pmb->precon->DonorCellPsi(prad,pl,pu,q_psi_,
                                          ql_psi_,qr_psi_);
            } else {
              pmb->precon->PiecewiseLinearPsi(prad,pl,pu,q_psi_,
                                          ql_psi_,qr_psi_);
            } 
 
            // psi flux
            if(nzeta > 0)
              pco->GetGeometryPsi(prad,k,j,i,n,g_psi_);
            else
              pco->GetGeometryPsi(prad,k,j,i,g_psi_);
            for(int m=0; m<npsi*2+1; ++m){
              int ang_num = n*2*npsi+m;
              Real g_coef = g_psi_(m);
              if(g_coef > 0)
                psi_flux_(k,j,i,ang_num) =  -prad->reduced_c * qr_psi_(m+NGHOST) 
                                        * g_coef;
              else if(g_coef < 0)
                psi_flux_(k,j,i,ang_num) =  -prad->reduced_c * ql_psi_(m+NGHOST) 
                                        * g_coef;
            }
          }// end nzeta
        }
      }
    }


  }// end ang_flag==1 and npsi > 0

  
}// end calculate_flux


void RadIntegrator::FluxDivergence(const Real wght, AthenaArray<Real> &ir_out)
{
  Radiation *prad=pmy_rad;
  MeshBlock *pmb=prad->pmy_block;

  AthenaArray<Real> &x1flux=prad->flux[X1DIR];
  AthenaArray<Real> &x2flux=prad->flux[X2DIR];
  AthenaArray<Real> &x3flux=prad->flux[X3DIR];

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  
  AthenaArray<Real> &x1area = x1face_area_, &x2area = x2face_area_,
                 &x2area_p1 = x2face_area_p1_, &x3area = x3face_area_,
                 &x3area_p1 = x3face_area_p1_, &vol = cell_volume_, &dflx = dflx_;

  AthenaArray<Real> &area_zeta = zeta_area_, &area_psi = psi_area_, 
                 &ang_vol = ang_vol_, &dflx_ang = dflx_ang_;
  int &nzeta = prad->nzeta, &npsi = prad->npsi;

  for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) {

      // calculate x1-flux divergence 
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
      for(int i=is; i<=ie; ++i){
        Real *flxr = &(x1flux(k,j,i+1,0));
        Real *flxl = &(x1flux(k,j,i,0));
        Real *flxn = &(dflx(i,0));
#pragma omp simd aligned(flxn,flxr,flxl:ALI_LEN)
        for(int n=0; n<prad->n_fre_ang; ++n){
          flxn[n] = (x1area(i+1) *flxr[n] - x1area(i)*flxl[n]);
        }// end n
      }// End i

     // calculate x2-flux
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
      for(int i=is; i<=ie; ++i){
        Real *flxr = &(x2flux(k,j+1,i,0));
        Real *flxl = &(x2flux(k,j,i,0));
        Real *flxn = &(dflx(i,0));
#pragma omp simd aligned(flxr,flxl,flxn:ALI_LEN) 
        for(int n=0; n<prad->n_fre_ang; ++n){
            flxn[n] += (x2area_p1(i)*flxr[n] - x2area(i)*flxl[n]);
          }
        }
      }// end nx2

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
      for(int i=is; i<=ie; ++i){
        Real *flxr = &(x3flux(k+1,j,i,0));
        Real *flxl = &(x3flux(k,j,i,0));
        Real *flxn = &(dflx(i,0));
#pragma omp simd aligned(flxn,flxr,flxl:ALI_LEN) 
        for(int n=0; n<prad->n_fre_ang; ++n){
            flxn[n] += (x3area_p1(i)*flxr[n] - x3area(i)*flxl[n]);
          }
        }
      }// end nx3
      // update variable with flux divergence
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for(int i=is; i<=ie; ++i){
        Real *iro = &(ir_out(k,j,i,0));
        Real *flxn = &(dflx(i,0));
#pragma omp simd aligned(iro,flxn:ALI_LEN)
        for(int n=0; n<prad->n_fre_ang; ++n){
          iro[n] = std::max(iro[n]-wght*flxn[n]/vol(i), TINY_NUMBER);
        }
      }

      // add angular flux
      if(prad->angle_flag == 1){
        for(int i=is; i<=ie; ++i){
          for(int n=0; n<prad->nang; ++n)
            dflx_ang(n) = 0.0;
          if(nzeta * npsi > 0){
            for(int m=0; m<2*npsi; ++m){
#pragma omp simd 
              for(int n=0; n<2*nzeta; ++n){
                int ang_num = n*2*npsi + m;
                int ang_num1 = (n+1)*2*npsi+m;
                dflx_ang(ang_num) += (area_zeta(m,n+1) * zeta_flux_(k,j,i,ang_num1)
                                     - area_zeta(m,n) * zeta_flux_(k,j,i,ang_num));
              }// end zeta angle
            }// end psi angles
            // now psi flux
            for(int n=0; n<2*nzeta; ++n){
              Real *flxn = &(dflx_ang(n*2*npsi));
              Real *areapsi = &(area_psi(n,0));
              Real *psiflx = &(psi_flux_(k,j,i,n*2*npsi));
#pragma omp simd aligned(flxn,areapsi,psiflx:ALI_LEN)
              for(int m=0; m<2*npsi; ++m){
                flxn[m] += (areapsi[m+1] * psiflx[m+1]
                          - areapsi[m] * psiflx[m]);                
              }
            }
          }else if(nzeta >0){// end if nzeta*npsi > 0
            Real *flxn = &(dflx_ang(0));
            Real *areazeta = &(area_zeta(0));
            Real *zetaflx = &(zeta_flux_(k,j,i,0));
#pragma omp simd aligned(flxn,areazeta,zetaflx:ALI_LEN)
            for(int n=0; n<2*nzeta; ++n){         
              flxn[n] += (areazeta[n+1] * zetaflx[n+1]
                        - areazeta[n] * zetaflx[n]);
            }// end zeta angle
          }else if(npsi > 0){
            Real *flxn = &(dflx_ang(0));
            Real *areapsi = &(area_psi(0));
            Real *psiflx = &(psi_flux_(k,j,i,0));
#pragma omp simd aligned(flxn,areapsi,psiflx:ALI_LEN) 
            for(int m=0; m<2*npsi; ++m){
              flxn[m] += (areapsi[m+1] * psiflx[m+1] 
                        - areapsi[m] * psiflx[m]);              
            }           
          }// end npsi > 0
          // apply the flux divergence back
          Real *iro = &(ir_out(k,j,i,0));
          Real *flxn = &(dflx_ang(0));
          Real *angv = &(ang_vol(0));
#pragma omp simd aligned(iro,flxn,angv:ALI_LEN)
          for(int n=0; n<prad->nang; ++n){
            iro[n] = std::max(iro[n]-wght*flxn[n]/angv[n], TINY_NUMBER);

          }// end angle
        }// end i
      }// end if angle_flag == 1
 
    }// end j
  }// End k
  
}

void RadIntegrator::GetTaufactor(const Real vx, const Real vy, const Real vz,
                                 const Real ds, const Real sigma, Real *factor)
{

   Real invcrat = 1.0/pmy_rad->crat;
   Real vel = vx*vx+vy*vy+vz*vz;
   Real taucell = taufact_ * ds * sigma;
   Real tausq = taucell * taucell * (vel*invcrat*invcrat);
   if(tausq > 1.e-3) tausq = 1.0 - exp(-tausq);
   (*factor) = tausq;

}
