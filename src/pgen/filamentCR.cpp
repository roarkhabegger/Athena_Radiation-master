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
#include "../fft/athena_fft.hpp"
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


//======================================================================================
/*! \file filamentCR.cpp
 *  problem for problem filamentary structure support through cosmic ray pressure 
 *  and examining multiphase CR turbulence
 *====================================================================================*/
//NOTES length scale is 1 pc
//      time scale is 1 Myr
//      density scale is m_p/cm^3
Real dens0, pres0, myGamma; // Initial hydro quantities
                       // dens0 in multiples of 1.67262192e-24 g / cm3
                       // pres0 in multiples of 1.59915640e-14 erg/cm^3
                       // myGamma is adiabatic index
Real phi; // background magnetic field angle wrt x1 axis, phi=pi/2 points B along x2 axis
Real beta;    // Ratio of gas pressure to magnetic pressure
Real beta_cr; // Ratio of gas pressure to cosmic ray pressure
Real vx; // background velocity along x axis

Real radC, widC; //Radius of cloud, and width of transistion. Both in length units 1pc
Real densC, presC, betaC, beta_crC; //Hydro quantities inside the cloud 

Real densProfile(Real x1, Real x2, Real x3);
Real presProfile(Real x1, Real x2, Real x3);
Real crProfile(Real x1, Real x2, Real x3);
Real magProfile(Real x1, Real x2, Real x3);

//cr Diffusion variables and function
Real sigmaParl, sigmaPerp; //CR diffusion 
                           //decouple parallel and perpendicular to the local magnetic field

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

Real densProfile(Real x1, Real x2, Real x3)
{
  Real cloud = 0.5*(1-tanh( (sqrt(SQR(x1) + SQR(x2) + SQR(x3)) - radC)/widC ));
  Real rho = 1.0 + cloud * (densC/dens0 - 1);
  return rho;
}

Real presProfile(Real x1, Real x2, Real x3)
{
  Real cloud = 0.5*(1-tanh( (sqrt(SQR(x1) + SQR(x2) + SQR(x3)) - radC)/widC ));
  Real pres = 1.0 + cloud * (presC/pres0-1);
  return pres;
}
Real crProfile(Real x1, Real x2, Real x3)
{
  Real cloud = 0.5*(1-tanh( (sqrt(SQR(x1) + SQR(x2) + SQR(x3)) - radC)/widC ));
  Real crP = 1.0 + cloud * (beta_cr/pres0 * presC/beta_crC - 1);
  return crP;
}
Real magProfile(Real x1, Real x2, Real x3)
{
  Real cloud =  0.5*(1-tanh( (sqrt(SQR(x1) + SQR(x2) + SQR(x3)) - radC)/widC ));
  Real myB = 1.0 + cloud * (sqrt(2*presC/betaC) * sqrt(beta/(2*pres0)) - 1);
  return myB;
}

// Define the function for Multi-Phase ODE Solver
Real Cooling(Real T, Real nH);
Real dTdt(Real T, Real nH);
Real AdaptiveODESolver(Real T, Real nH, Real dt);
Real CoolingTimeStep(MeshBlock *pmb);
const Real T_floor = 50;

// Set local (each cell independently) opacity funciton
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  if(CR_ENABLED){
   pcr->EnrollOpacityFunction(Diffusion);
  }
}

//Set up initial MESH data
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  myGamma = peos->GetGamma();
  Real gm1 = peos->GetGamma() - 1.0;

  // Load background variables
  pres0 = pin->GetReal("problem","pres0");
  dens0 = pin->GetReal("problem","dens0");
  beta = pin->GetReal("problem","beta");
  phi = (PI/180.0)pin->GetOrAddReal("problem","phi",0.0);
  beta_cr = pin->GetReal("problem","beta_cr");
  vx = pin->GetOrAddReal("problem","vx",0.0)

  // Load cloud variables
  presC = pin->GetOrAddReal("problem","presC",pres0);
  densC = pin->GetOrAddReal("problem","densC",dens0);
  betaC = pin->GetOrAddReal("problem","betaC",beta);
  beta_crC = pin->GetOrAddReal("problem","beta_crC",beta_cr);
  radC = pin->GetOrAddReal("problem","radC",5.0);
  widC = pin->GetOrAddReal("problem","radC",0.001);

  if(CR_ENABLED){
    //Load CR Variables
    sigmaPerp = pin->GetReal("cr","sigmaPerp");
    sigmaParl = pin->GetReal("cr","sigmaParl");
  }
  EnrollUserTimeStepFunction(CoolingTimeStep);
  turb_flag = pin->GetInteger("problem","turb_flag");
  if (turb_flag != 0) {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    ATHENA_ERROR(msg);
    return;
#endif
  }
}

//Setup initial mesh and variables
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real density = densProfile(x1,x2,x3);
        Real pressure = presProfile(x1,x2,x3);

        // set hydro variables
        phydro->u(IDN,k,j,i) = dens0*density;
        phydro->u(IM1,k,j,i) = dens0*density*vx;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = pres0*pressure/gm1+0.5*dens0*density*SQR(vx);
        
        if(CR_ENABLED){
          // get CR parameters
          Real CR = crProfile(x1,x2,x3);
          Real crp = pres0/beta_cr*CR; //*(1+ampCR*(1-x2/centerCR));
          pcr->u_cr(CRE,k,j,i) = 3.0*crp;
          pcr->u_cr(CRF1,k,j,i) = 0.0;
          pcr->u_cr(CRF2,k,j,i) = 0.0;
          pcr->u_cr(CRF3,k,j,i) = 0.0;
        }
      }// end i
    }
  }
  // initialize uniform interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    Real B0 = sqrt(2*pres0/beta);
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          Real x1 = pcoord->x1f(i);
          Real x2 = pcoord->x2v(j);
          Real x3 = pcoord->x3v(k);
          Real mag = magProfile(x1,x2,x3);
          pfield->b.x1f(k,j,i) = B0*mag*cos(phi);
        }
      }
    }
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2f(j);
          Real x3 = pcoord->x3v(k);
          Real mag = magProfile(x1,x2,x3);
          pfield->b.x2f(k,j,i) = B0*mag*sin(phi);
        }
      }
    }
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          Real x1 = pcoord->x1v(i);
          Real x2 = pcoord->x2v(j);
          Real x3 = pcoord->x3f(k);
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
          pcr->sigma_diff(1,k,j,i) = sigmaPerp/sqrt(2);
          pcr->sigma_diff(2,k,j,i) = sigmaPerp/sqrt(2);
        }
      }
    }// end k,j,i
  }// End CR
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkInLoop(ParameterInput *pin)
//  \brief Time-Intregator for adaptive Solver Method Solver
//========================================================================================
void MeshBlock::UserWorkInLoop() {
  const Real T_floor = 30.0;
  const Real kb  = 1.381e-16;
  const Real unit_length_in_cm_  = 3.086e+18;
  const Real unit_vel_in_cms_    = 1.0e5;
  const Real unit_density_in_nH_ = 1;
  const Real unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                             * unit_vel_in_cms_ * unit_vel_in_cms_;
  const Real unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  Real dt = pmy_mesh->dt*unit_time_in_s_;
  Real g = peos->GetGamma(); //Gamma
  Real nH,e,E_ergs,T,T_next,e_next;

  for (int k=ks-NGHOST; k<=ke+NGHOST; ++k) {
    for (int j=js-NGHOST; j<=je+NGHOST; ++j) {
#pragma omp simd
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
        nH     = phydro->u(IDN,k,j,i)*unit_density_in_nH_;
        e      = phydro->w(IPR,k,j,i)/(g-1.0);
        E_ergs = e* unit_E_in_cgs_ / nH;
        T      = E_ergs / (1.5*kb);
        T      = (T > T_floor) ? T : T_floor;
        T_next = AdaptiveODESolver(T,nH,dt);
        //std::cout<<"T = "<<T<<"T_next = "<<T_next<<std::endl;
        if (T_next < T_floor) {
          e_next = T_floor*(1.5*kb)*nH/unit_E_in_cgs_;
        }else{
          e_next = T_next*(1.5*kb)*nH/unit_E_in_cgs_;
        }
        phydro->u(IPR,k,j,i) += (e_next-e)*(g-1.0);
        phydro->w(IEN,k,j,i) += (e_next-e); 
      }
    }
  }
}

Real Cooling(Real T, Real nH){
  const Real HeatingRate = 2e-26;
  Real CoolingRate = 1e7*std::exp(-1.184e5/(T+1e3))+1.4e-2*std::sqrt(T)*std::exp(-92/T);
  Real dEdt = (T>T_floor)? nH*HeatingRate*(1.0 - nH*CoolingRate) : nH*HeatingRate;
  return nH*HeatingRate*(1.0 - nH*CoolingRate);
}

Real dTdt(Real T, Real nH){
  const Real kb  = 1.381e-16;
  Real dEdt = Cooling(T, nH);
  Real dTdt = dEdt/(kb*1.5);
  return dTdt;
}

//========================================================================================
//! \fn Real AdaptiveODESolver(Real T, Real nH, Real dt)
//  \brief Applying the Adaptive sub-cycle method O(h^2) to solve the cooling function 
//========================================================================================
Real AdaptiveODESolver(Real T, Real nH, Real dt) {
  const Real tol  = 1e-1;
  const Real tolmin = 1e-4;
  const Real minh = 1e-10;
  const int Maxiter = 100000;
  Real T0 = T;
  Real Tn,T_huan,err,t_reamin,s1,s2;
  Real  t = 0.0;
  Real  h = dt/2;
  int  iter = 0;
  while ( iter < Maxiter & t < dt){
    s1 = dTdt(T ,nH);
    Tn = T0 + h*s1;
    s2 = dTdt(Tn,nH);
    T_huan = T0 + h*(s1+s2)/2;
    err = std::abs(T_huan - Tn);
    if (err > tol ){
      // error is too large, recompute the cycle with smaller h
      h /= 2;
    }else if ( err < tolmin ){
      // take the result but increasing h by a factor of 2
      t += h;
      t_reamin = dt - t;
      T0 = T_huan;
      h  =  (t_reamin > 2*h )? 2*h : t_reamin;
    }else if ( tol > err > tolmin ){
      // take the result
      t += h;
      t_reamin = dt - t;
      T0 = T_huan;
      h  =  (t_reamin > h )? h : t_reamin;
    }else if (h<minh){
      std::stringstream msg;
      msg << "### FATAL ERROR in ProblemGenerator::AdaptiveODESolver" << std::endl
          << "h < minh, Solution is not convergence!" << std::endl
          << "T0 = "<< T <<", nH = "<<nH<<", dt ="<<dt << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    iter+=1;
  }

  if (iter == Maxiter){
    std::stringstream msg;
    msg << "### FATAL ERROR in ProblemGenerator::AdaptiveODESolver" << std::endl
        << "iter = Maxiter, Solution is not convergence!" << std::endl
        << "T0 = "<< T <<"nH = "<<nH<<", dt ="<<dt << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return Tn;
}

Real CoolingTimeStep(MeshBlock *pmb){
  const Real unit_length_in_cm_  = 3.086e+18;
  const Real unit_vel_in_cms_    = 1.0e5;
  const Real unit_density_in_nH_ = 1;
  const Real unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
                             * unit_vel_in_cms_ * unit_vel_in_cms_;
  const Real unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  const Real  g = 5.0/3.0;
  Real min_dt = 0.3*pmb->pcoord->dx1f(0)/std::abs(pmb->phydro->u(IVX,pmb->ke,pmb->je,pmb->ie));
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real   nH   = pmb->phydro->u(IDN,k,j,i)*unit_density_in_nH_;
        Real   ED   = pmb->phydro->w(IPR,k,j,i)/(g-1.0);
        Real E_ergs = ED * unit_E_in_cgs_ / nH;
        Real     T  =  E_ergs / (1.5*1.381e-16);
        Real Heating = 2e-26;
        Real Cooling = 2e-26*nH*(1e7*exp(-1.184e5/(T+ 1e3)) + 1.4e-2*sqrt(T)*exp(-92/T));
        Real dEdt    = Heating;
        if (T > T_floor) {
          dEdt = dEdt - Cooling;
        }
        Real cool_dt = std::abs(0.015*E_ergs/dEdt/unit_time_in_s_);
        if (min_dt > cool_dt){
          min_dt   = cool_dt;
        }
      }
    }
  }
  return min_dt;
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
