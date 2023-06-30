//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//  \brief Problem generator for turbulence generator
//

// C++ headers
#include <sstream>    // stringstream
#include <cmath>      // abs(), cos(), exp(), log(), NAN, pow(), sin(), sqrt()
#include <stdexcept>  // runtime_error
#include <ctime>
//#include <chrono>
#include <random>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../fft/athena_fft.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// Define the function for Multi-Phase ODE Solver
Real Cooling(Real T, Real nH);
Real dTdt(Real T, Real nH);
Real AdaptiveODESolver(Real T, Real nH, Real dt);
Real CoolingTimeStep(MeshBlock *pmb);
int cooling_flag;
const Real Heat    =  3.68962948e+01;
const Real Lamb1   =  1.34671476e+07;
const Real Lamb2   =  1.45740365e+01;
const Real T1a     =  9.77320931e+02;
const Real T1b     =  1.23815996e+01;
const Real T2      =  7.59404778e-01;
const Real T_floor =  4.12719988e-01;


Real sigmaParl, sigmaPerp; //CR diffusion 
                           //decouple parallel and perpendicular to the local magnetic field
Real crLoss;
void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr,
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);
void CRSource(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &u_cr);
void mySource(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons);


//========================================================================================
//! \fn void Mesh::InitUserMeshBlockData(ParameterInput *pin)
//  \brief
//========================================================================================
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  if(CR_ENABLED){
   pcr->EnrollOpacityFunction(Diffusion);
   pcr->EnrollUserCRSource(CRSource);
  }
}
//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {

  if(CR_ENABLED){
    //Load CR Variables
    sigmaPerp = pin->GetReal("cr","sigmaPerp");
    sigmaParl = pin->GetReal("cr","sigmaParl");
    crLoss = pin->GetOrAddReal("problem","crLoss",0.0);
  }
  cooling_flag = pin->GetInteger("problem","cooling");
  if (cooling_flag != 0) {
	  EnrollUserTimeStepFunction(CoolingTimeStep);
    EnrollUserExplicitSourceFunction(mySource);
  }
  // turb_flag is initialzed in the Mesh constructor to 0 by default;
  // turb_flag = 1 for decaying turbulence
  // turb_flag = 2 for driven turbulence
  turb_flag = pin->GetInteger("problem","turb_flag");
  if (turb_flag != 0) {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
#endif
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //dimensions of meshblock
  const int Nx = ie - is + 1;
  const int Ny = je - js + 1;
  const int Nz = ke - ks + 1;
  //read input parameters
  const Real nH = pin->GetReal("problem", "nH"); //density
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real pres = nH*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;
  
  const Real vx = pin->GetOrAddReal("problem", "vx", 0); //velocity x
  
  const Real invbeta = pin->GetOrAddReal("problem","invbeta",0.0);
  const Real b0 = sqrt(2*invbeta*pres); //mean field strength
  const Real angle = (PI/180.0)*pin->GetOrAddReal("problem","angle",0.0);
  const Real G0 = pin->GetOrAddReal("problem", "G0", 0.);

  const Real s_init = pin->GetOrAddReal("problem", "s_init", 0.);
  //mean and std of the initial gaussian profile

  const Real invbetaCR = pin->GetOrAddReal("problem","invbetaCR",0.0);
  const Real crpres = pres*invbetaCR;


  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        //density
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        //Real cloud = 0.5*(1-tanh( (sqrt(SQR(x1) + SQR(x2) + SQR(x3)) - rad_c)/wid_c ));
        phydro->u(IDN, k, j, i) = nH;//*(1.0 + cloud * (nH_c/nH - 1));
        //velocity, x, y, z direction
        phydro->u(IM1, k, j, i) = nH*vx;//*(1.0 + cloud * (nH_c/nH - 1));
        phydro->u(IM2, k, j, i) = 0.0;
        phydro->u(IM3, k, j, i) = 0.0;
        //energy
        if (NON_BAROTROPIC_EOS) {
          //phydro->u(IEN, k, j, i) = pres*(1.0 + cloud * (pres_c/pres-1))/gm1 + 0.5*nH*(1.0 + cloud * (nH_c/nH - 1))*SQR(vx);
          phydro->u(IEN, k, j, i) = pres/gm1 + 0.5*nH*SQR(vx);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          //phydro->u(IEN,k,j,i)+=0.5*SQR(b0*(1.0 + cloud * (b0_c/b0_c-1)));
          phydro->u(IEN,k,j,i)+=0.5*SQR(b0);
        }
        if(CR_ENABLED){
          pcr->u_cr(CRE,k,j,i) = 3.0*crpres;//*(1.0 + cloud * (crpres_c/crpres-1));
          pcr->u_cr(CRF1,k,j,i) = 0.0;
          pcr->u_cr(CRF2,k,j,i) = 0.0;
          pcr->u_cr(CRF3,k,j,i) = 0.0;
        }
      }
    }
  }

  //intialize magnetic field field
  if (MAGNETIC_FIELDS_ENABLED) {
    if (COORDINATE_SYSTEM == "cartesian") {

        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
            Real x1 = pcoord->x1f(i);
            Real x2 = pcoord->x2v(j);
            Real x3 = pcoord->x3v(k);
            pfield->b.x1f(k,j,i) = b0*std::cos(angle);
        }}}
        for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
            Real x1 = pcoord->x1v(i);
            Real x2 = pcoord->x2f(j);
            Real x3 = pcoord->x3v(k);
            //Real cloud = 0.5*(1-tanh( (sqrt(SQR(x1) + SQR(x2) + SQR(x3)) - rad_c)/wid_c ));
            pfield->b.x2f(k,j,i) = b0* std::sin(angle);//*(1.0 + cloud * (b0_c/b0_c-1));
        }}}
        for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
            // Real x1 = pcoord->x1v(i);
            // Real x2 = pcoord->x2v(j);
            // Real x3 = pcoord->x3f(k);
            //Real cloud = 0.5*(1-tanh( (sqrt(SQR(x1) + SQR(x2) + SQR(x3)) - rad_c)/wid_c ));
            pfield->b.x3f(k,j,i) = 0.0;//*(1.0 + cloud * (b0_c/b0_c-1)) ;
        }}}
    }else{
      std::stringstream msg;
    msg << "### FATAL ERROR in ProblemGenerator::MultiPhase Module" << std::endl
        << "Only support cartesian system!" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    }
  }
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {


}
//========================================================================================
//! \fn void Mesh::UserWorkInLoop(ParameterInput *pin)
//  \brief Time-Intregator for adaptive Solver Method Solver
//========================================================================================
void MeshBlock::UserWorkInLoop() {
//   const Real kb  = 1.381e-16;
//   const Real unit_length_in_cm_  = 3.086e+18;
//   const Real unit_vel_in_cms_    = 1.0e5;
//   const Real unit_density_in_nH_ = 1;
//   const Real unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
//                              * unit_vel_in_cms_ * unit_vel_in_cms_;
//   const Real unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
//   Real dt = pmy_mesh->dt*unit_time_in_s_;
//   Real g = peos->GetGamma(); //Gamma
//   Real nH,p,P_cgs,T,T_next,P_next;
//   if (cooling_flag != 0) {
//   for (int k=ks-NGHOST; k<=ke+NGHOST; ++k) {
//     for (int j=js-NGHOST; j<=je+NGHOST; ++j) {
// #pragma omp simd
//       for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
//         nH     = phydro->u(IDN,k,j,i)*unit_density_in_nH_;
//         p      = phydro->w(IPR,k,j,i);
//         P_cgs  = p * unit_E_in_cgs_ ;
//         // T      = E_ergs / (1.5*kb); //Why forced to be gamma=5/3 with 1.5 = 1/(g-1)
// 		    T      = P_cgs / (kb*nH);
//         T      = (T > T_floor) ? T : T_floor;
//         T_next = AdaptiveODESolver(T,nH,dt);
//         //std::cout<<"T = "<<T<<"T_next = "<<T_next<<std::endl;
//         if (T_next < T_floor) {
//           P_next = T_floor*kb*nH/unit_E_in_cgs_;
//         }else{
//           P_next = T_next*kb*nH/unit_E_in_cgs_;
//         }
//         phydro->w(IPR,k,j,i) += (P_next-P);
//         phydro->u(IEN,k,j,i) += (P_next-P)/(g-1); 
//       }
//     }
//   }
//   }
  return;
}

// Real Cooling(Real T, Real nH){
//   Real CoolingRate = 1e7*std::exp(-1.184e5/(T+1e3))+1.4e-2*std::sqrt(T)*std::exp(-92/T);
//   Real dEdt = (T>T_floor)? nH*HeatingRate*(1.0 - nH*CoolingRate) : nH*HeatingRate;
//   return dEdt;
// }

// Real dTdt(Real T, Real nH){
//   Real g = peos->GetGamma();
//   Real dEdt = Cooling(T, nH);
//   Real dTdt = dEdt/(kb)*(g-1);
//   return dTdt;
// }

//========================================================================================
//! \fn Real AdaptiveODESolver(Real T, Real nH, Real dt)
//  \brief Applying the Adaptive sub-cycle method O(h^2) to solve the cooling function 
//========================================================================================
// Real AdaptiveODESolver(Real T, Real nH, Real dt) {
//   const Real tol  = 1e-1;
//   const Real tolmin = 1e-4;
//   const Real minh = 1e-10;
//   const int Maxiter = 100000;
//   Real T0 = T;
//   Real Tn,T_huan,err,t_reamin,s1,s2;
//   Real  t = 0.0;
//   Real  h = dt/2;
//   int  iter = 0;
//   while ( iter < Maxiter & t < dt){
//     s1 = dTdt(T ,nH);
//     Tn = T0 + h*s1;
//     s2 = dTdt(Tn,nH);
//     T_huan = T0 + h*(s1+s2)/2;
//     err = std::abs(T_huan - Tn);
//     if (err > tol ){
//       // error is too large, recompute the cycle with smaller h
//       h /= 2;
//     }else if ( err < tolmin ){
//       // take the result but increasing h by a factor of 2
//       t += h;
//       t_reamin = dt - t;
//       T0 = T_huan;
//       h  =  (t_reamin > 2*h )? 2*h : t_reamin;
//     }else if ( tol > err > tolmin ){
//       // take the result
//       t += h;
//       t_reamin = dt - t;
//       T0 = T_huan;
//       h  =  (t_reamin > h )? h : t_reamin;
//     }else if (h<minh){
//       std::stringstream msg;
//       msg << "### FATAL ERROR in ProblemGenerator::AdaptiveODESolver" << std::endl
//           << "h < minh, Solution is not convergence!" << std::endl
//           << "T0 = "<< T <<", nH = "<<nH<<", dt ="<<dt << std::endl;
//       throw std::runtime_error(msg.str().c_str());
//     }
//     iter+=1;
//   }

//   if (iter == Maxiter){
//     std::stringstream msg;
//     msg << "### FATAL ERROR in ProblemGenerator::AdaptiveODESolver" << std::endl
//         << "iter = Maxiter, Solution is not convergence!" << std::endl
//         << "T0 = "<< T <<"nH = "<<nH<<", dt ="<<dt << std::endl;
//     throw std::runtime_error(msg.str().c_str());
//   }

//   return Tn;
// }

void mySource(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons){

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real d = prim(IDN,k,j,i);
        Real T = prim(IPR,k,j,i)/d;
        Real Lamb = Lamb1*exp(-1*T1a/(T + T1b)) + Lamb2*exp(-1*T2/T);
        Real dEdt = (T > T_floor) ? d*( d*Lamb - Heat ) : -1*d*Heat;
        cons(IEN,k,j,i) -= dEdt*dt;
      }
    }
  }
  return;

}


Real CoolingTimeStep(MeshBlock *pmb){

  Real min_dt=pmb->pcoord->dx1f(0)/pmb->pcr->vmax;
  Real g = pmb->peos->GetGamma();
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real d = pmb->phydro->w(IDN,k,j,i);
        Real T = pmb->phydro->w(IPR,k,j,i)/d;
        Real Lamb = Lamb1*exp(-1*T1a/(T + T1b)) + Lamb2*exp(-1*T2/T);
        Real dEdt = (T > T_floor) ? d*( d*Lamb - Heat ) : -1*d*Heat;
        Real cool_dt = abs( ( pmb->phydro->w(IPR,k,j,i)/(g-1) )/ dEdt);
        if (min_dt > cool_dt){
          min_dt   = cool_dt;
        }
      }
    }
  }

  // const Real unit_length_in_cm_  = 3.086e+18;
  // const Real unit_vel_in_cms_    = 1.0e5;
  // const Real unit_density_in_nH_ = 1;
  // const Real unit_E_in_cgs_ = 1.67e-24 * 1.4 * unit_density_in_nH_
  //                            * unit_vel_in_cms_ * unit_vel_in_cms_;
  // const Real unit_time_in_s_ = unit_length_in_cm_/unit_vel_in_cms_;
  // Real g = peos->GetGamma();
  // Real min_dt = 0.3*pmb->pcoord->dx1f(0)/std::abs(pmb->phydro->u(IVX,pmb->ke,pmb->je,pmb->ie));
  // for (int k=pmb->ks; k<=pmb->ke; ++k) {
  //   for (int j=pmb->js; j<=pmb->je; ++j) {
  //     for (int i=pmb->is; i<=pmb->ie; ++i) {
  //       Real   nH   = pmb->phydro->u(IDN,k,j,i)*unit_density_in_nH_;
  //       Real   ED   = pmb->phydro->w(IPR,k,j,i)/(g-1.0);
  //       Real E_ergs = ED * unit_E_in_cgs_ / nH;
  //       Real     T  =  E_ergs / (1.5*1.381e-16);
  //       Real Heating = 2e-26;
  //       Real Cooling = 2e-26*nH*(1e7*exp(-1.184e5/(T+ 1e3)) + 1.4e-2*sqrt(T)*exp(-92/T));
  //       Real dEdt    = Heating;
  //       if (T > T_floor) {
  //         dEdt = dEdt - Cooling;
  //       }
  //       Real cool_dt = std::abs(0.015*E_ergs/dEdt/unit_time_in_s_);
  //       if (min_dt > cool_dt){
  //         min_dt   = cool_dt;
  //       }
  //     }
  //   }
  // }
  return min_dt;
}


void CRSource(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &u_cr)
{ 
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
  #pragma omp simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        //CRLoss Term
        u_cr(CRE,k,j,i) -= crLoss*dt*u_cr(CRE,k,j,i);
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
