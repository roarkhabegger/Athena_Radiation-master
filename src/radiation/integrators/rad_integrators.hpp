#ifndef RADINTEGRATORS_HPP
#define RADINTEGRATORS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../radiation.hpp" // radiation
#include "../../task_list/task_list.hpp"

class MeshBlock;
class ParameterInput;
class Radiation;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer


class RadIntegrator {
  friend class Radiation;
public:
  RadIntegrator(Radiation *prad, ParameterInput *pin);
  ~RadIntegrator();
  
  Radiation *pmy_rad;
  
  void FluxDivergence(const Real wght, AthenaArray<Real> &ir_out);
    
  void CalculateFluxes(AthenaArray<Real> &w,
                       AthenaArray<Real> &ir, const int order);


  
  void AddSourceTerms(MeshBlock *pmb, const Real dt, AthenaArray<Real> &u,
        AthenaArray<Real> &w, AthenaArray<Real> &bcc, AthenaArray<Real> &ir);

  void AbsorptionScattering(AthenaArray<Real> &wmu_cm,
          AthenaArray<Real> &tran_coef, Real *sigma_a, Real *sigma_p,
          Real *sigma_ae, Real *sigma_s, Real dt, Real rho, Real &tgas,
          AthenaArray<Real> &ir_cm);
  
  void Compton(AthenaArray<Real> &wmu_cm,
          AthenaArray<Real> &tran_coef, Real *sigma_s,
          Real dt, Real rho, Real &tgas, AthenaArray<Real> &ir_cm);
  
  void LabToCom(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          Real *ir_lab, AthenaArray<Real> &ir_cm);
  
  void ComToLab(const Real vx, const Real vy, const Real vz,
                          Real *mux, Real *muy, Real *muz,
                          AthenaArray<Real> &ir_cm, Real *ir_lab);
  
  void ComAngle(const Real vx, const Real vy, const Real vz,
          Real mux, Real muy, Real muz, Real *mux0, Real *muy0, Real *muz0);

  
  void GetTaufactor(const Real vx, const Real vy, const Real vz,
                                 const Real ds, const Real sigma, Real *factor);

  void PredictVel(AthenaArray<Real> &ir, int k, int j, int i, Real dt, Real rho,
                  Real *vx, Real *vy, Real *vz);

  int rad_xorder; 
private:
  AthenaArray<Real> vel_, velx_,vely_,velz_;
  AthenaArray<Real> il_, ilb_, ir_;// for recontruction
                          // temporary array to store the flux, velocity
  AthenaArray<Real> temp_i1_, temp_i2_; // temporary array to store Div q
  AthenaArray<Real> vncsigma_, vncsigma2_, wmu_cm_, tran_coef_, ir_cm_;
  AthenaArray<Real> cm_to_lab_;
  AthenaArray<Real> g_zeta_, q_zeta_, ql_zeta_, qr_zeta_, zeta_flux_, zeta_area_;
  AthenaArray<Real> g_psi_, q_psi_, ql_psi_, qr_psi_, psi_flux_, psi_area_;
  AthenaArray<Real> dflx_ang_, ang_vol_;
                                    
 // temporary 1D array with size of nang
  Real taufact_;
  int compton_flag_; // flag to add simple Compton scattering
  int planck_flag_; // flag to add additional Planck absorption opacity
  int adv_flag_; // flag used to indicate whether separate
                 // advection flux from diffustion flux or not.

  int flux_correct_flag_; // flag to do second order flux crrection or not.
  Real tau_limit_; // the limit of optical depth sure.
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_, dflx_, cwidth2_, cwidth3_;

};

#endif // RADINTEGRATORS_HPP
