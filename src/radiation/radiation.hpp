#ifndef RADIATION_HPP_
#define RADIATION_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================


// C++ headers
#include <cstdint>     // int64_t
#include <functional>  // reference_wrapper
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/cc/radiation/bvals_rad.hpp"

class MeshBlock;
class ParameterInput;
class RadIntegrator;


//! \class Radiation
//  \brief radiation data and functions


// prototype for user-defined opacity function for radiative transfer

// Array indices for radiation moments
enum {IER=0, IFR1=1, IFR2=2, IFR3=3, IPR11=4, IPR22=5, IPR33=6, IPR12=7,
      IPR13=8, IPR23=9, IPR21=10, IPR31=11, IPR32=12};

enum {OPAS=0, OPAA=1, OPAP=2}; // scattering, absorption, Planck, opacity

class Radiation {
  friend class RadIntegrator;
public:
  Radiation(MeshBlock *pmb, ParameterInput *pin);
//  ~Radiation();
    
  AthenaArray<Real> ir, ir1, ir2; // radiation specific intensity
  AthenaArray<Real> rad_mom; // frequency integrated radiation moments
  AthenaArray<Real> rad_mom_cm; // co-moving frame Er, Frx, Fry, Frz
  AthenaArray<Real> sigma_s, sigma_a, sigma_ae; //   opacity
                       //sigma_a T and sigma_ae I
  AthenaArray<Real> sigma_planck; // absorption opacity to account for
                                 // the difference between Planck
                                // mean and Rosseland mean
  AthenaArray<Real> grey_sigma; // frequency integrated opacity
  AthenaArray<Real> mu, wmu; // angles and weight
  AthenaArray<Real> wfreq; // weight in frequency space
  
  AthenaArray<Real> flux[3]; // store transport flux, also need for refinement

  AthenaArray<Real> coarse_ir_;
  int refinement_idx{-1};
   
  Real prat, crat; // prat=aT^4/P_0, crat=c/c_s
  Real vmax;
  Real reduced_c; // reduced speed of light
  Real tunit, telectron; // gas temperature cgs unit,
                         // effective electron scattering temperature
  
  int nang, nfreq, noct, n_fre_ang; // n_fre_ang=nang*nfreq
  int angle_flag;
  // variables related to the angular space transport
  int nzeta, npsi;
  AthenaArray<Real> coszeta_v, zeta_v_full, zeta_f_full, dzeta_v, dzeta_f, 
                    coszeta_f, len_zeta;
  AthenaArray<Real> psi_v, psi_f, len_psi, psi_v_full, psi_f_full, 
                    dpsi_v, dpsi_f;


//  int ir_output; // the number of specific intensity to dump
//  AthenaArray<int> ir_index; // the array
//  AthenaArray<Real> dump_ir;

  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Fluid
  RadBoundaryVariable rad_bvar;
   
  
  RadIntegrator *pradintegrator;

    // The function pointer for the opacity
  OpacityFunc UpdateOpacity;
  
  int rotate_theta; // flag to rotate the boundary
  int rotate_phi;  
  int set_source_flag; // flag to add radiation source term or not


  // Functions   
  //Function in problem generators to update opacity
  void EnrollOpacityFunction(OpacityFunc MyOpacityFunction);


  
  //functin to calculate the radiation moments
  void CalculateMoment(AthenaArray<Real> &ir_in);
  void CalculateComMoment();
 
  void AngularGrid(int angle_flag, int nmu);
  void AngularGrid(int angle_flag, int nzeta, int npsi);

  void FrequencyGrid();

private:
  Real t_floor_; // temperature floor


  friend class BoundaryValues;
  

};

#endif // RADIATION_HPP_
