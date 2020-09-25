#ifndef BVALS_CC_RAD_BVALS_RAD_HPP_
#define BVALS_CC_RAD_BVALS_RAD_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_rad.hpp
//  \brief

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "../bvals_cc.hpp"


//----------------------------------------------------------------------------------------
//! \class CellCenteredBoundaryVariable
//  \brief

class RadBoundaryVariable : public CellCenteredBoundaryVariable {
 public:
  RadBoundaryVariable(MeshBlock *pmb,
                      AthenaArray<Real> *var_rad, AthenaArray<Real> *coarse_var,
                      AthenaArray<Real> *var_flux);
                                                // AthenaArray<Real> &prim);
  virtual ~RadBoundaryVariable() = default;

// functions unique implementation to radiation class
  void SendFluxCorrection() override;
  bool ReceiveFluxCorrection() override;
  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb) override;
  int LoadBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaries() override;
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) override;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) override;
  void PolarBoundarySingleAzimuthalBlock() override;

  // BoundaryPhysics: need to rotate the intensity
  void ReflectInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override;
  void ReflectOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override;
  void ReflectInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override;
  void ReflectOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override;
  void ReflectInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override;
  void ReflectOuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;

  void VacuumInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override;
  void VacuumOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override;
  void VacuumInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override;
  void VacuumOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override;
  void VacuumInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override;
  void VacuumOuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;


  void OutflowInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override;
  void OutflowOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override;
  void OutflowInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override;
  void OutflowOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override;
  void OutflowInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override;
  void OutflowOuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;
  
  //rotate the specific intensity by some amount
  //these functions are typically called with periodic boundary together
  void RotateHPi_InnerX2(Real time, Real dt, 
                      int il, int iu, int jl, int kl, int ku, int ngh);

  void RotateHPi_OuterX2(Real time, Real dt, 
                      int il, int iu, int ju, int kl, int ku, int ngh);

  void RotateHPi_InnerX3(Real time, Real dt, 
                      int il, int iu, int jl, int ju, int kl, int ngh);

  void RotateHPi_OuterX3(Real time, Real dt, 
                      int il, int iu, int jl, int ju, int ku, int ngh);

  void RotatePi_InnerX3(Real time, Real dt, 
                      int il, int iu, int jl, int ju, int kl, int ngh);
 
  void RotatePi_OuterX3(Real time, Real dt, 
                      int il, int iu, int jl, int ju, int ku, int ngh);
  
  void PolarWedgeInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;
  
  void PolarWedgeOuterX2(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;

private:
  AthenaArray<Real>  azimuthal_shift_rad_;


};

#endif // BVALS_CC_HYDRO_BVALS_HYDRO_HPP_
