//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_rad.cpp
//  \brief implements boundary functions for Hydro variables and utilities to manage
//  primitive/conservative variable relationship in a derived class of the
//  CellCenteredBoundaryVariable base class.

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../globals.hpp"
#include "bvals_rad.hpp"

//----------------------------------------------------------------------------------------
//! \class RadiationBoundaryFunctions

RadBoundaryVariable::RadBoundaryVariable(MeshBlock *pmb, 
    AthenaArray<Real> *var_rad, AthenaArray<Real> *coarse_var,
    AthenaArray<Real> *var_flux) :
    CellCenteredBoundaryVariable(pmb, var_rad, coarse_var, var_flux, 1){

    // the radiation array is (k,j,i,n)
    // the number of variables is GetDim1
    // All the other shared functions are initialized in CellCenteredBoundaryVariable
    // radiation specific functions need to be defined here

    azimuthal_shift_rad_.NewAthenaArray(pmb->ke + NGHOST + 2,nu_+1);

}



//radiation needs to send and receive flux in differnt order
//over ride the SendFluxCorrection function in the base type
void RadBoundaryVariable::SendFluxCorrection()
{

  MeshBlock *pmb = pmy_block_;
  Coordinates *pco = pmb->pcoord;

  // cache pointers to surface area arrays (BoundaryBase protected variable)
  AthenaArray<Real> &sarea0 = pbval_->sarea_[0];
  AthenaArray<Real> &sarea1 = pbval_->sarea_[1];

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.ni.type != NeighborConnect::face) break;
    if (bd_var_flcor_.sflag[nb.bufid] == BoundaryStatus::completed) continue;
    if (nb.snb.level == pmb->loc.level - 1) {
      int p = 0;
      Real *sbuf = bd_var_flcor_.send[nb.bufid];
      // x1 direction
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int i = pmb->is + (pmb->ie-pmb->is + 1)*nb.fid;
        if (pmb->block_size.nx3>1) { // 3D
          for (int k=pmb->ks; k<=pmb->ke; k+=2) {
            for (int j=pmb->js; j<=pmb->je; j+=2) {
              Real amm = pco->GetFace1Area(k,   j,   i);
              Real amp = pco->GetFace1Area(k,   j+1, i);
              Real apm = pco->GetFace1Area(k+1, j,   i);
              Real app = pco->GetFace1Area(k+1, j+1, i);
              Real tarea = amm + amp + apm + app;
              for (int nn=nl_; nn<=nu_; nn++) {
                sbuf[p++] = (x1flux(k  , j  , i, nn)*amm
                            + x1flux(k  , j+1, i, nn)*amp
                            + x1flux(k+1, j  , i, nn)*apm
                            + x1flux(k+1, j+1, i, nn)*app)/tarea;
              }
            }
          }
        } else if (pmb->block_size.nx2>1) { // 2D
          int k = pmb->ks;          
          for (int j=pmb->js; j<=pmb->je; j+=2) {
            Real am = pco->GetFace1Area(k, j,   i);
            Real ap = pco->GetFace1Area(k, j+1, i);
            Real tarea = am + ap;
            for (int nn=nl_; nn<=nu_; nn++) {
              sbuf[p++] = (x1flux(k, j  , i,nn)*am + x1flux(k, j+1, i,nn)*ap)/tarea;
            }
          }
        } else { // 1D
          int k = pmb->ks, j = pmb->js;
          for (int nn=nl_; nn<=nu_; nn++)
            sbuf[p++] = x1flux(k, j, i, nn);
        }
        // x2 direction
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int j = pmb->js + (pmb->je-pmb->js + 1)*(nb.fid & 1);
        if (pmb->block_size.nx3>1) { // 3D
          for (int k=pmb->ks; k<=pmb->ke; k+=2) {
            pco->Face2Area(k  , j, pmb->is, pmb->ie, sarea0);
            pco->Face2Area(k+1, j, pmb->is, pmb->ie, sarea1);
            for (int i=pmb->is; i<=pmb->ie; i+=2) {
              Real tarea = sarea0(i) + sarea0(i+1) + sarea1(i) + sarea1(i+1);
              for (int nn=nl_; nn<=nu_; nn++) {
                sbuf[p++] = (x2flux(k  , j, i  , nn)*sarea0(i  )
                            + x2flux(k  , j, i+1, nn)*sarea0(i+1)
                            + x2flux(k+1, j, i  , nn)*sarea1(i  )
                            + x2flux(k+1, j, i+1, nn)*sarea1(i+1))/tarea;
              }
            }
          }
        } else if (pmb->block_size.nx2>1) { // 2D
          int k = pmb->ks;
          pco->Face2Area(0, j, pmb->is ,pmb->ie, sarea0);
          for (int i=pmb->is; i<=pmb->ie; i+=2) {
            Real tarea = sarea0(i) + sarea0(i+1);
            for (int nn=nl_; nn<=nu_; nn++) {
              sbuf[p++] = (x2flux(k, j, i , nn)*sarea0(i  )
                          + x2flux(k, j, i+1, nn)*sarea0(i+1))/tarea;
            }
          }
        }
        // x3 direction - 3D onl_y
      } else if (nb.fid == BoundaryFace::inner_x3 || nb.fid == BoundaryFace::outer_x3) {
        int k = pmb->ks + (pmb->ke-pmb->ks + 1)*(nb.fid & 1);
        for (int j=pmb->js; j<=pmb->je; j+=2) {
          pco->Face3Area(k, j,   pmb->is, pmb->ie, sarea0);
          pco->Face3Area(k, j+1, pmb->is, pmb->ie, sarea1);
          for (int i=pmb->is; i<=pmb->ie; i+=2) {
            Real tarea = sarea0(i) + sarea0(i+1) + sarea1(i) + sarea1(i+1);
            for (int nn=nl_; nn<=nu_; nn++) {
              sbuf[p++] = (x3flux(k, j  , i, nn)*sarea0(i  )
                           + x3flux( k, j  , i+1, nn)*sarea0(i+1)
                           + x3flux(k, j+1, i  , nn)*sarea1(i  )
                           + x3flux(k, j+1, i+1, nn)*sarea1(i+1))/tarea;
            }
          }
        }
      }
      if (nb.snb.rank == Globals::my_rank) { // on the same node
        CopyFluxCorrectionBufferSameProcess(nb, p);
      }
#ifdef MPI_PARALLEL
      else
        MPI_Start(&(bd_var_flcor_.req_send[nb.bufid]));
#endif
      bd_var_flcor_.sflag[nb.bufid] = BoundaryStatus::completed;
    }
  }
  return;

}


//----------------------------------------------------------------------------------------
//! \fn bool CellCenteredBoundaryVariable::ReceiveFluxCorrection()
//  \brief Receive and apply the surface flux from the finer neighbor(s)

bool RadBoundaryVariable::ReceiveFluxCorrection() {
  MeshBlock *pmb = pmy_block_;
  bool bflag=true;

  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.ni.type != NeighborConnect::face) break;
    if (nb.snb.level == pmb->loc.level+1) {
      if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::completed) continue;
      if (bd_var_flcor_.flag[nb.bufid] == BoundaryStatus::waiting) {
        if (nb.snb.rank == Globals::my_rank) {// on the same process
          bflag = false;
          continue;
        }
#ifdef MPI_PARALLEL
        else { // NOLINT
          int test;
          MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test,
                     MPI_STATUS_IGNORE);
          MPI_Test(&(bd_var_flcor_.req_recv[nb.bufid]), &test, MPI_STATUS_IGNORE);
          if (!static_cast<bool>(test)) {
            bflag = false;
            continue;
          }
          bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::arrived;
        }
#endif
      }
      // boundary arrived; apply flux correction
      int p = 0;
      Real *rbuf=bd_var_flcor_.recv[nb.bufid];
      if (nb.fid == BoundaryFace::inner_x1 || nb.fid == BoundaryFace::outer_x1) {
        int il = pmb->is + (pmb->ie - pmb->is)*nb.fid+nb.fid;
        int jl = pmb->js, ju = pmb->je, kl = pmb->ks, ku = pmb->ke;
        if (nb.ni.fi1 == 0) ju -= pmb->block_size.nx2/2;
        else          jl += pmb->block_size.nx2/2;
        if (nb.ni.fi2 == 0) ku -= pmb->block_size.nx3/2;
        else          kl += pmb->block_size.nx3/2;
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju; j++){
            for (int nn=nl_; nn<=nu_; nn++) 
              x1flux(k,j,il,nn) = rbuf[p++];    
          }
        }
      } else if (nb.fid == BoundaryFace::inner_x2 || nb.fid == BoundaryFace::outer_x2) {
        int jl = pmb->js + (pmb->je - pmb->js)*(nb.fid & 1) + (nb.fid & 1);
        int il = pmb->is, iu = pmb->ie, kl = pmb->ks, ku = pmb->ke;
        if (nb.ni.fi1 == 0) iu -= pmb->block_size.nx1/2;
        else          il += pmb->block_size.nx1/2;
        if (nb.ni.fi2 == 0) ku -= pmb->block_size.nx3/2;
        else          kl += pmb->block_size.nx3/2;
        for (int k=kl; k<=ku; k++) {
          for (int i=il; i<=iu; i++){
            for (int nn=nl_; nn<=nu_; nn++) 
              x2flux(k,jl,i,nn) = rbuf[p++];
          
          }
        }
      } else if (nb.fid == BoundaryFace::inner_x3 || nb.fid == BoundaryFace::outer_x3) {
        int kl = pmb->ks + (pmb->ke - pmb->ks)*(nb.fid & 1) + (nb.fid & 1);
        int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je;
        if (nb.ni.fi1 == 0) iu -= pmb->block_size.nx1/2;
        else          il += pmb->block_size.nx1/2;
        if (nb.ni.fi2 == 0) ju -= pmb->block_size.nx2/2;
        else          jl += pmb->block_size.nx2/2;
        for (int j=jl; j<=ju; j++) {
          for (int i=il; i<=iu; i++){
            for (int nn=nl_; nn<=nu_; nn++) 
              x3flux(kl,j,i,nn) = rbuf[p++];            
          }
        }
      }
      bd_var_flcor_.flag[nb.bufid] = BoundaryStatus::completed;
    }
  }
  return bflag;
}


//----------------------------------------------------------------------------------------
//! \fn int RadBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the same level

int RadBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                              const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;

  si = (nb.ni.ox1 > 0) ? (pmb->ie - NGHOST + 1) : pmb->is;
  ei = (nb.ni.ox1 < 0) ? (pmb->is + NGHOST - 1) : pmb->ie;
  sj = (nb.ni.ox2 > 0) ? (pmb->je - NGHOST + 1) : pmb->js;
  ej = (nb.ni.ox2 < 0) ? (pmb->js + NGHOST - 1) : pmb->je;
  sk = (nb.ni.ox3 > 0) ? (pmb->ke - NGHOST + 1) : pmb->ks;
  ek = (nb.ni.ox3 < 0) ? (pmb->ks + NGHOST - 1) : pmb->ke;
  int p = 0;
  AthenaArray<Real> &var = *var_cc;
  BufferUtility::PackData(var, buf, sk, ek, nl_, nu_, si, ei, sj, ej, p);

  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int RadBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the coarser level

int RadBoundaryVariable::LoadBoundaryBufferToCoarser(Real *buf,
                                                              const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  MeshRefinement *pmr = pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cn = NGHOST - 1;
  AthenaArray<Real> &var = *var_cc;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  si = (nb.ni.ox1 > 0) ? (pmb->cie - cn) : pmb->cis;
  ei = (nb.ni.ox1 < 0) ? (pmb->cis + cn) : pmb->cie;
  sj = (nb.ni.ox2 > 0) ? (pmb->cje - cn) : pmb->cjs;
  ej = (nb.ni.ox2 < 0) ? (pmb->cjs + cn) : pmb->cje;
  sk = (nb.ni.ox3 > 0) ? (pmb->cke - cn) : pmb->cks;
  ek = (nb.ni.ox3 < 0) ? (pmb->cks + cn) : pmb->cke;

  int p = 0;
  // function overload to do restriction for radiation variabbles
  pmr->RestrictCellCenteredValues(var, coarse_var, -1, nl_, nu_, si, ei, sj, ej, sk, ek);
  BufferUtility::PackData(coarse_var, buf, sk, ek, nl_, nu_, si, ei, sj, ej, p);
  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int CellCenteredBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered boundary buffers for sending to a block on the finer level

int RadBoundaryVariable::LoadBoundaryBufferToFiner(Real *buf,
                                                            const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cn = pmb->cnghost - 1;
  AthenaArray<Real> &var = *var_cc;

  si = (nb.ni.ox1 > 0) ? (pmb->ie - cn) : pmb->is;
  ei = (nb.ni.ox1 < 0) ? (pmb->is + cn) : pmb->ie;
  sj = (nb.ni.ox2 > 0) ? (pmb->je - cn) : pmb->js;
  ej = (nb.ni.ox2 < 0) ? (pmb->js + cn) : pmb->je;
  sk = (nb.ni.ox3 > 0) ? (pmb->ke - cn) : pmb->ks;
  ek = (nb.ni.ox3 < 0) ? (pmb->ks + cn) : pmb->ke;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1)   si += pmb->block_size.nx1/2 - pmb->cnghost;
    else            ei -= pmb->block_size.nx1/2 - pmb->cnghost;
  }
  if (nb.ni.ox2 == 0 && pmb->block_size.nx2 > 1) {
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) sj += pmb->block_size.nx2/2 - pmb->cnghost;
      else          ej -= pmb->block_size.nx2/2 - pmb->cnghost;
    } else {
      if (nb.ni.fi2 == 1) sj += pmb->block_size.nx2/2 - pmb->cnghost;
      else          ej -= pmb->block_size.nx2/2 - pmb->cnghost;
    }
  }
  if (nb.ni.ox3 == 0 && pmb->block_size.nx3 > 1) {
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) sk += pmb->block_size.nx3/2 - pmb->cnghost;
      else          ek -= pmb->block_size.nx3/2 - pmb->cnghost;
    } else {
      if (nb.ni.fi2 == 1) sk += pmb->block_size.nx3/2 - pmb->cnghost;
      else          ek -= pmb->block_size.nx3/2 - pmb->cnghost;
    }
  }

  int p = 0;
  BufferUtility::PackData(var, buf, sk, ek, nl_, nu_, si, ei, sj, ej, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryVariable::SetBoundaries()
//  \brief set the boundary data

void RadBoundaryVariable::SetBoundaries() {
  MeshBlock *pmb = pmy_block_;
  int mylevel = pmb->loc.level;
  for (int n=0; n<pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.snb.level == mylevel)
      SetBoundarySameLevel(bd_var_.recv[nb.bufid], nb);
    else if (nb.snb.level < mylevel) // only sets the prolongation buffer
      SetBoundaryFromCoarser(bd_var_.recv[nb.bufid], nb);
    else
      SetBoundaryFromFiner(bd_var_.recv[nb.bufid], nb);
    bd_var_.flag[nb.bufid] = BoundaryStatus::completed; // completed
  }

  if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar ||
      pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar)
    PolarBoundarySingleAzimuthalBlock();

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::SetBoundarySameLevel(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set cell-centered boundary received from a block on the same level

void RadBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                        const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> &var = *var_cc;

  if (nb.ni.ox1 == 0)     si = pmb->is,        ei = pmb->ie;
  else if (nb.ni.ox1 > 0) si = pmb->ie + 1,      ei = pmb->ie + NGHOST;
  else              si = pmb->is - NGHOST, ei = pmb->is - 1;
  if (nb.ni.ox2 == 0)     sj = pmb->js,        ej = pmb->je;
  else if (nb.ni.ox2 > 0) sj = pmb->je + 1,      ej = pmb->je + NGHOST;
  else              sj = pmb->js - NGHOST, ej = pmb->js - 1;
  if (nb.ni.ox3 == 0)     sk = pmb->ks,        ek = pmb->ke;
  else if (nb.ni.ox3 > 0) sk = pmb->ke + 1,      ek = pmb->ke + NGHOST;
  else              sk = pmb->ks - NGHOST, ek = pmb->ks - 1;

  int p = 0;
  // no need to flip for radiation
  if (nb.polar) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
        for (int i=si; i<=ei; ++i) {
 #pragma omp simd linear(p)
          for (int n=nl_; n<=nu_; ++n) {
 //           Real sign = 1.0;
 //           if (flip_across_pole_ != nullptr) 
 //           	sign = flip_across_pole_[n] ? -1.0 : 1.0;
            var(k,j,i,n) = buf[p++];
          }// nu
        }//i
      }//j
    }// k
  } else {
    BufferUtility::UnpackData(buf, var, sk, ek, nl_, nu_, si, ei, sj, ej, p);
  }

  return;
}



//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
//                                                                const NeighborBlock& nb)
//  \brief Set cell-centered prolongation buffer received from a block on a coarser level

void RadBoundaryVariable::SetBoundaryFromCoarser(Real *buf,
                                                          const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cng = pmb->cnghost;
  AthenaArray<Real> &coarse_var = *coarse_buf;

  if (nb.ni.ox1 == 0) {
    si = pmb->cis, ei = pmb->cie;
    if ((pmb->loc.lx1 & 1LL) == 0LL) ei += cng;
    else                             si -= cng;
  } else if (nb.ni.ox1 > 0)  {
    si = pmb->cie + 1,   ei = pmb->cie + cng;
  } else {
    si = pmb->cis - cng, ei = pmb->cis - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = pmb->cjs, ej = pmb->cje;
    if (pmb->block_size.nx2 > 1) {
      if ((pmb->loc.lx2 & 1LL) == 0LL) ej += cng;
      else                             sj -= cng;
    }
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->cje + 1,   ej = pmb->cje + cng;
  } else {
    sj = pmb->cjs - cng, ej = pmb->cjs - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = pmb->cks, ek = pmb->cke;
    if (pmb->block_size.nx3 > 1) {
      if ((pmb->loc.lx3 & 1LL) == 0LL) ek += cng;
      else                             sk -= cng;
    }
  } else if (nb.ni.ox3 > 0)  {
    sk = pmb->cke + 1,   ek = pmb->cke + cng;
  } else {
    sk = pmb->cks - cng, ek = pmb->cks - 1;
  }

  int p = 0;
  if (nb.polar) {
//      Real sign = 1.0;
//      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
        for (int i=si; i<=ei; ++i){
#pragma omp simd linear(p)
          for (int n=nl_; n<=nu_; ++n) {
            coarse_var(k,j,i,n) = buf[p++];
          }
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, coarse_var, sk, ek, nl_, nu_, si, ei, sj, ej, p);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::SetBoundaryFromFiner(Real *buf,
//                                                              const NeighborBlock& nb)
//  \brief Set cell-centered boundary received from a block on a finer level

void RadBoundaryVariable::SetBoundaryFromFiner(Real *buf,
                                                        const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_cc;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;

  if (nb.ni.ox1 == 0) {
    si = pmb->is, ei = pmb->ie;
    if (nb.ni.fi1 == 1)   si += pmb->block_size.nx1/2;
    else            ei -= pmb->block_size.nx1/2;
  } else if (nb.ni.ox1 > 0) {
    si = pmb->ie + 1,      ei = pmb->ie + NGHOST;
  } else {
    si = pmb->is - NGHOST, ei = pmb->is - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = pmb->js, ej = pmb->je;
    if (pmb->block_size.nx2 > 1) {
      if (nb.ni.ox1 != 0) {
        if (nb.ni.fi1 == 1) sj += pmb->block_size.nx2/2;
        else          ej -= pmb->block_size.nx2/2;
      } else {
        if (nb.ni.fi2 == 1) sj += pmb->block_size.nx2/2;
        else          ej -= pmb->block_size.nx2/2;
      }
    }
  } else if (nb.ni.ox2 > 0) {
    sj = pmb->je + 1,      ej = pmb->je + NGHOST;
  } else {
    sj = pmb->js - NGHOST, ej = pmb->js - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = pmb->ks, ek = pmb->ke;
    if (pmb->block_size.nx3 > 1) {
      if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
        if (nb.ni.fi1 == 1) sk += pmb->block_size.nx3/2;
        else          ek -= pmb->block_size.nx3/2;
      } else {
        if (nb.ni.fi2 == 1) sk += pmb->block_size.nx3/2;
        else          ek -= pmb->block_size.nx3/2;
      }
    }
  } else if (nb.ni.ox3 > 0) {
    sk = pmb->ke + 1,      ek = pmb->ke + NGHOST;
  } else {
    sk = pmb->ks - NGHOST, ek = pmb->ks - 1;
  }

  int p = 0;
  if (nb.polar) {
//     Real sign=1.0;
//      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
    for (int k=sk; k<=ek; ++k) {
      for (int j=ej; j>=sj; --j) {
        for (int i=si; i<=ei; ++i){
#pragma omp simd linear(p)
          for (int n=nl_; n<=nu_; ++n) {
            var(k,j,i,n) = buf[p++];
          }//n
        }//i
      }//j
    }// k
  } else {
    BufferUtility::UnpackData(buf, var,  sk, ek, nl_, nu_, si, ei, sj, ej, p);
  }
  return;
}




//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::PolarBoundarySingleAzimuthalBlock()
// \brief polar boundary edge-case: single MeshBlock spans the entire azimuthal (x3) range

void RadBoundaryVariable::PolarBoundarySingleAzimuthalBlock() {
  MeshBlock *pmb = pmy_block_;

  if (pmb->loc.level  ==  pmy_mesh_->root_level && pmy_mesh_->nrbx3 == 1
      && pmb->block_size.nx3 > 1) {
    AthenaArray<Real> &var = *var_cc;
    if (pbval_->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k){
            for (int n=nl_; n<=nu_; ++n) {
              azimuthal_shift_rad_(k,n) = var(k,j,i,n);
            }
          }
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half + NGHOST) ? 1 : -1) * nx3_half;
            for (int n=nl_; n<=nu_; ++n) {
              var(k,j,i,n) = azimuthal_shift_rad_(k_shift,n);
            }//n
          }//k
        }// end i
      }// end j
    }// end inner_x2 polar

    if (pbval_->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) {
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
        for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i) {
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k){
            for (int n=nl_; n<=nu_; ++n) {
              azimuthal_shift_rad_(k,n) = var(k,j,i,n);
            }// end n
          }// end k
          for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
            int k_shift = k;
            k_shift += (k < (nx3_half + NGHOST) ? 1 : -1) * nx3_half;
            for(int n=nl_; n<=nu_; ++n){
              var(k,j,i,n) = azimuthal_shift_rad_(k_shift,n);
            }// end n
          }// end k
          
        }// end i
      }// end j
    }//end outer_x2 polar
  }// end if nx3 > 1
  return;
}




