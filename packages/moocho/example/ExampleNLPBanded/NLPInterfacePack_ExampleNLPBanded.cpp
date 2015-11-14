// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <math.h>

#include "NLPInterfacePack_ExampleNLPBanded.hpp"
#include "DenseLinAlgPack_PermVecMat.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_Assert.hpp"

namespace NLPInterfacePack {

// Constructors / initializers

ExampleNLPBanded::ExampleNLPBanded(
  size_type     nD
  ,size_type    nI
  ,size_type    bw
  ,size_type    mU
  ,size_type    mI
  ,value_type   xo
  ,value_type   xDl
  ,value_type   xDu
  ,value_type   xIl
  ,value_type   xIu
  ,value_type   hl
  ,value_type   hu
  ,bool         nlp_selects_basis
  ,value_type   diag_scal
  ,value_type   diag_vary
  ,bool         sym_basis
  ,value_type   f_offset
  ,value_type   co
  ,bool         ignore_constraints
  )
  :is_initialized_(false)
  ,nlp_selects_basis_(nlp_selects_basis)
  ,basis_selection_was_given_(false)
  ,has_var_bounds_(false)
  ,f_offset_(f_offset)
  ,nD_(nD)
  ,nI_(nI)
  ,bw_(bw)
  ,mU_(mU)
  ,mI_(mI)
  ,ignore_constraints_(ignore_constraints)
  ,diag_scal_(diag_scal)
  ,diag_vary_(diag_vary)
  ,fu_( sym_basis ? 3 : 6 )
{
#ifdef TEUCHOS_DEBUG	
  const char msg_err_head[] = "ExampleNLPBanded::ExampleNLPBanded(...) : Error";
  TEUCHOS_TEST_FOR_EXCEPTION(
    nD <= 0, std::invalid_argument
    ,msg_err_head<<"!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    nI <= 0 || nD < nI, std::invalid_argument
    ,msg_err_head<<"!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    bw < 1 || nD < bw, std::invalid_argument
    ,msg_err_head<<"!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    mU < 0, std::invalid_argument
    ,msg_err_head<<"!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    mI < 0, std::invalid_argument
    ,msg_err_head<<"!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    mU != 0, std::invalid_argument
    ,msg_err_head<<", can't handle undecomposed equalities yet!" );
#endif
  Gc_orig_nz_ = nD_ * ( 1 + 2*(bw_-1) + nI_ ); // Overestimate, ToDo: compute the exact value!
  Gh_orig_nz_ = 2*mI_;
  //
  xinit_orig_.resize(nD_ + nI_);
  xl_orig_.resize(xinit_orig_.dim());
  xu_orig_.resize(xinit_orig_.dim());
  hl_orig_.resize(mI_);
  hu_orig_.resize(mI_);
  co_orig_.resize(nD_ + mU_);
  //
  xinit_orig_             = xo;
  xl_orig_(1,nD_)         = xDl;
  xu_orig_(1,nD_)         = xDu;
  xl_orig_(nD_+1,nD_+nI_) = xIl;
  xu_orig_(nD_+1,nD_+nI_) = xIu;
  if( mI_ ) {
    hl_orig_    = hl;
    hu_orig_    = hu;
  }
  co_orig_ = co;
  //
  const value_type inf = NLP::infinite_bound();
  if( xDl > -inf || xDu < +inf || xIl > -inf || xIu < +inf )
    has_var_bounds_ = true;
  else
    has_var_bounds_ = false;
}

// Overridden public members from NLP

void ExampleNLPBanded::initialize(bool test_setup)
{
  if(is_initialized_) {
    NLPSerialPreprocessExplJac::initialize(test_setup);
    return;
  }
  // Nothing to initialize?
  NLPSerialPreprocessExplJac::initialize(test_setup);
  is_initialized_ = true;
}

bool ExampleNLPBanded::is_initialized() const
{
  return is_initialized_;
}

value_type ExampleNLPBanded::max_var_bounds_viol() const
{
  return +1e+20; // Functions defined everywhere!
}

// Overridden protected methods from NLPSerialPreprocess

bool ExampleNLPBanded::imp_nlp_has_changed() const
{
  return !is_initialized_;
}

size_type ExampleNLPBanded::imp_n_orig() const
{
  return nD_ + nI_;
}

size_type ExampleNLPBanded::imp_m_orig() const
{
  if(ignore_constraints_) return 0;
  return nD_ + mU_;
}

size_type ExampleNLPBanded::imp_mI_orig() const
{
  if(ignore_constraints_) return 0;
  return mI_;
}

const DVectorSlice ExampleNLPBanded::imp_xinit_orig() const
{
  return xinit_orig_();
}

bool ExampleNLPBanded::imp_has_var_bounds() const
{
  return has_var_bounds_;
}

const DVectorSlice ExampleNLPBanded::imp_xl_orig() const
{
  return xl_orig_();
}

const DVectorSlice ExampleNLPBanded::imp_xu_orig() const
{
  return xu_orig_();
}

const DVectorSlice ExampleNLPBanded::imp_hl_orig() const
{
  return hl_orig_();
}

const DVectorSlice ExampleNLPBanded::imp_hu_orig() const
{
  return hu_orig_();
}

void ExampleNLPBanded::imp_calc_f_orig(
  const DVectorSlice            &x_full
  ,bool                        newx
  ,const ZeroOrderInfoSerial   &zero_order_info
  ) const
{
  inform_new_point(newx);
  const DVectorSlice x_orig = x_full(1,imp_n_orig());
  *zero_order_info.f = f_offset_ + ( 1.0 / 2.0 ) * DenseLinAlgPack::dot( x_orig, x_orig );
}

void ExampleNLPBanded::imp_calc_c_orig(
  const DVectorSlice            &x_full
  ,bool                        newx
  ,const ZeroOrderInfoSerial   &zero_order_info
  ) const
{
  inform_new_point(newx);
  if(c_orig_updated_)
    return; // c(x) is already computed in *zero_order_info.c
  TEUCHOS_TEST_FOR_EXCEPT( !( zero_order_info.c ) );
  DVector
    &c  = *zero_order_info.c;
  const size_type
    num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
    I_remainder = nD_ % nI_;
  size_type j = 0;
  const value_type
    ds_alpha = nD_ > 1 ? diag_scal_ * (diag_vary_ - 1.0)/(nD_ - 1.0) : 0.0,
    ds_beta = diag_scal_;
  for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
    const size_type num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
    for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
      ++j;
      const size_type
        klu = ( j - bw_     >= 0   ? bw_-1 : j-1   ),
        kuu = ( j + bw_ - 1 <= nD_ ? bw_-1 : nD_-j );
      const value_type
        ds_j = ds_alpha * (j-1) + ds_beta;
      value_type
        &c_j = (c(j) = ds_j * x_full(j));
      {for( size_type k = 1; k <= klu; ++k ) {
        c_j -= (3.0 / k) * x_full(j-k);
      }}
      {for( size_type k = 1; k <= kuu; ++k ) {
        c_j -= (fu_ / k) * x_full(j+k);
      }}
      const value_type
        term = x_full(nD_ + q_i) + 1;
      c_j *= (term * term);
      c_j += co_orig_(j);
    }
  }
  c_orig_updated_ = true;
}

void ExampleNLPBanded::imp_calc_h_orig(
  const DVectorSlice            &x_full
  ,bool                        newx
  ,const ZeroOrderInfoSerial   &zero_order_info
  ) const
{
  inform_new_point(newx);
  DVector
    &h  = *zero_order_info.h;
  const size_type
    num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
    I_remainder = nD_ % nI_;
  size_type jI = 0;
  for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
    const value_type  x_q = x_full(nD_ + q_i);
    const size_type   num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
    for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
      ++jI;
      if( jI > mI_ ) return;
      h(jI) = x_full(jI) - x_q;
    }
  }
 }

void ExampleNLPBanded::imp_calc_Gf_orig(
  const DVectorSlice            &x_full
  ,bool                        newx
  ,const ObjGradInfoSerial     &obj_grad_info
  ) const
{
  inform_new_point(newx);
  const Range1D var_orig(1,imp_n_orig());
  (*obj_grad_info.Gf)(var_orig) = x_full(var_orig);
}

bool ExampleNLPBanded::imp_get_next_basis(
  IVector      *var_perm_full
  ,IVector     *equ_perm_full
  ,size_type   *rank_full
  ,size_type   *rank
  )
{
  if(basis_selection_was_given_)
    return false; // Already gave this basis selection.
  // Select the first nD variables as basis variables which gives
  // a nice banded matrix for the basis matrix C.
  // Also, if the general inequality constraints are begin
  // converted to equalities with slacks, make the slack variables
  // basic variables also (after the nD variables).
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT( !( var_perm_full ) );
  TEUCHOS_TEST_FOR_EXCEPT( !( equ_perm_full ) );
  TEUCHOS_TEST_FOR_EXCEPT( !( rank ) );
#endif
  const size_type    n_orig = nD_    + nI_;
  const size_type    n_full = n_orig + mI_;
  const size_type    m_full = nD_    + mI_;
  var_perm_full->resize(n_full);
  equ_perm_full->resize(m_full);
  if( mI_ ) {
    index_type k, i_perm = 1;;
    // basic variables
    for( k = 1; k <= nD_; ++k, ++i_perm )  // Put xD variables first
      (*var_perm_full)(i_perm) = k;
    for( k = 1; k <= mI_; ++k, ++i_perm )  // Put slacks after xD
      (*var_perm_full)(i_perm) = n_orig + k;
    // non-basic variables
    for( k = 1; k <= nI_; ++k, ++i_perm )  // Followed by nI
      (*var_perm_full)(i_perm) = nD_ + k;
  }
  else {
    DenseLinAlgPack::identity_perm( var_perm_full );
  }
  DenseLinAlgPack::identity_perm( equ_perm_full );
  // Count the number of fixed basic variables to reduce
  // the rank of the basis.
  index_type num_fixed = 0;
  for( index_type k = 1; k <= nD_; ++k ) {
    if( xl_orig_(k) == xu_orig_(k) )
      ++num_fixed;
  }
  if( mI_ ) {
    for( index_type k = 1; k <= mI_; ++k ) {
      if( hl_orig_(k) == hu_orig_(k) )
        ++num_fixed;
    }
  }
  // Set the rank
  *rank_full = m_full;
  *rank      = m_full - num_fixed;
  basis_selection_was_given_ = true;
  return true;
}

void ExampleNLPBanded::imp_report_orig_final_solution(
  const DVectorSlice      &x_orig
  ,const DVectorSlice     *lambda_orig
  ,const DVectorSlice     *lambdaI_orig
  ,const DVectorSlice     *nu_orig
  ,bool                   is_optimal
  )
{
  // ToDo: Do something with the final soltuion?
}

bool ExampleNLPBanded::nlp_selects_basis() const
{
  return nlp_selects_basis_;
}

// Overridden protected methods from NLPSerialPreprocessExplJac

size_type ExampleNLPBanded::imp_Gc_nz_orig() const
{
  return Gc_orig_nz_;
}

size_type ExampleNLPBanded::imp_Gh_nz_orig() const
{
  return Gh_orig_nz_;
}

void ExampleNLPBanded::imp_calc_Gc_orig(
  const DVectorSlice& x_full, bool newx
  , const FirstOrderExplInfo& first_order_expl_info
  ) const
{
  inform_new_point(newx);
  // Compute c(x) if not already (will compute in temp if needed)
  this->imp_calc_c_orig( x_full, newx, zero_order_orig_info() );
  TEUCHOS_TEST_FOR_EXCEPT( !( first_order_expl_info.c ) );
  DVector
    &c = *first_order_expl_info.c;
  // Get references/pointers to data for Gc to be computed/updated.
  index_type
    &Gc_nz = *first_order_expl_info.Gc_nz;
  value_type
    *Gc_val = &(*first_order_expl_info.Gc_val)[0];
  index_type
    *Gc_ivect = ( first_order_expl_info.Gc_ivect
            ? &(*first_order_expl_info.Gc_ivect)[0] : NULL ),
    *Gc_jvect = ( first_order_expl_info.Gc_jvect
            ? &(*first_order_expl_info.Gc_jvect)[0] : NULL );
  TEUCHOS_TEST_FOR_EXCEPT( !(  (Gc_ivect != NULL) == (Gc_jvect != NULL)  ) );
  // Set nonzeros for Gc (in sorted compressed column format w.r.t., i.e. grouped by constraints)
  const size_type
    num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
    I_remainder = nD_ % nI_;
  Gc_nz = 0;
  size_type j = 0;
  const value_type
    ds_alpha = nD_ > 1 ? diag_scal_ * (diag_vary_ - 1.0)/(nD_ - 1.0) : 0.0,
    ds_beta = diag_scal_;
  for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
    const value_type
      x_q = x_full(nD_ + q_i),
      x_q_term = (x_q + 1) * (x_q + 1);
    const size_type num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
    for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
      ++j;
      const size_type
        klu = ( j - bw_     >= 0   ? bw_-1 : j-1   ),
        kuu = ( j + bw_ - 1 <= nD_ ? bw_-1 : nD_-j );
      const value_type
        ds_j = ds_alpha * (j-1) + ds_beta;
      //
      {for( index_type k = klu; k >= 1; --k ) {
        ++Gc_nz;
        *Gc_val++ = -3.0 / k * x_q_term;
        if(Gc_ivect) {
          *Gc_ivect++ = j - k;
          *Gc_jvect++ = j;
        }
      }}
      //
      ++Gc_nz;
      *Gc_val++ = ds_j * x_q_term;
      if(Gc_ivect) {
        *Gc_ivect++ = j;
        *Gc_jvect++ = j;
      }
      //
      {for( index_type k = 1; k <= kuu; ++k ) {
        ++Gc_nz;
        *Gc_val++ = -fu_ / k * x_q_term;
        if(Gc_ivect) {
          *Gc_ivect++ = j + k;
          *Gc_jvect++ = j;
        }
      }}
      //
      ++Gc_nz;
      *Gc_val++ = 2.0 * (c(j) - co_orig_(j)) / (x_q + 1);
      if(Gc_ivect) {
        *Gc_ivect++ = nD_ + q_i;
        *Gc_jvect++ = j;
      }
    }
  }
}

void ExampleNLPBanded::imp_calc_Gh_orig(
  const DVectorSlice& x_full, bool newx
  , const FirstOrderExplInfo& first_order_expl_info
  ) const
{
  inform_new_point(newx);
  // Get references/pointers to data for Gh to be computed/updated.
  index_type
    &Gh_nz = *first_order_expl_info.Gh_nz;
  value_type
    *Gh_val = &(*first_order_expl_info.Gh_val)[0];
  index_type
    *Gh_ivect = ( first_order_expl_info.Gh_ivect
            ? &(*first_order_expl_info.Gh_ivect)[0] : NULL ),
    *Gh_jvect = ( first_order_expl_info.Gh_jvect
            ? &(*first_order_expl_info.Gh_jvect)[0] : NULL );
  TEUCHOS_TEST_FOR_EXCEPT( !(  (Gh_ivect != NULL) == (Gh_jvect != NULL)  ) );
  // Set nonzeros for Gh (in sorted compressed column format w.r.t., i.e. grouped by constraints)
  const size_type
    num_I_per_D = nD_ / nI_,  // Integer division (rounds down)
    I_remainder = nD_ % nI_;
  Gh_nz = 0;
  size_type jI = 0;
  for( size_type q_i = 1; q_i <= nI_; ++q_i ) {
    const size_type   nD_q_i = nD_ + q_i;
    const size_type   num_I_per_D_local = num_I_per_D + ( q_i <= I_remainder ? 1 : 0 );
    for( size_type q_k = 0; q_k < num_I_per_D_local; ++q_k ) {
      ++jI;
      if( jI > mI_ ) return;
      // w.r.t. x(jI)
      ++Gh_nz;
      *Gh_val++ = 1.0;
      if(Gh_ivect) {
        *Gh_ivect++ = jI;
        *Gh_jvect++ = jI;
      }
      // w.r.t. x(nD+q(jI))
      ++Gh_nz;
      *Gh_val++ = -1.0;
      if(Gh_ivect) {
        *Gh_ivect++ = nD_q_i;
        *Gh_jvect++ = jI;
      }
    }
  }
}

// private

void ExampleNLPBanded::assert_is_initialized() const
{
  TEUCHOS_TEST_FOR_EXCEPT(true); //  ToDo: Implemenet!
}

void ExampleNLPBanded::inform_new_point(bool newx) const
{
  if(newx) {
    c_orig_updated_ = false;
  }
}

}	// end namespace NLPInterfacePack
