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

#include <assert.h>

#include <iostream> // Debug only
#include "DenseLinAlgPack_PermOut.hpp"

#include <algorithm>
#include <sstream>
#include <limits>
#include <stdio.h>
#include <fstream>

#include "NLPInterfacePack_NLPSerialPreprocess.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack_PermutationSerial.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_IVector.hpp"
#include "DenseLinAlgPack_PermVecMat.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace LinAlgOpPack {
  using AbstractLinAlgPack::Vp_StV;
}

namespace NLPInterfacePack {

// NLPSerialPreprocess

// Static public members

value_type
NLPSerialPreprocess::fixed_var_mult()
{
  return std::numeric_limits<DenseLinAlgPack::DVector::value_type>::max()-100; // Don't know what to use?
}

// Constructors / nitializers

NLPSerialPreprocess::NLPSerialPreprocess(
  )
  :initialized_(false)
  ,force_xinit_in_bounds_(true)
    ,scale_f_(1.0)
  ,basis_selection_num_(0)

{}

// Overridden public members from NLP

void NLPSerialPreprocess::force_xinit_in_bounds(bool force_xinit_in_bounds)
{
  force_xinit_in_bounds_ = force_xinit_in_bounds;
}

bool NLPSerialPreprocess::force_xinit_in_bounds() const
{
  return force_xinit_in_bounds_;
}

void NLPSerialPreprocess::initialize(bool test_setup)
{
  namespace mmp = MemMngPack;

  const value_type inf_bnd = NLP::infinite_bound();

  basis_selection_num_ = 0;

  if( initialized_  && !imp_nlp_has_changed() ) {
    // The subclass NLP has not changed so we can just
    // slip this preprocessing.
    NLPObjGrad::initialize(test_setup);
    return;
  }

  // Get the dimensions of the original problem

  n_orig_  = imp_n_orig();
  m_orig_  = imp_m_orig();   // This may be zero!
  mI_orig_ = imp_mI_orig();  // This may be zero!
  
  // Get the dimensions of the full problem

  n_full_  = n_orig_ + mI_orig_;
  m_full_  = m_orig_ + mI_orig_;

  // Initialize the storage for the intermediate quanities
  
  xinit_full_.resize(n_full_);
  xl_full_.resize(n_full_);
  xu_full_.resize(n_full_);
  x_full_.resize(n_full_);
  c_orig_.resize(m_orig_);
  h_orig_.resize(mI_orig_);
  Gf_full_.resize(n_full_);
  var_full_to_fixed_.resize(n_full_);
  equ_perm_.resize(m_full_);
  inv_equ_perm_.resize(m_full_);
  space_c_.initialize(m_full_);
  space_c_breve_.initialize(m_orig_);
  space_h_breve_.initialize(mI_orig_);
  factory_P_var_   = Teuchos::rcp( new Teuchos::AbstractFactoryStd<Permutation,PermutationSerial>() );
  factory_P_equ_   = Teuchos::rcp( new Teuchos::AbstractFactoryStd<Permutation,PermutationSerial>() );

  // Intialize xinit_full_, xl_full_ and xu_full_ for the initial point which will set the
  // fixed elements which will not change during the optimization.
  xinit_full_(1,n_orig_)  = imp_xinit_orig();
  xl_full_(1,n_orig_)     = imp_xl_orig();
  xu_full_(1,n_orig_)     = imp_xu_orig();
  if( n_full_ > n_orig_ ) { // Include slack varaibles
    xinit_full_(n_orig_+1,n_full_)  = 0.0;
    xl_full_(n_orig_+1,n_full_)     = imp_hl_orig();
    xu_full_(n_orig_+1,n_full_)     = imp_hu_orig();
  }

  const bool has_var_bounds = imp_has_var_bounds() || n_full_ > n_orig_;

  // Force the initial point in bounds if it is not.
  if( force_xinit_in_bounds() && has_var_bounds ) {
    AbstractLinAlgPack::force_in_bounds(
      VectorMutableDense( xl_full_(), Teuchos::null )
      ,VectorMutableDense( xu_full_(), Teuchos::null )
      ,&VectorMutableDense( x_full_(), Teuchos::null )
      );
  }
  
  // Determine which variables are fixed by bounds!
  size_type
    xl_nz     = 0,
    xu_nz     = 0,
    num_bnd_x = 0;
  if( has_var_bounds ) {
    // Determine which variables are fixed by bounds and
    // adjust the bounds if needed.
    DVector::iterator
      xl_full		= xl_full_.begin(),
      xu_full		= xu_full_.begin();
    n_ = 0;
    size_type num_fixed = 0;
    for(int i = 1; i <= n_full_; ++i, ++xl_full, ++xu_full) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        *xl_full > *xu_full, InconsistantBounds
        ,"NLPSerialPreprocess::initialize() : Error, Inconsistant bounds: xl_full("
        << i << ") > xu_full(" << i << ")" ); 
      if(*xl_full == *xu_full) {
        //
        // Fixed between bounds
        //
        var_full_to_fixed_(n_full_ - num_fixed) = i;
        num_fixed++;
      }
      else {
        //
        // Not Fixed between bounds
        //
        // Adjust the bounds if needed
        *xl_full = *xl_full < -inf_bnd ? -inf_bnd : *xl_full;
        *xu_full = *xu_full > +inf_bnd ? +inf_bnd : *xu_full;
        //
        n_++;
        var_full_to_fixed_(n_) = i;
        // Check if xl is bounded
        if(	*xl_full != -inf_bnd )
          ++xl_nz;
        // Check if xu is bounded
        if(	*xu_full != inf_bnd )
          ++xu_nz;
        if( *xl_full != -inf_bnd || *xu_full != inf_bnd )
          ++num_bnd_x;
      }
    }
  }
  else {
    // None of the variables are fixed by bounds because there are no bounds
    n_ = n_full_;
    DenseLinAlgPack::identity_perm( &var_full_to_fixed_ );
  }

//	std::cerr << "n_ =" << n_ << std::endl;
//	std::cerr << "var_full_to_fixed_ =\n" << var_full_to_fixed_;
  
  num_bounded_x_ = num_bnd_x;

  // Validate that we still have a solvable problem
  TEUCHOS_TEST_FOR_EXCEPTION(
    n_ < m_full_, InvalidInitialization
    ,"NLPSerialPreprocess::initialize() : Error, after removing fixed "
    << "variables, n = " << n_ << " < m = " << m_full_
    << ", and the NLP is over determined and can not be solved!" );

  // Initialize inverse of var_full_to_fixed_
  DenseLinAlgPack::inv_perm( var_full_to_fixed_, &inv_var_full_to_fixed_ );

//	std::cerr << "inv_var_full_to_fixed_ =\n"  << inv_var_full_to_fixed_;

  var_perm_.resize(n_);
  space_x_.initialize(n_);
  
  // Resize xinit, xl, xu, hl and hu
  xinit_.initialize(n_);
  xl_.initialize(n_);
  xu_.initialize(n_);
  if(mI_orig_) {
    hl_breve_.initialize(mI_orig_);
    hu_breve_.initialize(mI_orig_);
  }

  if( m_full_ ) {
    // Get the first basis
    if( !nlp_selects_basis() ) {
      // The NLP is not selecting the first basis so set to the initial basis to
      // the indentity permutations and assume full column rank for Gc.
      DenseLinAlgPack::identity_perm(&var_perm_);
//			std::cerr << "var_perm_ =\n" << var_perm_;
      DenseLinAlgPack::identity_perm(&equ_perm_);
//			std::cerr << "equ_perm_ =\n" << equ_perm_;
      DenseLinAlgPack::identity_perm(&inv_equ_perm_);
//			std::cerr << "inv_equ_perm_ =\n" << inv_equ_perm_;
      r_ = m_full_;
      var_from_full( xinit_full_().begin(), xinit_.set_vec().begin() );
      if(has_var_bounds) {
        var_from_full( xl_full_().begin(), xl_.set_vec().begin() );
        var_from_full( xu_full_().begin(), xu_.set_vec().begin() );
        do_force_xinit_in_bounds();
      }
      else {
        xl_ = -inf_bnd;
        xu_ = +inf_bnd;
      }
    }
    else {
      // The nlp subclass is selecting the first basis.
      
      // make intialized_ true temporaraly so you can call get_next_basis()
      // and assert_initialized() called in it will not throw an exception.
      initialized_ = true;
      
      try {
        size_type rank;
        const bool 
           get_next_basis_return = get_next_basis_remove_fixed( &var_perm_, &equ_perm_, &rank );
        TEUCHOS_TEST_FOR_EXCEPTION(
          !get_next_basis_return, std::logic_error
          ,"NLPSerialPreprocess::initialize():  "
          " If nlp_selects_basis() is true then imp_get_next_basis() "
          " must return true for the first call" );
        assert_and_set_basis( var_perm_, equ_perm_, rank );
//				std::cerr << "var_perm_ =\n" << var_perm_;
//				std::cerr << "equ_perm_ =\n" << equ_perm_;
      }
      catch(...) {
        // In case an exception was thrown I don't want to leave #this#
        // in an inconsistant state.
        initialized_ = false;
        throw;
      }
      
      initialized_ = false;	// resize to false to continue initialization
    }
  }
  else {
    DenseLinAlgPack::identity_perm(&var_perm_);
    r_ = 0;
    var_from_full( xinit_full_().begin(), xinit_.set_vec().begin() );
    if(has_var_bounds) {
      var_from_full( xl_full_().begin(), xl_.set_vec().begin() );
      var_from_full( xu_full_().begin(), xu_.set_vec().begin() );
      do_force_xinit_in_bounds();
    }
    else {
      xl_ = -inf_bnd;
      xu_ = +inf_bnd;
    }
  }
  
//	std::cerr << "n_full_ = " << n_full_ << std::endl;
//	std::cerr << "n_ = " << n_ << std::endl;
//	std::cerr << "var_full_to_fixed_ =\n" << var_full_to_fixed_;
//	std::cerr << "inv_var_full_to_fixed_ =\n"  << inv_var_full_to_fixed_;
//	std::cerr << "var_perm_ =\n" << var_perm_;
//	std::cerr << "equ_perm_ =\n" << equ_perm_;

  // If you get here then the initialization went Ok.
  NLPObjGrad::initialize(test_setup);
  initialized_ = true;
}

bool NLPSerialPreprocess::is_initialized() const
{
  return initialized_;
}

size_type NLPSerialPreprocess::n() const 
{
  assert_initialized();
  return n_;
}

size_type NLPSerialPreprocess::m() const 
{	
  assert_initialized(); 
  return m_full_;
}

NLP::vec_space_ptr_t NLPSerialPreprocess::space_x() const
{
  namespace mmp = MemMngPack;
  return Teuchos::rcp(&space_x_,false);
}

NLP::vec_space_ptr_t NLPSerialPreprocess::space_c() const
{
  namespace mmp = MemMngPack;
  return ( m_full_ ? Teuchos::rcp(&space_c_,false) : Teuchos::null );
}

size_type NLPSerialPreprocess::num_bounded_x() const 
{
  return num_bounded_x_;
}

const Vector& NLPSerialPreprocess::xl() const 
{
  assert_initialized();
  return xl_;
}

const Vector& NLPSerialPreprocess::xu() const 
{
  assert_initialized();
  return xu_;
}

const Vector& NLPSerialPreprocess::xinit() const 
{
  assert_initialized();
  return xinit_;
}

void NLPSerialPreprocess::get_init_lagrange_mult(
  VectorMutable*   lambda
  ,VectorMutable*  nu
  ) const
{
  assert_initialized();
  // ToDo: Get subclass to define what these are!
  if(lambda)
    *lambda   = 0.0;
  if(nu)
    *nu      = 0.0;
}

void NLPSerialPreprocess::scale_f( value_type scale_f )
{
  assert_initialized();
  scale_f_ = scale_f;
}

value_type NLPSerialPreprocess::scale_f() const
{
  assert_initialized();
  return scale_f_;
}

void NLPSerialPreprocess::report_final_solution(
  const Vector&    x
  ,const Vector*   lambda
  ,const Vector*   nu
  ,bool            is_optimal
  )
{
  assert_initialized();
  // set x_full
  VectorDenseEncap  x_d(x);
  DVector x_full(n_full_);
  x_full = x_full_;	// set any fixed components (as well as the others at first)
  var_to_full( x_d().begin(), x_full().begin() );	// set the nonfixed components
  // set lambda_full
  DVector lambda_full;
  if( lambda ) {
    VectorDenseEncap lambda_d(*lambda);
    DVectorSlice      lambda = lambda_d();
    lambda_full.resize(m_full_);
    for(size_type j = 1; j <= m_full_; j++)
      lambda_full(equ_perm_(j)) = lambda(j);
  }
  // set nu_full
  DVector nu_full(n_full_);
  if(nu) {
    nu_full = 0.0; // We don't give lagrange multipliers for fixed varaibles!
    // ToDo: Define a special constrant for multiplier values for fixed variables 
    VectorDenseEncap nu_d(*nu);
    var_to_full( nu_d().begin(), nu_full().begin() );	// set the nonfixed components
  }
  // Report the final solution
  DVectorSlice
    lambda_orig   = lambda && m_orig_ ? lambda_full(1,m_orig_) : DVectorSlice(),
    lambdaI_orig  = ( lambda && m_full_ > m_orig_ 
              ? lambda_full(m_orig_+1,m_full_)
              : DVectorSlice() ),
    nu_orig       = nu ? nu_full(1,n_orig_) : DVectorSlice();
  imp_report_orig_final_solution(
    x_full()
    ,lambda_orig.dim()  ? &lambda_orig  : NULL
    ,lambdaI_orig.dim() ? &lambdaI_orig : NULL
    ,nu_orig.dim()      ? &nu_orig      : NULL
    ,is_optimal
    );
}

size_type NLPSerialPreprocess::ns() const
{
  assert_initialized();
  return mI_orig_;
}

NLP::vec_space_ptr_t
NLPSerialPreprocess::space_c_breve() const
{
  namespace mmp = MemMngPack;
  assert_initialized();
  return ( m_orig_ ? Teuchos::rcp(&space_c_breve_,false) : Teuchos::null );
} 
NLP::vec_space_ptr_t
NLPSerialPreprocess::space_h_breve() const
{
  namespace mmp = MemMngPack;
  assert_initialized();
  return ( mI_orig_ ? Teuchos::rcp(&space_h_breve_,false) : Teuchos::null );
}

const Vector& NLPSerialPreprocess::hl_breve() const
{
  assert_initialized();
  return hl_breve_;
}

const Vector& NLPSerialPreprocess::hu_breve() const
{
  assert_initialized();
  return hu_breve_;
}

const Permutation& NLPSerialPreprocess::P_var() const
{
  assert_initialized();
  return P_var_;
}

const Permutation& NLPSerialPreprocess::P_equ() const
{
  assert_initialized();
  return P_equ_;
}

// Overridden public members from NLPVarReductPerm

const NLPVarReductPerm::perm_fcty_ptr_t
NLPSerialPreprocess::factory_P_var() const
{
  return factory_P_var_;
}

const NLPVarReductPerm::perm_fcty_ptr_t
NLPSerialPreprocess::factory_P_equ() const
{
  return factory_P_equ_;
}

Range1D NLPSerialPreprocess::var_dep() const
{
  assert_initialized();
  return r_ ? Range1D(1,r_) : Range1D::Invalid;
}

Range1D NLPSerialPreprocess::var_indep() const
{
  assert_initialized();
  return Range1D(r_+1,n_);
}

Range1D NLPSerialPreprocess::equ_decomp() const
{
  assert_initialized();
  return r_ ? Range1D(1,r_) : Range1D::Invalid;
}

Range1D NLPSerialPreprocess::equ_undecomp() const
{
  assert_initialized();
  return r_ < m_full_ ? Range1D(r_+1,m_full_) : Range1D::Invalid;
}

bool NLPSerialPreprocess::nlp_selects_basis() const
  {
  // Check if the user has supplied a basis from a file
  char ind[17];
  sprintf(ind, "%d", basis_selection_num_);
  std::string fname = "basis_";
  fname += ind;
  fname += ".sel";

  std::ifstream basis_file(fname.c_str());
  if (basis_file)
    {
    return true;
    }

  return false;
  }

bool NLPSerialPreprocess::get_next_basis(
  Permutation*  P_var,   Range1D* var_dep
  ,Permutation* P_equ,   Range1D* equ_decomp
  )
{
  assert_initialized();
  size_type rank = 0;
  const bool 
    get_next_basis_return = get_next_basis_remove_fixed( &var_perm_, &equ_perm_, &rank );
  if(get_next_basis_return)
    assert_and_set_basis( var_perm_, equ_perm_, rank );
  else
    return false; // The NLP subclass did not have a new basis to give us!
  this->get_basis(P_var,var_dep,P_equ,equ_decomp);
  return true;
}

void NLPSerialPreprocess::set_basis(
  const Permutation   &P_var,   const Range1D  &var_dep
  ,const Permutation  *P_equ,   const Range1D  *equ_decomp
  )
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  TEUCHOS_TEST_FOR_EXCEPTION(
    (m_full_ > 0 && (P_equ == NULL || equ_decomp == NULL))
    ,std::invalid_argument
    ,"NLPSerialPreprocess::set_basis(...) : Error!" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    m_full_ > 0 && var_dep.size() != equ_decomp->size()
    ,InvalidBasis
    ,"NLPSerialPreprocess::set_basis(...) : Error!" );
  // Get the concrete types
  const PermutationSerial
    &P_var_s   = dyn_cast<const PermutationSerial>(P_var),
    *P_equ_s   = m_full_  ? &dyn_cast<const PermutationSerial>(*P_equ)   : NULL;
  // Get the underlying permutation vectors
  Teuchos::RCP<IVector>
    var_perm   = Teuchos::rcp_const_cast<IVector>(P_var_s.perm()),
    equ_perm   = ( m_full_
             ? Teuchos::rcp_const_cast<IVector>(P_equ_s->perm())
             : Teuchos::null );
  TEUCHOS_TEST_FOR_EXCEPTION(
    (m_full_ > 0 && equ_perm.get() == NULL)
    ,std::invalid_argument
    ,"NLPSerialPreprocess::set_basis(...) : Error, P_equ is not initialized properly!" );
  // Set the basis
  assert_and_set_basis( *var_perm, *equ_perm, var_dep.size() );
}

void NLPSerialPreprocess::get_basis(
  Permutation*  P_var,   Range1D* var_dep
  ,Permutation* P_equ,   Range1D* equ_decomp
  ) const
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  assert_initialized();
  TEUCHOS_TEST_FOR_EXCEPTION(
    P_var == NULL || var_dep == NULL
    || (m_full_ > 0 && (P_equ == NULL || equ_decomp == NULL))
    ,std::invalid_argument
    ,"NLPSerialPreprocess::get_basis(...) : Error!" );
  // Get the concrete types
  PermutationSerial
    &P_var_s   = dyn_cast<PermutationSerial>(*P_var),
    *P_equ_s   = m_full_  ? &dyn_cast<PermutationSerial>(*P_equ)   : NULL;
  // Get the underlying permutation vectors
  Teuchos::RCP<IVector>
    var_perm   = Teuchos::rcp_const_cast<IVector>(P_var_s.perm()),
    equ_perm   = ( m_full_
             ? Teuchos::rcp_const_cast<IVector>(P_equ_s->perm())
             : Teuchos::null );
  // Allocate permutation vectors if none allocated yet or someone else has reference to them
  if( var_perm.get() == NULL || var_perm.total_count() > 2 ) // P_var reference and my reference
    var_perm = Teuchos::rcp( new IVector(n_) );
  if( m_full_ && ( equ_perm.get() == NULL || equ_perm.total_count() > 2 ) ) // P_equ reference and my reference
    equ_perm = Teuchos::rcp( new IVector(m_full_) );
  // Copy the basis selection
  (*var_perm)   = var_perm_;
  (*var_dep)    = Range1D(1,r_);
  if(m_full_) {
    (*equ_perm)   = equ_perm_;
    (*equ_decomp) = Range1D(1,r_);
  }
  // Reinitialize the Permutation arguments.
  P_var_s.initialize( var_perm, Teuchos::null, true );  // Allocate the inverse permuation as well!
  if(m_full_)
    P_equ_s->initialize( equ_perm, Teuchos::null, true );
}

// Overridden protected members from NLP

void NLPSerialPreprocess::imp_calc_f(
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    &zero_order_info
  ) const
{
  assert_initialized();
  VectorDenseEncap  x_d(x);
  set_x_full( x_d(), newx, &x_full_() );
  imp_calc_f_orig( x_full(), newx, zero_order_orig_info() );
  *zero_order_info.f = scale_f_ * f_orig_;
}

void NLPSerialPreprocess::imp_calc_c(
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    &zero_order_info
  ) const
{
  assert_initialized();
  VectorDenseEncap  x_d(x);
  set_x_full( x_d(), newx, &x_full_() );
  if( m_orig_ )
    imp_calc_c_orig( x_full(), newx, zero_order_orig_info() );
  if( mI_orig_ )
    imp_calc_h_orig( x_full(), newx, zero_order_orig_info() );
  VectorDenseMutableEncap  c_d(*zero_order_info.c);
  equ_from_full(
    m_orig_   ? c_orig_()                     : DVectorSlice()
    ,mI_orig_ ? h_orig_()                     : DVectorSlice()
    ,mI_orig_ ? x_full()(n_orig_+1,n_full_)   : DVectorSlice() // s_orig
    ,&c_d()
    );
}

void NLPSerialPreprocess::imp_calc_c_breve(
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    &zero_order_info_breve
  ) const
{
  assert_initialized();
  VectorDenseEncap x_d(x);
  set_x_full( x_d(), newx, &x_full_() );
  if( m_orig_ )
    imp_calc_c_orig( x_full(), newx, zero_order_orig_info() );
  if( mI_orig_ )
    imp_calc_h_orig( x_full(), newx, zero_order_orig_info() );
  VectorDenseMutableEncap  c_breve_d(*zero_order_info_breve.c);
  c_breve_d() = c_orig_();
}

void NLPSerialPreprocess::imp_calc_h_breve(
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    &zero_order_info_breve
  ) const
{
  // If this function gets called then this->mI() > 0 must be true
  // which means that convert_inequ_to_equ must be false!
  assert_initialized();
  VectorDenseEncap  x_d(x);
  set_x_full( x_d(), newx, &x_full_() );
  imp_calc_h_orig( x_full(), newx, zero_order_orig_info() );
  VectorDenseMutableEncap  h_breve_d(*zero_order_info_breve.h);
  h_breve_d() = h_orig_(); // Nothing fancy right now
}

// Overridden protected members from NLPObjGrad

void NLPSerialPreprocess::imp_calc_Gf(
  const Vector            &x
  ,bool                   newx
  ,const ObjGradInfo      &obj_grad_info
  ) const
{
  using DenseLinAlgPack::Vt_S;
  assert_initialized();
  VectorDenseEncap  x_d(x);
  set_x_full( x_d(), newx, &x_full_() );
  if( n_full_ > n_orig_ ) Gf_full_(n_orig_+1,n_full_) = 0.0; // Initialize terms for slacks to zero!
  imp_calc_Gf_orig( x_full(), newx, obj_grad_orig_info() );
  VectorDenseMutableEncap  Gf_d(*obj_grad_info.Gf);
  var_from_full( Gf_full_.begin(), Gf_d().begin() );
  Vt_S( &Gf_d(), scale_f_ );
}

// protected members

bool NLPSerialPreprocess::imp_get_next_basis(
  IVector      *var_perm_full
  ,IVector     *equ_perm_full
  ,size_type   *rank_full
  ,size_type   *rank
  )
{
  return false; // default is that the subclass does not select the basis
}

void NLPSerialPreprocess::assert_initialized() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !initialized_, UnInitialized
    ,"NLPSerialPreprocess : The nlp has not been initialized yet" );
}

void NLPSerialPreprocess::set_x_full(
  const DVectorSlice& x, bool newx
  ,DVectorSlice* x_full
  ) const
{
  DenseLinAlgPack::assert_vs_sizes(x.dim(),n_);
  if(newx)
    var_to_full(x.begin(), x_full->begin());
}

void NLPSerialPreprocess::var_from_full(
  DVectorSlice::const_iterator   vec_full
  ,DVectorSlice::iterator        vec
  ) const
{
//	std::cout << "\nvar_from_full(...) : ";
  for(size_type i = 1; i <= n_; i++) {
    *vec++ = vec_full[var_full_to_fixed_(var_perm_(i)) - 1];
//		std::cout
//			<< "\ni = " << i
//			<< "\nvar_perm_(i) = " << var_perm_(i)
//			<< "\nvar_full_to_fixed_(var_perm_(i)) = " << var_full_to_fixed_(var_perm_(i))
//			<< "\nvec_full[var_full_to_fixed_(var_perm_(i)) - 1] = " << vec_full[var_full_to_fixed_(var_perm_(i)) - 1]
//			<< "\nvec[i] = " << *(vec-1) << "\n\n";
  }		
}

void NLPSerialPreprocess::var_to_full( DVectorSlice::const_iterator vec, DVectorSlice::iterator vec_full ) const
{
  for(size_type i = 1; i <= n_; i++)
    vec_full[var_full_to_fixed_(var_perm_(i)) - 1] = *vec++;
}

void NLPSerialPreprocess::equ_from_full(
  const DVectorSlice   &c_orig
  ,const DVectorSlice  &h_orig
  ,const DVectorSlice  &s_orig
  ,DVectorSlice        *c_full
  ) const
{
  size_type i;
  // c_full = [ c_orig; h_orig - s_orig ]
   for(i = 1; i <= m_orig_; ++i)
    (*c_full)(inv_equ_perm_(i)) = c_orig(i);
   for(i = 1; i <= mI_orig_; ++i)
    (*c_full)(inv_equ_perm_(m_orig_+i)) = h_orig(i) - s_orig(i);
}

// private members

bool NLPSerialPreprocess::get_next_basis_remove_fixed(
  IVector* var_perm, IVector* equ_perm, size_type* rank
  )
{
  IVector var_perm_full(n_full_);
  equ_perm->resize(m_full_);
  size_type rank_full, rank_fixed_removed;
  if( imp_get_next_basis( &var_perm_full, equ_perm, &rank_full, &rank_fixed_removed ) ) {
//		std::cerr << "var_perm_full =\n"  << var_perm_full;
//		std::cerr << "equ_perm =\n"  << *equ_perm;
//		std::cerr << "rank_full = "  << rank_full << std::endl;
    //
    // The NLP subclass has another basis to select
    //
    // Translate the basis by removing variables fixed by bounds.
    // This is where it is important that var_perm_full is
    // sorted in assending order for basis and non-basis variables
    //
    // This is a little bit of a tricky algorithm.  We have to
    // essentially loop through the set of basic and non-basic
    // variables, remove fixed variables and adjust the indexes
    // of the remaining variables.  Since the set of indexes in
    // the basic and non-basis sets are sorted, this will not
    // be too bad of an algorithm.

    // Iterator for the fixed variables that we are to remove
    IVector::const_iterator     fixed_itr     = var_full_to_fixed_.begin() + n_;
    IVector::const_iterator     fixed_end     = var_full_to_fixed_.end();

    // Iterator for the basis variables
    IVector::iterator           basic_itr  = var_perm_full.begin();
    IVector::iterator           basic_end  = basic_itr + rank_full;

    // Iterator for the non-basis variables
    IVector::iterator           nonbasic_itr  = basic_end;
    IVector::iterator           nonbasic_end  = var_perm_full.end();
    
    // Count the number of fixed basic and non-basic variables
    index_type
      count_fixed          = 0,
      count_basic_fixed    = 0,
      count_nonbasic_fixed = 0;

    // Loop through all of the fixed variables and remove and compact
    for( ; fixed_itr != fixed_end; ++fixed_itr ) {
      const index_type
        next_fixed = ( fixed_itr +1 != fixed_end ? *(fixed_itr+1) : n_full_+1);
      // Bring the basic and nonbasic variables up to this fixed variable
      for( ; *basic_itr < *fixed_itr; ++basic_itr )
        *(basic_itr - count_basic_fixed) = *basic_itr - count_fixed;
      for( ; *nonbasic_itr < *fixed_itr; ++nonbasic_itr )
        *(nonbasic_itr - count_nonbasic_fixed) = *nonbasic_itr - count_fixed;
      // Update the count of the fixed variables
      if( *basic_itr == *fixed_itr ) {
        ++count_basic_fixed;
        ++basic_itr;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPT( !( *nonbasic_itr == *fixed_itr ) ); // If basic was not fixed then nonbasic better be!
        ++count_nonbasic_fixed;
        ++nonbasic_itr;

      }
      ++count_fixed;
      // Now update the indexes until the next fixed variable
      for( ; *basic_itr < next_fixed; ++basic_itr )
        *(basic_itr - count_basic_fixed) = *basic_itr - count_fixed;
      for( ; *nonbasic_itr < next_fixed; ++nonbasic_itr )
        *(nonbasic_itr - count_nonbasic_fixed) = *nonbasic_itr - count_fixed;
    }
    TEUCHOS_TEST_FOR_EXCEPT( !( count_fixed == n_full_ - n_ ) ); // Basic check

    var_perm->resize(n_);
    std::copy(
      var_perm_full.begin()
      ,var_perm_full.begin() + rank_fixed_removed
      ,var_perm->begin()
      );
    std::copy(
      var_perm_full.begin() + rank_full
      ,var_perm_full.begin() + rank_full + (n_-rank_fixed_removed)
      ,var_perm->begin() + rank_fixed_removed
      );
    *rank = rank_fixed_removed;
    return true;
  }
  else {
  
  // try to find the file giving the basis...
  char ind[17];
  sprintf(ind, "%d", basis_selection_num_);
  std::string fname = "basis_";
  fname += ind;
  fname += ".sel";
  basis_selection_num_++;

  std::ifstream basis_file(fname.c_str());
  if (basis_file)
    {
    // try to read the basis file
    std::string tags;

    int n;
    basis_file >> tags;
    TEUCHOS_TEST_FOR_EXCEPTION(
      tags != "n", std::logic_error
      ,"Incorrect basis file format - \"n\" expected, \"" << tags << "\" found.");
    basis_file >> n;
    TEUCHOS_TEST_FOR_EXCEPTION(
      n <= 0, std::logic_error
      , "Incorrect basis file format - n > 0 \"" << n << "\" found.");

    int m;
    basis_file >> tags;
    TEUCHOS_TEST_FOR_EXCEPTION(
      tags != "m", std::logic_error
      ,"Incorrect basis file format - \"m\" expected, \"" << tags << "\" found.");
    basis_file >> m;
    TEUCHOS_TEST_FOR_EXCEPTION(
      m > n , std::logic_error
      ,"Incorrect basis file format - 0 < m <= n expected, \"" << m << "\" found.");
    
    int r;
    basis_file >> tags;
    TEUCHOS_TEST_FOR_EXCEPTION(
      tags != "rank", std::logic_error
      ,"Incorrect basis file format - \"rank\" expected, \"" << tags << "\" found.");
    basis_file >> r;
    TEUCHOS_TEST_FOR_EXCEPTION(
      r > m, std::logic_error
      ,"Incorrect basis file format - 0 < rank <= m expected, \"" << r << "\" found.");		
    if (rank)
      { *rank = r; }

    // var_permutation
    basis_file >> tags;
    TEUCHOS_TEST_FOR_EXCEPTION(
      tags != "var_perm", std::logic_error
      ,"Incorrect basis file format -\"var_perm\" expected, \"" << tags << "\" found.");
    var_perm->resize(n);
    {for (int i=0; i < n; i++)
      {
      int var_index;
      basis_file >> var_index;
      TEUCHOS_TEST_FOR_EXCEPTION(
        var_index < 1 || var_index > n, std::logic_error
        ,"Incorrect basis file format for var_perm: 1 <= indice <= n expected, \"" << n << "\" found.");
      (*var_perm)[i] = var_index;
      }}

    // eqn_permutation
    basis_file >> tags;
    TEUCHOS_TEST_FOR_EXCEPTION(
      tags != "equ_perm", std::logic_error
      ,"Incorrect basis file format -\"equ_perm\" expected, \"" << tags << "\" found.");
    equ_perm->resize(r);
    {for (int i=0; i < r; i++)
      {
      int equ_index;
      basis_file >> equ_index;
      TEUCHOS_TEST_FOR_EXCEPTION(
        equ_index < 1 || equ_index > m, std::logic_error
        ,"Incorrect basis file format for equ_perm: 1 <= indice <= m expected, \"" << m << "\" found.");
      (*equ_perm)[i] = equ_index;
      }}

    return true;
    }
  }

  return false;
}

void NLPSerialPreprocess::assert_and_set_basis(
  const IVector& var_perm, const IVector& equ_perm, size_type rank
  )
{
  namespace mmp = MemMngPack;

  // Assert that this is a valid basis and set the internal basis.  Also repivot 'xinit', 
  // 'xl', and 'xu'.

  const value_type inf_bnd = NLPSerialPreprocess::infinite_bound();

  // Assert the preconditions
  TEUCHOS_TEST_FOR_EXCEPTION(
    var_perm.size() != n_ || equ_perm.size() != m_full_, std::length_error
    ,"NLPSerialPreprocess::set_basis():  The sizes "
    "of the permutation vectors does not match the size of the NLP" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    rank > m_full_, InvalidBasis
    ,"NLPSerialPreprocess::set_basis():  The rank "
    "of the basis can not be greater that the number of constraints" );
    
  // Set the basis
  r_ = rank;
  if( &var_perm_ != &var_perm )
    var_perm_ = var_perm;
  if( &equ_perm_ != &equ_perm )
    equ_perm_ = equ_perm;
  DenseLinAlgPack::inv_perm( equ_perm_, &inv_equ_perm_ );

  var_from_full( xinit_full_().begin(), xinit_.set_vec().begin() );
  if(num_bounded_x_) {
    var_from_full( xl_full_().begin(), xl_.set_vec().begin() );
    var_from_full( xu_full_().begin(), xu_.set_vec().begin() );
    do_force_xinit_in_bounds();
  }
  else {
    xl_ = -NLP::infinite_bound();
    xu_ = +NLP::infinite_bound();
  }
  P_var_.initialize(Teuchos::rcp(new IVector(var_perm)),Teuchos::null);
  P_equ_.initialize(Teuchos::rcp(new IVector(equ_perm)),Teuchos::null);
}

void NLPSerialPreprocess::assert_bounds_on_variables() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(imp_has_var_bounds() || n_full_ > n_orig_), NLP::NoBounds
    ,"There are no bounds on the variables for this NLP" );
}

void NLPSerialPreprocess::do_force_xinit_in_bounds()
{
  AbstractLinAlgPack::force_in_bounds( xl_, xu_, &xinit_ );
}

} // end namespace NLPInterfacePack
