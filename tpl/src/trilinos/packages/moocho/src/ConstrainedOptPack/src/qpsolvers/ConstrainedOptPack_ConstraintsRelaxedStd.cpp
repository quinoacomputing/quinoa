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
//
// ToDo: 12/29/00: Consider scaling when determining if a
// constraint is violated.  We should consider the size of
// ||d(e[j](x))/d(x)||inf but this is expensive to compute
// given the current interfaces.  We really need to rectify
// this!
//

#include <assert.h>

#include <limits>

#include "ConstrainedOptPack_ConstraintsRelaxedStd.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "Teuchos_Assert.hpp"

namespace {

ConstrainedOptPack::EBounds
convert_bnd_type( int bnd_type )
{
  switch(bnd_type) {
    case -1:
      return ConstrainedOptPack::LOWER;
    case 0:
      return ConstrainedOptPack::EQUALITY;
    case +1:
      return ConstrainedOptPack::UPPER;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return ConstrainedOptPack::LOWER; // Never be called
}

// Get an element from a sparse vector and return zero if it does not exist
AbstractLinAlgPack::value_type get_sparse_element(
  const AbstractLinAlgPack::SpVectorSlice& v
  ,AbstractLinAlgPack::size_type i
  )
{
  const AbstractLinAlgPack::SpVectorSlice::element_type
    *ele_ptr = v.lookup_element(i);
  return ele_ptr ? ele_ptr->value() : 0.0;
}

}	// end namespace

namespace ConstrainedOptPack {
namespace QPSchurPack {

// members for ConstraintsRelaxedStd

ConstraintsRelaxedStd::ConstraintsRelaxedStd()
  :inequality_pick_policy_(ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY)
  ,etaL_(0.0)
  ,dL_(NULL)
  ,dU_(NULL)
  ,eL_(NULL)
  ,eU_(NULL)
  ,Ed_(NULL)
  ,last_added_j_(0)
  ,last_added_bound_(0.0)
  ,last_added_bound_type_(FREE)
  ,next_undecomp_f_k_(0)
{}

void ConstraintsRelaxedStd::initialize(
  const VectorSpace::space_ptr_t   &space_d_eta
  ,value_type                      etaL
  ,const Vector              *dL
  ,const Vector              *dU
  ,const MatrixOp              *E
  ,BLAS_Cpp::Transp                trans_E
  ,const Vector              *b
  ,const Vector              *eL
  ,const Vector              *eU
  ,const MatrixOp              *F
  ,BLAS_Cpp::Transp                trans_F
  ,const Vector              *f
  ,size_type                       m_undecomp
  ,const size_type                 j_f_undecomp[]
  ,VectorMutable             *Ed
  ,bool                            check_F
  ,value_type                      bounds_tol
  ,value_type                      inequality_tol
  ,value_type                      equality_tol
  )
{
  size_type
    nd   = space_d_eta->dim() - 1,
    m_in = 0,
    m_eq = 0;

  TEUCHOS_TEST_FOR_EXCEPT( !(  m_undecomp == (F ? f->dim() : 0)  ) ); // ToDo: support decomposed equalities in future.

  // Validate that the correct sets of constraints are selected
  TEUCHOS_TEST_FOR_EXCEPTION(
    dL && !dU, std::invalid_argument
    ,"ConstraintsRelaxedStd::initialize(...) : Error, "
    "If dL!=NULL then dU!=NULL must also be true." );
  TEUCHOS_TEST_FOR_EXCEPTION(
    E && ( !b || !eL || !eU ), std::invalid_argument
    ,"ConstraintsRelaxedStd::initialize(...) : Error, "
    "If E!=NULL then b!=NULL, eL!=NULL and eU!=NULL must also be true." );
  TEUCHOS_TEST_FOR_EXCEPTION(
    F && !f, std::invalid_argument
    ,"ConstraintsRelaxedStd::initialize(...) : Error, "
    "If F!=NULL then f!=NULL must also be true." );

  // Validate input argument sizes
  if(dL) {
    const size_type dL_dim = dL->dim(), dU_dim = dU->dim();
    TEUCHOS_TEST_FOR_EXCEPTION(
      dL_dim != nd, std::invalid_argument
      ,"ConstraintsRelaxedStd::initialize(...) : Error, "
      "dL.dim() != d->dim()." );
    TEUCHOS_TEST_FOR_EXCEPTION(
      dU_dim != nd, std::invalid_argument 
      ,"ConstraintsRelaxedStd::initialize(...) : Error, "
      "dU.dim() != d->dim()." );
  }
  if(E) {
    const size_type
      E_rows = E->rows(),
      E_cols = E->cols(),
      opE_cols = BLAS_Cpp::cols( E_rows, E_cols, trans_E ),
      b_dim = b->dim(),
      eL_dim = eL->dim(),
      eU_dim = eU->dim(),
      Ed_dim = Ed ? Ed->dim() : 0;
    m_in = BLAS_Cpp::rows( E_rows, E_cols, trans_E );
    TEUCHOS_TEST_FOR_EXCEPTION(
      opE_cols != nd, std::invalid_argument
      ,"ConstraintsRelaxedStd::initialize(...) : Error, "
      "op(E).cols() != nd." );
    TEUCHOS_TEST_FOR_EXCEPTION(
      b_dim != m_in, std::invalid_argument
      ,"ConstraintsRelaxedStd::initialize(...) : Error, "
      "b->dim() != op(E).rows()." );
    TEUCHOS_TEST_FOR_EXCEPTION(
      eL_dim != m_in, std::invalid_argument
      ,"ConstraintsRelaxedStd::initialize(...) : Error, "
      "eL->dim() != op(E).rows()." );
    TEUCHOS_TEST_FOR_EXCEPTION(
      eU_dim != m_in, std::invalid_argument
      ,"ConstraintsRelaxedStd::initialize(...) : Error, "
      "eU->dim() != op(E).rows()." );
    TEUCHOS_TEST_FOR_EXCEPTION(
      Ed && Ed_dim != m_in, std::invalid_argument
      ,"ConstraintsRelaxedStd::initialize(...) : Error, "
      "Ed->dim() != op(E).rows()." );
  }
  if(F) {
    const size_type
      F_rows = F->rows(),
      F_cols = F->cols(),
      opF_cols = BLAS_Cpp::cols( F_rows, F_cols, trans_F ),
      f_dim = f->dim();
    m_eq = BLAS_Cpp::rows( F_rows, F_cols, trans_F );
    TEUCHOS_TEST_FOR_EXCEPTION(
      opF_cols != nd, std::invalid_argument
      ,"QPSolverRelaxed::solve_qp(...) : Error, "
      "op(F).cols() != nd." );
    TEUCHOS_TEST_FOR_EXCEPTION(
      f_dim != m_eq, std::invalid_argument
      ,"QPSolverRelaxed::solve_qp(...) : Error, "
      "f->dim() != op(F).rows()." );
  }
  
  // Initialize other members
  A_bar_.initialize(
    space_d_eta,m_in,m_eq,E,trans_E,b,F,trans_F,f,m_undecomp,j_f_undecomp);
  etaL_				= etaL;
  dL_					= dL;
  dU_					= dU;
  eL_					= eL;
  eU_					= eU;
  Ed_					= Ed;
  check_F_			= check_F;
  bounds_tol_			= bounds_tol;
  inequality_tol_		= inequality_tol;
  equality_tol_		= equality_tol;
  last_added_j_		= 0;	// No cached value.
  next_undecomp_f_k_	= m_undecomp ? 1 : 0; // Check the first undecomposed equality

}

const ConstraintsRelaxedStd::MatrixConstraints&
ConstraintsRelaxedStd::A_bar_relaxed() const
{
  return A_bar_;
}

// Overridden from Constraints

size_type ConstraintsRelaxedStd::n() const
{
  return A_bar_.nd() + 1;
}

size_type ConstraintsRelaxedStd::m_breve() const
{
  return A_bar_.m_in() + A_bar_.m_eq();
}

const MatrixOp& ConstraintsRelaxedStd::A_bar() const
{
  return A_bar_;
}

void ConstraintsRelaxedStd::pick_violated_policy( EPickPolicy pick_policy )
{
  switch(pick_policy) {
    case ANY_VIOLATED:
      inequality_pick_policy_ = ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY;
      break;
    case MOST_VIOLATED:
      inequality_pick_policy_ = ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
}

Constraints::EPickPolicy
ConstraintsRelaxedStd::pick_violated_policy() const
{
  switch(inequality_pick_policy_) {
    case ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY:
      return ANY_VIOLATED;
    case ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY:
      return ANY_VIOLATED;
    case ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY:
      return MOST_VIOLATED;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return ANY_VIOLATED;	// will never be executed
}

void ConstraintsRelaxedStd::pick_violated(
  const DVectorSlice& x_in, size_type* j_viol, value_type* constr_val
  ,value_type* viol_bnd_val, value_type* norm_2_constr, EBounds* bnd, bool* can_ignore
  ) const
{
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  using AbstractLinAlgPack::max_inequ_viol;
  using LinAlgOpPack::V_VmV;

  TEUCHOS_TEST_FOR_EXCEPTION(
    x_in.dim() != A_bar_.nd()+1, std::length_error
    ,"ConstraintsRelaxedStd::pick_violated(...) : Error, "
    "x is the wrong length" );

  const size_type
    nd = A_bar_.nd();

  // Get a version of x = [ d; eta ] in the correct vector object
  VectorSpace::vec_mut_ptr_t
    x = A_bar_.space_cols().create_member();
  (VectorDenseMutableEncap(*x)()) = x_in;
  VectorSpace::vec_mut_ptr_t
    d = x->sub_view(1,nd);
  const value_type
    eta = x->get_ele(nd+1);

  bool Ed_computed = false;

  // //////////////////////////////////////////////
  // Check the equality constraints first
  if( check_F_ && A_bar_.F() ) {
    TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Update below code!
/*
    // The basic strategy here is to go through all of the equality
    // constraints first and add all of the ones that are violated by
    // more than the set tolerance.  Those that are not sufficiently
    // violated are passed over but they are remembered for later.
    // Only when all of the constraints have been gone through once
    // will those passed over constraints be considered and then only
    // if ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY is selected.
    const GenPermMatrixSlice& P_u = A_bar_.P_u();
    size_type e_k_mat_row, e_k_mat_col = 1; // Mapping matrix for e(j)
    GenPermMatrixSlice e_k_mat; // ToDo: Change to EtaVector
    if( next_undecomp_f_k_ <= P_u.nz() ) {
      // This is our first pass through the undecomposed equalities.
      GenPermMatrixSlice::const_iterator P_u_itr, P_u_end; // only used if P_u is not identity
      if( !P_u.is_identity() ) {
        P_u_itr = P_u.begin() + (next_undecomp_f_k_ - 1);
        P_u_end = P_u.end();
      }
      for( ; next_undecomp_f_k_ <=  P_u.nz(); ) {
        size_type k = 0;
        if( !P_u.is_identity() ) {
          k = P_u_itr->row_i();
          ++P_u_itr;
        }
        else {
          k = next_undecomp_f_k_;
        }
        ++next_undecomp_f_k_;
        // evaluate the constraint e(k)'*[op(F)*d + (1-eta)*f]
        value_type
          Fd_k = 0.0,
          f_k = (*A_bar_.f())(k);
        e_k_mat.initialize(
          A_bar_.m_eq(),1,1,0,0,GPMSIP::BY_ROW_AND_COL
          ,&(e_k_mat_row=k),&e_k_mat_col,false );
        DVectorSlice Fd_k_vec(&Fd_k,1);
        AbstractLinAlgPack::Vp_StPtMtV(   // ToDo: Use transVtMtV(...) instead!
          &Fd_k_vec, 1.0, e_k_mat, BLAS_Cpp::trans
          ,*A_bar_.F(), A_bar_.trans_F(), d, 0.0 );
        const value_type
          err = ::fabs(Fd_k + (1.0 - eta)*f_k) / (1.0 + ::fabs(f_k));
        if( err > equality_tol() ) {
          *j_viol         = nd + 1 + A_bar_.m_in() + k;
          *constr_val     = Fd_k - eta*f_k;
          *norm_2_constr  = 1.0; // ToDo: Compute it some how?
          *bnd            = EQUALITY;
          *viol_bnd_val   = -f_k;
          *can_ignore     = false; // Given this careful algorithm this should be false
          // cache this
          last_added_j_			= *j_viol;
          last_added_bound_type_	= *bnd;
          last_added_bound_		= *viol_bnd_val;
          return;
        }
        else {
          passed_over_equalities_.push_back(k); // remember for later
        }
      }
    }
    else if(
      inequality_pick_policy() == ADD_MOST_VIOLATED_BOUNDS_AND_INEQUALITY
      && passed_over_equalities_.size() > 0
      )
    {
      // Now look through the list of passed over equalities and see which one
      // is violated.  If a passed over constraint is added to the active set
      // then it is removed from this list.
      TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement!
    }
*/
  }

  // /////////////////////////////////////////////
  // Find the most violated variable bound.

  size_type       max_bound_viol_j       = 0;
  value_type      max_bound_viol         = 0.0;
  value_type      max_bound_d_viol       = 0.0;
  value_type      max_bound_dLU_viol     = 0.0;
  int             max_bound_viol_type    = -2;
  if( dL_ && ( dL_->nz() || dU_->nz() ) ) {
    // dL <= d <= dU
    max_inequ_viol(
      *d, *dL_, *dU_
      ,&max_bound_viol_j, &max_bound_viol
      ,&max_bound_d_viol, &max_bound_viol_type, &max_bound_dLU_viol 
      );
    if(  max_bound_viol > bounds_tol_ ) {
      // Set the return values
      *j_viol         = max_bound_viol_j;
      *constr_val     = max_bound_d_viol;
      *norm_2_constr  = 1.0; // This is correct
      *bnd            = convert_bnd_type(max_bound_viol_type);
      *viol_bnd_val   = max_bound_dLU_viol;
      *can_ignore     = false;
    }
    else {
      max_bound_viol_j = 0;	// No variable bounds sufficiently violated.
    }
  }

  if( (	inequality_pick_policy_ == ADD_BOUNDS_THEN_MOST_VIOLATED_INEQUALITY
      ||	inequality_pick_policy_ == ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY )
    && max_bound_viol_j
    )
  {
    // A variable bound has been violated so lets just return it!
    last_added_j_           = *j_viol;
    last_added_bound_type_  = *bnd;
    last_added_bound_       = *viol_bnd_val;
    return;
  }

  // /////////////////////////////////////////////
  // Check the general inequalities

  size_type       max_inequality_viol_j        = 0;
  value_type      max_inequality_viol          = 0.0;
  value_type      max_inequality_e_viol        = 0.0;
  value_type      max_inequality_eLU_viol      = 0.0;
  int             max_inequality_viol_type     = -2;

  if( inequality_pick_policy_ == ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY ) {
    // Find the first general inequality violated by more than
    // the defined tolerance.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error
      ,"ConstraintsRelaxedStd::pick_violated(...) : Error,\n"
      "The option ADD_BOUNDS_THEN_FIRST_VIOLATED_INEQUALITY has not been implemented yet\n" );
  }
  else {
    // Find the most violated inequality constraint
    if( A_bar_.m_in() &&  ( eL_->nz() || eU_->nz() ) ) {
      // e = op(E)*d - b*eta
      VectorSpace::vec_mut_ptr_t e = eL_->space().create_member();
      LinAlgOpPack::V_MtV( e.get(), *A_bar_.E(), A_bar_.trans_E(), *d );
      if(Ed_) {
        *Ed_        = *e;
        Ed_computed = true;
      }
      LinAlgOpPack::Vp_StV( e.get(), -eta, *A_bar_.b() );
      // eL <= e <= eU
      max_inequ_viol(
        *e, *eL_, *eU_
        ,&max_inequality_viol_j, &max_inequality_viol
        ,&max_inequality_e_viol, &max_inequality_viol_type, &max_inequality_eLU_viol 
        );
      if( max_inequality_viol > max_bound_viol
        && max_inequality_viol > inequality_tol_ )
      {
        *j_viol         = max_inequality_viol_j + nd + 1; // offset into A_bar
        *constr_val     = max_inequality_e_viol;
        *norm_2_constr  = 1.0;  // This is not correct!
        *bnd            = convert_bnd_type(max_inequality_viol_type);
        *viol_bnd_val   = max_inequality_eLU_viol;
        *can_ignore     = false;
      }
      else {
        max_inequality_viol_j = 0; // No general inequality constraints sufficiently violated.
      }
    }
  }

  if( max_bound_viol_j || max_inequality_viol_j ) {
    // One of the constraints has been violated so just return it.
    last_added_j_           = *j_viol;
    last_added_bound_type_  = *bnd;
    last_added_bound_       = *viol_bnd_val;
    return;
  }

  // If we get here then no constraint was found that violated any of the tolerances.
  if( Ed_ && !Ed_computed ) {
    // Ed = op(E)*d
    LinAlgOpPack::V_MtV( Ed_, *A_bar_.E(), A_bar_.trans_E(), *d );
  }
  *j_viol			= 0;
  *constr_val		= 0.0;
  *viol_bnd_val	= 0.0;
  *norm_2_constr	= 0.0;
  *bnd			= FREE;	 // Meaningless
  *can_ignore		= false; // Meaningless
}

void ConstraintsRelaxedStd::ignore( size_type j )
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"ConstraintsRelaxedStd::ignore(...) : Error, "
    "This operation is not supported yet!" );
}

value_type ConstraintsRelaxedStd::get_bnd( size_type j, EBounds bnd ) const
{
  const value_type inf = std::numeric_limits<value_type>::max();

  TEUCHOS_TEST_FOR_EXCEPTION(
    j > A_bar_.cols(), std::range_error
    ,"ConstraintsRelaxedStd::get_bnd(j,bnd) : Error, "
    "j = " << j << " is not in range [1," << A_bar_.cols() << "]" );

  // See if this is the last constraint we added to the active set.
  if( j == last_added_j_ && bnd == last_added_bound_type_ ) {
    return last_added_bound_;
  }

  // Lookup the bound! (sparse lookup)
  size_type j_local = j;
  const SpVectorSlice::element_type *ele_ptr = NULL;
  if( j_local <= A_bar_.nd() ) {
    if(dL_) {
      switch( bnd ) {
        case EQUALITY:
        case LOWER:
          return dL_->get_ele(j_local);
        case UPPER:
          return dU_->get_ele(j_local);
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true);
      }
    }
    else {
      return ( bnd == LOWER ? -1.0 : +1.0 ) * inf;
    }
  }
  else if( (j_local -= A_bar_.nd()) <= 1 ) {
    switch( bnd ) {
      case EQUALITY:
      case LOWER:
        return etaL_;
      case UPPER:
        return +inf;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
  else if( (j_local -= 1) <= A_bar_.m_in() ) {
    switch( bnd ) {
      case EQUALITY:
      case LOWER:
        return eL_->get_ele(j_local);
      case UPPER:
        return eU_->get_ele(j_local);
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
  else if( (j_local -= A_bar_.m_in()) <= A_bar_.m_eq() ) {
    switch( bnd ) {
      case EQUALITY:
      case LOWER:
      case UPPER:
        return -A_bar_.f()->get_ele(j_local);
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
  return 0.0;	// will never be executed!
}

void ConstraintsRelaxedStd::cache_last_added(
  size_type last_added_j, value_type last_added_bound
  ,EBounds last_added_bound_type
  ) const
{
  last_added_j_           = last_added_j;
  last_added_bound_       = last_added_bound;
  last_added_bound_type_  = last_added_bound_type;
}

// members for ConstraintsRelaxedStd::MatrixConstraints

ConstraintsRelaxedStd::MatrixConstraints::MatrixConstraints()
  :nd_(0)
  ,m_in_(0)
  ,m_eq_(0)
  ,E_(NULL)
  ,trans_E_(BLAS_Cpp::no_trans)
  ,b_(NULL)
  ,F_(NULL)
  ,trans_F_(BLAS_Cpp::no_trans)
  ,f_(NULL)
  ,space_cols_(Teuchos::null)
  ,space_rows_(NULL,0)
{}

void ConstraintsRelaxedStd::MatrixConstraints::initialize(
  const VectorSpace::space_ptr_t   &space_d_eta
  ,const size_type                 m_in
  ,const size_type                 m_eq
  ,const MatrixOp                  *E
  ,BLAS_Cpp::Transp                trans_E
  ,const Vector                    *b
  ,const MatrixOp                  *F
  ,BLAS_Cpp::Transp                trans_F
  ,const Vector                    *f
  ,size_type                       m_undecomp
  ,const size_type                 j_f_undecomp[]
  )
{
  namespace mmp = MemMngPack;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;

  const size_type nd = space_d_eta->dim() - 1;

  // Setup P_u
  const bool test_setup = true; // Todo: Make this an argument!
  if( m_undecomp > 0 && f->dim() > m_undecomp ) {
    P_u_row_i_.resize(m_undecomp);
    P_u_col_j_.resize(m_undecomp);
    const size_type
      *j_f_u = j_f_undecomp;  // This is sorted by row!
    row_i_t::iterator
      row_i_itr = P_u_row_i_.begin();
    col_j_t::iterator
      col_j_itr = P_u_col_j_.begin();
    for( size_type i = 1; i <= m_undecomp; ++i, ++j_f_u, ++row_i_itr, ++col_j_itr ) {
      *row_i_itr = *j_f_u; // This is sorted in asscending order!
      *col_j_itr = i;
    }
    P_u_.initialize(nd,m_undecomp,m_undecomp,0,0,GPMSIP::BY_ROW_AND_COL
      ,&P_u_row_i_[0],&P_u_col_j_[0],test_setup);
  }
  else if( m_undecomp > 0) { // Must be == f->dim()
    // Set to identity
    P_u_.initialize(m_undecomp,m_undecomp,GenPermMatrixSlice::IDENTITY_MATRIX);
  }

  // space_cols_
  space_cols_ = space_d_eta;

  // space_rows_
  VectorSpace::space_ptr_t  row_spaces[3];
  int num_row_spaces = 1;
  row_spaces[0] = space_d_eta;
  if(m_in)
    row_spaces[num_row_spaces++] = Teuchos::rcp(
      trans_E == BLAS_Cpp::no_trans ? &E->space_cols() : &E->space_rows()
      ,false
      );
  if(m_eq) {
    VectorSpace::space_ptr_t
      vs = Teuchos::rcp(
        trans_F == BLAS_Cpp::no_trans ? &F->space_cols() : &F->space_rows()
        ,false
        );
    if(m_undecomp)
      vs = vs->space(P_u_,BLAS_Cpp::trans);
    row_spaces[num_row_spaces++] = vs;
  }
  space_rows_.initialize( row_spaces, num_row_spaces, space_d_eta->small_vec_spc_fcty() );

  // Set the rest of the members
  nd_       = space_d_eta->dim() - 1;
  m_in_     = m_in;
  m_eq_     = m_eq;
  E_        = E;
  trans_E_  = trans_E;
  b_        = b;
  F_        = F;
  trans_F_  = trans_F;
  f_        = f;

}

// Overridden from MatrixOp

const VectorSpace& ConstraintsRelaxedStd::MatrixConstraints::space_cols() const
{
  return *space_cols_;
}

const VectorSpace& ConstraintsRelaxedStd::MatrixConstraints::space_rows() const
{
  return space_rows_;
}

MatrixOp& ConstraintsRelaxedStd::MatrixConstraints::operator=(
  const MatrixOp& M
  )
{
  // ToDo: Finish me
  TEUCHOS_TEST_FOR_EXCEPT(true);
  return *this;
}

/* 10/25/00 I don't think we need this function yet!
void ConstraintsRelaxedStd::MatrixConstraints::Mp_StPtMtP(
  DMatrixSlice* C, value_type a
  ,const GenPermMatrixSlice& P1, BLAS_Cpp::Transp P1_trans
  ,BLAS_Cpp::Transp M_trans
  ,const GenPermMatrixSlice& P2, BLAS_Cpp::Transp P2_trans
  )const
{
  using BLAS_Cpp::trans_not;
  using AbstractLinAlgPack::dot;
  using AbstractLinAlgPack::Vp_StMtV;
  using AbstractLinAlgPack::Vp_StPtMtV;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;

  //	
  //	A_bar = [  I   0  op(E')   op(F')  ]
  //	        [  0   1   -b'      -f'    ]
  //

  const size_type
    E_start = nd() + 1 + 1,
    F_start = E_start + m_in(),
    F_end	= F_start + m_eq() - 1;
  const Range1D
    d_rng = Range1D(1,nd()),
    E_rng = m_in() ? Range1D(E_start,F_start-1) : Range1D(),
    F_rng = m_eq() ? Range1D(F_start,F_end) : Range1D();

  // For this to work (as shown below) we need to have P1 sorted by
  // row if op(P1) = P1' or sorted by column if op(P1) = P1.
  // Also, we must have P1 sorted by
  // row if op(P2) = P2 or sorted by column if op(P2) = P2'
  // If P1 or P2 are not sorted properly, we will just use the default
  // implementation of this operation.
  if( 	( P1.ordered_by() == GPMSIP::BY_ROW && P1_trans == BLAS_Cpp::no_trans )
      || 	( P1.ordered_by() == GPMSIP::BY_COL && P1_trans == BLAS_Cpp::trans )
      || 	( P2.ordered_by() == GPMSIP::BY_ROW && P2_trans == BLAS_Cpp::trans )
      || 	( P2.ordered_by() == GPMSIP::BY_COL && P2_trans == BLAS_Cpp::no_trans ) )
  {
    // Call the default implementation
    MatrixOp::Vp_StPtMtV(C,a,P1,P1_trans,M_trans,P2,P2_trans);
    return;
  }

  if( M_trans == BLAS_Cpp::no_trans ) {
    //
    // C += a * op(P1) * A_bar * op(P2)
    //
    //    = a * [ op(P11)  op(P12) ] * [ I  0  op(E')  op(F') ] * [ op(P21) ]
    //                                 [ 0  1    -b'    -f'   ]   [ op(P22) ]
    //                                                            [ op(P23) ]
    //                                                            [ op(P24) ]
    //
    // C +=   a*op(P11)*op(P21) + a*op(P21)*op(P22)
    //      + a*op(P11)*op(E')*op(P23) - a*op(P12)*b'*op(P23)
    //      + a*op(P11)*op(F')*op(P24) - a*op(P12)*f'*op(P24)
    //      

    TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: Implement This!

  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true);	// ToDo: Finish This!
  }
}
*/

void ConstraintsRelaxedStd::MatrixConstraints::Vp_StMtV(
  VectorMutable* y, value_type a, BLAS_Cpp::Transp trans_rhs1
  ,const Vector& x, value_type b
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  !F_ || P_u_.cols() == f_->dim()  ) ); // ToDo: Add P_u when needed!

  namespace mmp = MemMngPack;
  using BLAS_Cpp::trans_not;
  using AbstractLinAlgPack::dot;
  using LinAlgOpPack::Vt_S;
  using LinAlgOpPack::Vp_StV;
  using LinAlgOpPack::Vp_StMtV;

  // ToDo: Replace with proper check!
//	LinAlgOpPack::Vp_MtV_assert_sizes(y->dim(),rows(),cols(),trans_rhs1,x.dim());

  //	
  //	A_bar = [  I   0  op(E')   op(F')  ]
  //	        [  0   1   -b'      -f'    ]
  //

  const size_type
    E_start = nd() + 1 + 1,
    F_start = E_start + m_in(),
    F_end	= F_start + m_eq() - 1;
  const Range1D
    d_rng = Range1D(1,nd()),
    E_rng = m_in() ? Range1D(E_start,F_start-1) : Range1D(),
    F_rng = m_eq() ? Range1D(F_start,F_end) : Range1D();

  // y = b * y
  Vt_S( y, b );
  
  if( trans_rhs1 == BLAS_Cpp::no_trans ) {
    //
    // y += a* A_bar * x
    // 
    //   += a * [ I  0  op(E')  op(F') ] * [ x1 ]
    //          [ 0  1   -b'     -f'   ]   [ x2 ]
    //                                     [ x3 ]
    //                                     [ x4 ]
    //
    // [ y1 ]  += [ a * x1 + a * op(E') * x3 + a * op(F') * x4 ]
    // [ y2 ]     [ a * x2 - a * b' * x3     - a * f' * x4     ]
    //
    VectorMutable::vec_mut_ptr_t
      y1 = y->sub_view(d_rng);
    value_type
      y2 = y->get_ele(nd()+1);
    Vector::vec_ptr_t
      x1 = x.sub_view(d_rng);
    const value_type
      x2 = x.get_ele(nd()+1);
    Vector::vec_ptr_t
      x3 = m_in() ? x.sub_view(E_rng) : Teuchos::null,
      x4 = m_eq() ? x.sub_view(F_rng) : Teuchos::null;
    
    // [ y1 ]  += [ a * x1 + a * op(E') * x3 + a * op(F') * x4 ]
    Vp_StV( y1.get(), a, *x1 );
    if( m_in() )
      Vp_StMtV( y1.get(), a, *E(), trans_not( trans_E() ), *x3 );
    if( m_eq() )
      Vp_StMtV( y1.get(), a, *F(), trans_not( trans_F() ), *x4 );
    // [ y2 ]  += [ a * x2 - a * b' * x3     - a * f' * x4     ]
    y2 += a * x2;
    if( m_in() )
      y2 += - a * dot( *this->b(), *x3 );
    if( m_eq() )
      y2 += - a * dot( *f(), *x4 );
    y->set_ele(nd()+1,y2);
  }
  else if ( trans_rhs1 == BLAS_Cpp::trans ) {
    //
    // y += a* A_bar' * x
    // 
    //   += a * [ I       0 ] * [ x1 ]
    //          [ 0       1 ]   [ x2 ]
    //          [ op(E)  -b ]
    //          [ op(F)  -f ]
    //
    // [ y1 ]    [ a * x1                        ]
    // [ y2 ]    [                + a * x2       ]
    // [ y3 ] += [ a * op(E) * x1 - a * b * x2   ]
    // [ y4 ]    [ a * op(F) * x1 - a * f * x2   ]
    //
    VectorMutable::vec_mut_ptr_t
      y1 = y->sub_view(d_rng);
    value_type
      y2 = y->get_ele(nd()+1);
    VectorMutable::vec_mut_ptr_t
      y3 = m_in() ? y->sub_view(E_rng) : Teuchos::null,
      y4 = m_eq() ? y->sub_view(F_rng) : Teuchos::null;
    Vector::vec_ptr_t
      x1 = x.sub_view(d_rng);
    const value_type
      x2 = x.get_ele(nd()+1);
    // y1 += a * x1
    Vp_StV( y1.get(), a, *x1 );
    // y2 += a * x2
    y2 += a * x2;
    y->set_ele(nd()+1,y2);
    // y3 += a * op(E) * x1 - (a*x2) * b
    if( m_in() ) {
      Vp_StMtV( y3.get(), a, *E(), trans_E(), *x1 );
      Vp_StV( y3.get(), - a * x2, *this->b() );
    }
    // y4 += a * op(F) * x1 - (a*x2) * f
    if( m_eq() ) {
      Vp_StMtV( y4.get(), a, *F(), trans_F(), *x1 );
      Vp_StV( y4.get(), - a * x2, *f() );
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true);	// Invalid trans value
  }
}

void ConstraintsRelaxedStd::MatrixConstraints::Vp_StPtMtV(
  VectorMutable* y_out, value_type a
  ,const GenPermMatrixSlice& P, BLAS_Cpp::Transp P_trans
  ,BLAS_Cpp::Transp M_trans
  ,const SpVectorSlice& x, value_type beta
  ) const
{
  MatrixOp::Vp_StPtMtV(y_out,a,P,P_trans,M_trans,x,beta); // ToDo: Update below code!
  
/*

  TEUCHOS_TEST_FOR_EXCEPT( !(  !F_ || P_u_.cols() == f_->dim()  ) ); // ToDo: Add P_u when needed!

  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using BLAS_Cpp::trans_not;
   using DenseLinAlgPack::dot;
  using AbstractLinAlgPack::Vp_StMtV;
  using AbstractLinAlgPack::Vp_StPtMtV;
  namespace GPMSIP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;

  AbstractLinAlgPack::VectorDenseMutableEncap y_d(*y_out);
  DVectorSlice *y = &y_d();

  DenseLinAlgPack::Vp_MtV_assert_sizes(
    y->dim(),P.rows(),P.cols(),P_trans
    ,BLAS_Cpp::rows( rows(), cols(), M_trans) );
  DenseLinAlgPack::Vp_MtV_assert_sizes(
    BLAS_Cpp::cols( P.rows(), P.cols(), P_trans)
    ,rows(),cols(),M_trans,x.dim());

  //	
  //	A_bar = [  I   0  op(E')   op(F')  ]
  //	        [  0   1   -b'      -f'    ]
  //

  const size_type
    E_start = nd() + 1 + 1,
    F_start = E_start + m_in(),
    F_end	= F_start + m_eq() - 1;
  const Range1D
    d_rng = Range1D(1,nd()),
    E_rng = m_in() ? Range1D(E_start,F_start-1) : Range1D(),
    F_rng = m_eq() ? Range1D(F_start,F_end) : Range1D();

  // For this to work (as shown below) we need to have P sorted by
  // row if op(P) = P' or sorted by column if op(P) = P.  If
  // P is not sorted properly, we will just use the default
  // implementation of this operation.
  if( 	( P.ordered_by() == GPMSIP::BY_ROW && P_trans == BLAS_Cpp::no_trans )
      || 	( P.ordered_by() == GPMSIP::BY_COL && P_trans == BLAS_Cpp::trans )
    ||  ( P.ordered_by() == GPMSIP::UNORDERED ) )
  {
    // Call the default implementation
    //MatrixOp::Vp_StPtMtV(y,a,P,P_trans,M_trans,x,beta);
    TEUCHOS_TEST_FOR_EXCEPT(true);
    return;
  }

  if( M_trans == BLAS_Cpp::no_trans ) {
    //
    // y = beta*y + a * op(P) * A_bar * x
    //
    //   = beta * y
    //   
    //    + a * [op(P1)  op(P2) ] * [ I  0  op(E')  op(F') ] * [ x1 ]
    //                              [ 0  1   -b'     -f'   ]   [ x2 ]
    //                                                         [ x3 ]
    //                                                         [ x4 ]
    //
    // y = beta*y + a*op(P1)*x1 + a*op(P1)*op(E')*x3 + a*op(P1)*op(F')*x4
    //     + a*op(P2)*x2 - a*op(P2)*b'*x3 - a*op(P2)*f'*x4
    //
    // Where:
    //   op(P1) = op(P)(:,1:nd)
    //   op(P2) = op(P)(:,nd+1:nd+1)
    //

    const GenPermMatrixSlice
      P1 = ( P.is_identity() 
           ? GenPermMatrixSlice( nd(), nd(), GenPermMatrixSlice::IDENTITY_MATRIX )
           : P.create_submatrix(d_rng,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
        ),
      P2 = ( P.is_identity()
           ? GenPermMatrixSlice(
             P_trans == no_trans ? nd() : 1
             , P_trans == no_trans ? 1 : nd()
             , GenPermMatrixSlice::ZERO_MATRIX )
           : P.create_submatrix(Range1D(nd()+1,nd()+1),P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
        );

    const SpVectorSlice
      x1 = x(d_rng);
    const value_type
      x2 = get_sparse_element(x,nd()+1);
    const SpVectorSlice
      x3 = m_in() ? x(E_rng) : SpVectorSlice(NULL,0,0,0),
      x4 = m_eq() ? x(F_rng) : SpVectorSlice(NULL,0,0,0);

    // y = beta*y + a*op(P1)*x1
    Vp_StMtV( y, a, P1, P_trans, x1, beta );
    // y += a*op(P1)*op(E')*x3
    if( m_in() && P1.nz() )
      LinAlgOpPack::Vp_StPtMtV( y, a, P1, P_trans, *E(), trans_not(trans_E()), x3 );
    // y += a*op(P1)*op(F')*x4
    if( m_eq() && P1.nz() )
      LinAlgOpPack::Vp_StPtMtV( y, a, P1, P_trans, *F(), trans_not(trans_F()), x4 );
    //
    // y += a*op(P2)*x2 - a*op(P2)*b'*x3 - a*op(P2)*f'*x4
    //   += a * op(P2) * ( x2 + b'*x3 - f'*x4 )
    //   
    //   ==>
    //   
    // y(i) +=  a * ( x2 - b'*x3 - f'*x4 )
    //   
    if( P2.nz() ){
      TEUCHOS_TEST_FOR_EXCEPT( !( P2.nz() == 1 ) );
      const size_type
        i = P_trans == BLAS_Cpp::no_trans
          ? P2.begin()->row_i() : P2.begin()->col_j();
      value_type
        &y_i = (*y)(i) += a * x2;
      if(m_in())
        y_i += -a * dot(*this->b(),x3);
      if(m_eq())
        y_i += -a * dot(*this->f(),x4);
    }
  }
  else if ( M_trans == BLAS_Cpp::trans ) {
    //
    // y = beta*y + a*op(P)*A_bar'*x
    // 
    //   = beta*y
    //   
    //    + a * [ P1  P2  P3  P4 ] * [ I       0 ] * [ x1 ]
    //                               [ 0       1 ]   [ x2 ]
    //                               [ op(E)  -b ]
    //                               [ op(F)  -f ]
    //
    // y = beta*y + a*P1*x1 + a*P2*x2 + a*P3*op(E)*x1 - a*P3*b*x2
    //     + a*P4*op(F)*x1 - a*P4*f*x2
    //
    // Where:
    //   P1 = op(P)(:,1:nd)
    //   P2 = op(P)(:,nd+1:nd+1)
    //   P3 = op(P)(:,nd+1+1:nd+1+m_in)
    //   P4 = op(P)(:,nd+1+m_in+1:nd+1+m_in+m_eq)
    //

    TEUCHOS_TEST_FOR_EXCEPT( !(  !P.is_identity()  ) ); // We should never have this!

    size_type off = 0;
    const GenPermMatrixSlice
      P1 = P.create_submatrix(Range1D(off+1,off+nd())
                  ,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL);
    off += nd();
    const GenPermMatrixSlice
      P2 = P.create_submatrix(Range1D(off+1,off+1)
                  ,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL);
    off += 1;
    const GenPermMatrixSlice
      P3 = m_in()
        ? P.create_submatrix(Range1D(off+1,off+m_in())
                   ,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
        : GenPermMatrixSlice();
    off += m_in();
    const GenPermMatrixSlice
      P4 = m_eq()
        ? P.create_submatrix(Range1D(off+1,off+m_eq())
                   ,P_trans==trans?GPMSIP::BY_ROW:GPMSIP::BY_COL)
        : GenPermMatrixSlice();

    const SpVectorSlice
      x1 = x(d_rng);
    const value_type
      x2 = get_sparse_element(x,nd()+1);

    // y = beta*y + a*op(P1)*x1
    Vp_StMtV( y, a, P1, P_trans, x1, beta );
    // y += a*op(P2)*x2
    if( P2.nz() ){
      TEUCHOS_TEST_FOR_EXCEPT( !( P2.nz() == 1 ) );
      (*y)( P_trans == BLAS_Cpp::no_trans
          ? P2.begin()->row_i() : P2.begin()->col_j() )
        += a * x2;
    }
    if( m_in() && P3.nz() ) {
      // y += a*P3*op(E)*x1
      LinAlgOpPack::Vp_StPtMtV( y, a, P3, P_trans, *E(), trans_E(), x1 );
      // y += (-a*x2)*P3*b
      AbstractLinAlgPack::Vp_StMtV( y, - a * x2, P3, P_trans, *this->b() );
    }
    if( m_eq() && P4.nz() ) {
      // y += a*P4*op(F)*x1
      LinAlgOpPack::Vp_StPtMtV( y, a, P4, P_trans, *F(), trans_F(), x1 );
      // y += (-a*x2)*P4*f
      AbstractLinAlgPack::Vp_StMtV( y, - a * x2, P4, P_trans, *this->f() );
    }
    
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true);	// Invalid trans value
  }
*/
}

} // end namespace QPSchurPack 
} // end namespace ConstrainedOptPack 
