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

#include <fstream>
#include <algorithm>

#include "AbstractLinAlgPack_DirectSparseSolverDense.hpp"
#include "AbstractLinAlgPack_VectorDenseEncap.hpp"
#include "DenseLinAlgLAPack.hpp"
#include "DenseLinAlgPack_PermVecMat.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace {

// My swap function
template<class T>
inline
void my_swap( T* v1, T* v2 )
{
  T tmp = *v1;
  *v1 = *v2;
  *v2 = tmp;
}

// A cast to const is needed because the standard does not return a reference from
// valarray<>::operator[]() const.
template <class T>
std::valarray<T>& cva(const std::valarray<T>& va )
{
  return const_cast<std::valarray<T>&>(va);
}

} // end namespace

namespace AbstractLinAlgPack {

//
// Implementation of DirectSparseSolver(Imp) interface using dense LAPACK routines.
//

// //////////////////////////////////////////////////
// DirectSparseSolverDense::BasisMatrixDense

// Overridden from BasisMatrixImp

Teuchos::RCP<DirectSparseSolverImp::BasisMatrixImp>
DirectSparseSolverDense::BasisMatrixDense::create_matrix() const
{
  return Teuchos::rcp(new BasisMatrixDense);
}

void DirectSparseSolverDense::BasisMatrixDense::V_InvMtV(
  VectorMutable* y, BLAS_Cpp::Transp M_trans, const Vector& x
  ) const 
{
  using Teuchos::dyn_cast;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  size_type k;

  // Get concrete objects
  const FactorizationStructureDense
    &fs = dyn_cast<const FactorizationStructureDense>(*this->get_fact_struc());
  const FactorizationNonzerosDense
    &fn = dyn_cast<const FactorizationNonzerosDense>(*this->get_fact_nonzeros());

  VectorDenseMutableEncap  yd(*y);
  VectorDenseEncap         xd(x);

  TEUCHOS_TEST_FOR_EXCEPTION(
    yd().dim() != xd().dim(), std::invalid_argument
    ,"DirectSparseSolverDense::BasisMatrixDense::V_InvMtV(...) : Error, "
    " y.dim() = " << yd().dim() << " != x.dim() = " << xd().dim() << "!"
    );

  // Get temp storage for rhs and solution to communicate with xGESTRS
  Workspace<value_type>   B_store(wss,xd().dim());
  DMatrixSlice  B(&B_store[0],B_store.size(),B_store.size(),B_store.size(),1);

  //
  // Now we must permute the rhs or solution vectors based on our own
  // permutation fn.basis_perm_.
  //
  // xGETRF(...) factored the transpose of the basis matrix C' = Ct = P*L*U
  // where the permtuation P is stored in the array fn.basis_perm_ where
  //
  //    q = P * v
  //
  // is given by
  //
  //    q(i) = v(fn.basis_perm_(i)), for i = 1...n
  //
  // and q = P' * v is given by
  //
  //    q(fn.basis_perm_(i)) = v(i), for i = 1...n
  //
  // The system we are solving is therefore:
  //
  //   C * y = x   =>  U'*L'*P'*y = x
  //   
  //        for M_trans == no_trans
  //
  //   C'* y = x   =>  P*L*U*y = x   =>  L*U*y = P'*x 
  //
  //        for M_trans == trans
  // 

  // Copy rsh
  if( M_trans == BLAS_Cpp::trans && fn.rect_analyze_and_factor_ ) {
    // b = P'*x =
    DVectorSlice b = B.col(1);
//		DenseLinAlgPack::inv_perm_ele(xd(),fn.basis_perm_,&b);
    DenseLinAlgPack::perm_ele(xd(),fn.basis_perm_,&b);
  }
  else {
    B.col(1) = xd();
  }

  // Solve
  DenseLinAlgLAPack::getrs(
    fn.LU_(1,fs.rank_,1,fs.rank_), &cva(fn.ipiv_)[0], BLAS_Cpp::trans_not(M_trans)
    ,&B
    );

  // Copy solution
  if( M_trans == BLAS_Cpp::no_trans  && fn.rect_analyze_and_factor_ ) {
    // y = P*b = P*(P'*y)
    const DVectorSlice b = B.col(1);
//		DenseLinAlgPack::perm_ele(b,fn.basis_perm_,&yd());
    DenseLinAlgPack::inv_perm_ele(b,fn.basis_perm_,&yd());
  }
  else {
    yd() = B.col(1);
  }

}

// //////////////////////////////////////////////////
// DirectSparseSolverDense::FactorizationStructureDense

DirectSparseSolverDense::FactorizationStructureDense::FactorizationStructureDense()
{}

// //////////////////////////////////////////////////
// DirectSparseSolverDense

// Constructors/initializers

DirectSparseSolverDense::DirectSparseSolverDense()
{}

// Overridden from DirectSparseSolver

const DirectSparseSolver::basis_matrix_factory_ptr_t
DirectSparseSolverDense::basis_matrix_factory() const
{
  namespace mmp = MemMngPack;
  return Teuchos::rcp(new Teuchos::AbstractFactoryStd<BasisMatrix,BasisMatrixDense>());
}

void DirectSparseSolverDense::estimated_fillin_ratio(
  value_type estimated_fillin_ratio
  )
{
  // We ignore this!
}

// Overridden from DirectSparseSolverImp

const Teuchos::RCP<DirectSparseSolver::FactorizationStructure>
DirectSparseSolverDense::create_fact_struc() const
{
  return Teuchos::rcp(new FactorizationStructureDense);
}

const Teuchos::RCP<DirectSparseSolverImp::FactorizationNonzeros>
DirectSparseSolverDense::create_fact_nonzeros() const
{
  return Teuchos::rcp(new FactorizationNonzerosDense);
}

void DirectSparseSolverDense::imp_analyze_and_factor(
  const AbstractLinAlgPack::MatrixConvertToSparse   &A
  ,FactorizationStructure                         *fact_struc
  ,FactorizationNonzeros                          *fact_nonzeros
  ,DenseLinAlgPack::IVector                            *row_perm
  ,DenseLinAlgPack::IVector                            *col_perm
  ,size_type                                      *rank
  ,std::ostream                                   *out
  )
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  typedef MatrixConvertToSparse MCTS;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  if(out)
    *out << "\nUsing LAPACK xGETRF to analyze and factor a new matrix ...\n";

  // Get the concrete factorization and nonzeros objects
  FactorizationStructureDense
    &fs = dyn_cast<FactorizationStructureDense>(*fact_struc);
  FactorizationNonzerosDense
    &fn = dyn_cast<FactorizationNonzerosDense>(*fact_nonzeros);

  // Get the dimensions of things.
  const index_type
    m = A.rows(),
    n = A.cols(),
    nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

  // Validate input
  TEUCHOS_TEST_FOR_EXCEPTION(
    n <= 0 || m <= 0 || m > n, std::invalid_argument
    ,"DirectSparseSolverDense::imp_analyze_and_factor(...) : Error!" );

  // Extract the matrix in coordinate format
  Workspace<value_type>   a_val(wss,nz);
  Workspace<index_type>   a_row_i(wss,nz);
  Workspace<index_type>   a_col_j(wss,nz);
  A.coor_extract_nonzeros(
    MCTS::EXTRACT_FULL_MATRIX
    ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
    ,nz
    ,&a_val[0]
    ,nz
    ,&a_row_i[0]
    ,&a_col_j[0]
    );

  //
  // Fill the matrix LU = A' so that xGETRF will pivot by row to find
  // the basis.
  //
  // Here we will form the factor of A' = P*L*U where L will be
  // a n x m upper trapizodial matrix containing the factor lower
  // triangular factor in LU(1:rank,1:rank) and junk below this.
  //
  // Note that xGETRF() pivots by row so if any dependent columns
  // are found they will always be the last few columns.
  //

  // Resize the storage
  fn.LU_.resize(n,m);
  fn.ipiv_.resize(n);

  // Add in the nonzero entires transposed (allows for multiple entries with same
  // row and column indexes).
  fn.LU_ = 0.0;
  for( size_type k = 0; k < nz; ++k )
    fn.LU_(a_col_j[k],a_row_i[k]) += a_val[k];

  //
  // Have xGETRF factor this matrix.
  //

  DenseLinAlgLAPack::getrf( &fn.LU_(), &fn.ipiv_[0], &fs.rank_ );

  // Remember the dimensions
  fs.m_  = m;
  fs.n_  = n;
  fs.nz_ = nz;

  //
  // At this point it is important to understand exactly what
  // ipiv() represents.  Each entry in ipiv(k) represents a row
  // interchange A(k) <=> A(ipiv(k)).  Therefore, we have to
  // do the same row interchanges to the identity permutation
  // of col_perm to form the permutation array expected by the
  // DSS interface.
  //

  // Form fs.col_perm_
  fs.col_perm_.resize(n);
  DenseLinAlgPack::identity_perm(&fs.col_perm_);
  Workspace<index_type> col_perm_unsorted(wss,fs.rank_);
  if( m == n && n == fs.rank_ ) {
    // Leave fs.col_perm_ as identity
    fn.rect_analyze_and_factor_ = false;
  }
  else {
    fn.rect_analyze_and_factor_ = true;
    // Permute fs.col_perm_ and save them in col_perm_unsorted
    for( index_type k = 1; k <= fs.rank_; ++k ) {
      my_swap( &fs.col_perm_(k), &fs.col_perm_(fn.ipiv_[k-1]) );
      col_perm_unsorted[k-1] = fs.col_perm_(k);
    }
    // Sort the basis selection
    std::sort(&(fs.col_perm_)[0]           , &(fs.col_perm_)[0] + fs.rank_ );
    std::sort(&(fs.col_perm_)[0] + fs.rank_, &(fs.col_perm_)[0] + n        );
  }

  // Form the inverse permutation
  fs.inv_col_perm_.resize(n);
  DenseLinAlgPack::inv_perm( fs.col_perm_, &fs.inv_col_perm_ );

  if( !(m == n && n == fs.rank_) ) {
    // Form fn.basis_perm_ and set fs.ipiv_ to identity
    fn.basis_perm_.resize(fs.rank_);
    for( size_type k = 1; k <= fs.rank_; ++k ) {
      fn.basis_perm_(k) = fs.inv_col_perm_(col_perm_unsorted[k-1]);
      fn.ipiv_[k-1] = k;
    }
  }

  // Copy the data to the output
  row_perm->resize(m);
  col_perm->resize(n);
  *rank = fs.rank_;
  DenseLinAlgPack::identity_perm(row_perm);
  *col_perm = fs.col_perm_;

}

void DirectSparseSolverDense::imp_factor(
  const AbstractLinAlgPack::MatrixConvertToSparse   &A
  ,const FactorizationStructure                   &fact_struc
  ,FactorizationNonzeros                          *fact_nonzeros
  ,std::ostream                                   *out
  )
{
  namespace mmp = MemMngPack;
  using Teuchos::dyn_cast;
  typedef MatrixConvertToSparse MCTS;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  if(out)
    *out << "\nUsing LAPACK xGETRF to refactor the basis matrix ...\n";

  // Get the concrete factorization and nonzeros objects
  const FactorizationStructureDense
    &fs = dyn_cast<const FactorizationStructureDense>(fact_struc);
  FactorizationNonzerosDense
    &fn = dyn_cast<FactorizationNonzerosDense>(*fact_nonzeros);

  // Get the dimensions of things.
  const index_type
    m = A.rows(),
    n = A.cols(),
    nz = A.num_nonzeros( MCTS::EXTRACT_FULL_MATRIX ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM );

  // Validate input
  TEUCHOS_TEST_FOR_EXCEPTION(
    (m != fs.m_ || n != fs.n_ || nz != fs.nz_), std::invalid_argument
    ,"DirectSparseSolverDense::imp_factor(...): Error!"
    );

  // Extract the matrix in coordinate format
  Workspace<value_type>   a_val(wss,nz);
  Workspace<index_type>   a_row_i(wss,nz);
  Workspace<index_type>   a_col_j(wss,nz);
  A.coor_extract_nonzeros(
    MCTS::EXTRACT_FULL_MATRIX
    ,MCTS::ELEMENTS_ALLOW_DUPLICATES_SUM
    ,nz
    ,&a_val[0]
    ,nz
    ,&a_row_i[0]
    ,&a_col_j[0]
    );
  
  //
  // Fill the matrix LU = B so that xGETRF will pivot by row to find
  // the basis.  Here B is the basis matrix from A'.
  //
  // Here we will form the factor of B = P*L*U where L will be
  // a rank x rank lower triangular.
  //

  // Resize the storage
  fn.rect_analyze_and_factor_ = false;
  fn.LU_.resize(fs.rank_,fs.rank_);
  fn.ipiv_.resize(fs.rank_);

  // Copy only the basis entries (transposed)
  fn.LU_ = 0.0;
  for( size_type k = 0; k < nz; ++k ) {
    const index_type B_i = fs.inv_col_perm_(a_col_j[k]);
    const index_type B_j = a_row_i[k];
    if( B_i <= fs.rank_ && B_j <= fs.rank_ )
      fn.LU_(B_i,B_j) = a_val[k];
  }

  //
  // Have xGETRF factor this matrix.
  //

  FortranTypes::f_int B_rank = 0;
  DenseLinAlgLAPack::getrf( &fn.LU_(), &fn.ipiv_[0], &B_rank );

  TEUCHOS_TEST_FOR_EXCEPTION(
    B_rank != fs.rank_, FactorizationFailure
    ,"DirectSparseSolverDense::imp_factor(...): Error, the basis matrix is no "
    "longer full rank with B_rank = " << B_rank << " != fs.rank = " << fs.rank_ << "!"
    );

}

// private

}	// end namespace AbstractLinAlgPack 
