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

#include <assert.h>

#include <limits>

#include "ConstrainedOptPack_MatrixSymAddDelBunchKaufman.hpp"
#include "DenseLinAlgLAPack.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_DMatrixOp.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "DenseLinAlgPack_delete_row_col.hpp"

namespace ConstrainedOptPack {

MatrixSymAddDelBunchKaufman::MatrixSymAddDelBunchKaufman()
  : S_size_(0), S_indef_(false), fact_updated_(false), fact_in1_(true), inertia_(0,0,0)
{}

void MatrixSymAddDelBunchKaufman::pivot_tols( PivotTolerances pivot_tols )
{
  S_chol_.pivot_tols(pivot_tols);
}

MatrixSymAddDelUpdateable::PivotTolerances MatrixSymAddDelBunchKaufman::pivot_tols() const
{
  return S_chol_.pivot_tols();
}

// Overridden from MatrixSymAddDelUpdateableWithOpNonsingular

const MatrixSymOpNonsing& MatrixSymAddDelBunchKaufman::op_interface() const
{
  return *this;
}

MatrixSymAddDelUpdateable& MatrixSymAddDelBunchKaufman::update_interface()
{
  return *this;
}

const MatrixSymAddDelUpdateable& MatrixSymAddDelBunchKaufman::update_interface() const
{
  return *this;
}

// Overridden from MatrixSymAddDelUpdateable

void MatrixSymAddDelBunchKaufman::initialize(
  value_type         alpha
  ,size_type         max_size
  )
{
  try {
    // Resize the storage if we have to
    if( S_store1_.rows() < max_size+1 && S_store1_.cols() < max_size+1 )
      S_store1_.resize(max_size+1,max_size+1);
    fact_in1_ = true;
    // Start out with a p.d. or n.d. matrix and maintain the original
    S_chol_.init_setup(&S_store1_(),Teuchos::null,0,true,true,true,false,0.0);
    S_chol_.initialize(alpha,max_size);
    // Set the state variables:
    S_size_  = 1;
    S_indef_ = false; // fact_updated_, fact_in1_ and inertia are meaningless!
  }
  catch(...) {
    S_size_ = 0;
    throw;
  }
}

void MatrixSymAddDelBunchKaufman::initialize(
  const DMatrixSliceSym      &A
  ,size_type         max_size
  ,bool              force_factorization
  ,Inertia           expected_inertia
  ,PivotTolerances   pivot_tols
  )
{
  using BLAS_Cpp::upper;
  using BLAS_Cpp::lower;
  using DenseLinAlgPack::tri_ele;
  using DenseLinAlgPack::nonconst_tri_ele;
  typedef MatrixSymAddDelUpdateable  MSADU;
  typedef MSADU::Inertia   Inertia;

  bool                throw_exception = false; // If true then throw exception
  std::ostringstream  omsg;                    // Will be set if an exception has to be thrown.
  value_type          gamma;                   // ...

  const size_type
    n = A.rows();

  // Validate proper usage of inertia parameter
  TEUCHOS_TEST_FOR_EXCEPT( ! ( expected_inertia.zero_eigens == Inertia::UNKNOWN
    || expected_inertia.zero_eigens == 0 ) );
  
  try {
    // Resize the storage if we have to
    if( S_store1_.rows() < max_size+1 && S_store1_.cols() < max_size+1 )
      S_store1_.resize(max_size+1,max_size+1);
    fact_in1_ = true;
    // See if the client claims that the matrix is p.d. or n.d.
    const bool not_indefinite
      = ( ( expected_inertia.neg_eigens == 0 && expected_inertia.pos_eigens == n )
        || ( expected_inertia.neg_eigens == n && expected_inertia.pos_eigens == 0 ) );
    // Initialize the matrix
    if( not_indefinite ) {
      // The client says that the matrix is p.d. or n.d. so
      // we will take their word for it.
      S_chol_.init_setup(&S_store1_(),Teuchos::null,0,true,true,true,false,0.0);
      try {
        S_chol_.initialize(A,max_size,force_factorization,expected_inertia,pivot_tols);
      }
      catch(const MSADU::WarnNearSingularUpdateException& except) {
        throw_exception = true; // Throw this same exception at the end!
        omsg << except.what();
        gamma = except.gamma;
      }
      // Set the state variables:
      S_size_       = n;
      S_indef_      = false; // fact_updated_, fact_in1_ and inertia are meaningless!
    }
    else {
      //
      // The client did not say that the matrix was p.d. or n.d. so we
      // must assume that it is indefinite.
      //
      bool fact_in1 = !fact_in1_;
      // Set the new matrix in the unused factorization location
      DenseLinAlgPack::assign(	&DU(n,fact_in1), tri_ele( A.gms(), A.uplo() ) );
      // Update the factorization in place
      try {
        factor_matrix( n, fact_in1 );
      }
      catch( const DenseLinAlgLAPack::FactorizationException& excpt ) {
        omsg
          << "MatrixSymAddDelBunchKaufman::initialize(...): "
          << "Error, the matrix A is singular:\n"
          << excpt.what();
        throw SingularUpdateException( omsg.str(), -1.0 );
      }
      // Compute the inertia and validate that it is correct.
      Inertia inertia;
      throw_exception = compute_assert_inertia(
        n,fact_in1,expected_inertia,"initialize",pivot_tols
        ,&inertia,&omsg,&gamma);
      // If the client did not know the inertia of the
      // matrix but it turns out to be p.d. or n.d. then modify the
      // DU factor appropriatly and switch to cholesky factorization.
      if( inertia.pos_eigens == n || inertia.neg_eigens == n ) {
        TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo Implememnt this!
      }
      else {
        // Indead the new matrix is indefinite
        // Set the original matrix now
        DenseLinAlgPack::assign(	&DenseLinAlgPack::nonconst_tri_ele( S(n).gms(), BLAS_Cpp::lower)
                  , tri_ele( A.gms(), A.uplo() ) );
        // Update the state variables:
        S_size_       = n;
        S_indef_      = true;
        fact_updated_ = true;
        fact_in1_     = fact_in1;
        inertia_      = inertia;
      }
    }
  }
  catch(...) {
    S_size_ = 0;
    throw;
  }
  if( throw_exception )
    throw WarnNearSingularUpdateException(omsg.str(),gamma);
}

size_type MatrixSymAddDelBunchKaufman::max_size() const
{
  return S_store1_.rows() -1;  // The center diagonal is used for workspace
}

MatrixSymAddDelUpdateable::Inertia
MatrixSymAddDelBunchKaufman::inertia() const
{
  return S_indef_ ? inertia_ : S_chol_.inertia();
}

void MatrixSymAddDelBunchKaufman::set_uninitialized()
{
  S_size_ = 0;
}

void MatrixSymAddDelBunchKaufman::augment_update(
  const DVectorSlice  *t
  ,value_type        alpha
  ,bool              force_refactorization
  ,EEigenValType     add_eigen_val
  ,PivotTolerances   pivot_tols
  )
{
  using Teuchos::implicit_ref_cast;
  using BLAS_Cpp::no_trans;
  using DenseLinAlgPack::norm_inf;
  using AbstractLinAlgPack::transVtInvMtV;
  typedef MatrixSymAddDelUpdateable  MSADU;

  assert_initialized();

  // Validate the input
  if (implicit_ref_cast<const MatrixBase>(*this).rows() >= max_size()) {
    throw MaxSizeExceededException(
      "MatrixSymAddDelBunchKaufman::augment_update(...) : "
      "Error, the maximum size would be exceeded." );
  }
  if (t && t->dim() != S_size_) {
    throw std::length_error(
      "MatrixSymAddDelBunchKaufman::augment_update(...): "
      "Error, t.dim() must be equal to this->rows()." );
  }
  if (!(add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN
      || add_eigen_val != MSADU::EIGEN_VAL_ZERO))
  {
    throw SingularUpdateException(
      "MatrixSymAddDelBunchKaufman::augment_update(...): "
      "Error, the client has specified a singular update in add_eigen_val.", -1.0 );
  }

  // Try to perform the update
  bool throw_exception = false; // If true then throw exception
  std::ostringstream omsg;      // Will be set if an exception has to be thrown.
  value_type gamma;             // ...
  if( !S_indef_ ) {
    // The current matrix is positive definite or negative definite.  If the
    // new matrix is still p.d. or n.d. when we are good to go.  We must first
    // check if the client has specified that the new matrix will be indefinite
    // and if so then we will take his/her word for it and do the indefinite
    // updating.
    MSADU::Inertia inertia = S_chol_.inertia();
    const bool
      new_S_not_indef
      = ( ( inertia.pos_eigens > 0
        && ( add_eigen_val == MSADU::EIGEN_VAL_POS || add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN ) )
        || ( inertia.neg_eigens > 0
        && ( add_eigen_val == MSADU::EIGEN_VAL_NEG || add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN ) )
        );
    if( !new_S_not_indef ) {
      // We must resort to an indefinite factorization
    }
    else {
      // The new matrix could/should still be p.d. or n.d. so let's try it!
      bool update_successful = false;
      try {
        S_chol_.augment_update(t,alpha,force_refactorization,add_eigen_val,pivot_tols);
        update_successful = true;
      }
      catch(const MSADU::WrongInertiaUpdateException&) {
        if( add_eigen_val != MSADU::EIGEN_VAL_UNKNOWN )
          throw; // The client specified the new inertia and it was wrong so throw execepiton.
        // If the client did not know that inertia then we can't fault them so we will
        // proceed unto the indefinite factorization.
      }
      catch(const MSADU::WarnNearSingularUpdateException& except) {
        throw_exception = true; // Throw this same exception at the end!
        update_successful = true;
        omsg << except.what();
        gamma = except.gamma;
      }
      if( update_successful ) {
        ++S_size_;     // Now we can resize the matrix
        if(throw_exception)
          throw MSADU::WarnNearSingularUpdateException(omsg.str(),gamma);
        else
          return;
      }
    }
  }
  //
  // If we get here then we have no choice but the perform an indefinite factorization.
  //
  // ToDo: (2/2/01): We could also try some more fancy updating
  // procedures and try to get away with some potentially unstable
  // methods in the interest of time savings.
  //
  const size_type n     = S_size_;
  DMatrixSlice  S_bar = this->S(n+1).gms();
  //
  // Validate that the new matrix will be nonsingular.
  //
  // Given any nonsingular matrix S (even unsymmetric) it is easy to show that
  // beta = alpha - t'*inv(S)*t != 0.0 if [ S, t; t', alpha ] is nonsingular.
  // This is a cheap O(n^2) way to check that the update is nonsingular before
  // we go through the O(n^3) refactorization.
  // In fact, the sign of beta even tells us the sign of the new eigen value
  // of the updated matrix even before we perform the refactorization.
  // If the current matrix is not factored then we will just skip this
  // test and let the full factorization take place to find this out.
  // We could even cheat a little and actually perform the update with this
  // diagonal beta and then compute the update to the factors U our selves
  // in O(n^2) time.
  //
  if( force_refactorization && fact_updated_ ) {
    const value_type
      beta        = alpha - ( t ? transVtInvMtV(*t,*this,no_trans,*t) : 0.0 ),
      abs_beta    = std::fabs(beta),
      nrm_D_diag  = norm_inf(DU(n,fact_in1_).gms().diag()); // ToDo: Consider 2x2 blocks also!
    gamma = abs_beta / nrm_D_diag;
    // Check gamma
    const bool
      correct_inertia = (
        add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN
        || beta > 0.0 && add_eigen_val == MSADU::EIGEN_VAL_POS
        || beta < 0.0 && add_eigen_val == MSADU::EIGEN_VAL_NEG
        );
    PivotTolerances use_pivot_tols = S_chol_.pivot_tols();
    if( pivot_tols.warning_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.warning_tol = pivot_tols.warning_tol;
    if( pivot_tols.singular_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.singular_tol = pivot_tols.singular_tol;
    if( pivot_tols.wrong_inertia_tol != PivotTolerances::UNKNOWN )
      use_pivot_tols.wrong_inertia_tol = pivot_tols.wrong_inertia_tol;
    throw_exception = (
      gamma == 0.0
      || correct_inertia && gamma <= use_pivot_tols.singular_tol
      || !correct_inertia
      );
    // Create message and throw exceptions
    std::ostringstream onum_msg;
    if(throw_exception) {
      onum_msg
        << "gamma = |alpha-t'*inv(S)*t)|/||diag(D)||inf = |" <<  beta << "|/" << nrm_D_diag
        << " = " << gamma;
      omsg
        << "MatrixSymAddDelBunchKaufman::augment_update(...): ";
      if( correct_inertia && gamma <= use_pivot_tols.singular_tol ) {
        omsg
          << "Singular update!\n" << onum_msg.str() << " <= singular_tol = "
          << use_pivot_tols.singular_tol;
        throw SingularUpdateException( omsg.str(), gamma );
      }
      else if( !correct_inertia && gamma <= use_pivot_tols.singular_tol ) {
        omsg
          << "Singular update!\n" << onum_msg.str() << " <= wrong_inertia_tol = "
          << use_pivot_tols.wrong_inertia_tol;
        throw SingularUpdateException( omsg.str(), gamma );
      }
      
      else if( !correct_inertia ) {
        omsg
          << "Indefinite update!\n" << onum_msg.str() << " >= wrong_inertia_tol = "
          << use_pivot_tols.wrong_inertia_tol;
        throw WrongInertiaUpdateException( omsg.str(), gamma );
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPT(true); // Only local programming error?
      }
    }
  }

  // Add new row to the lower part of the original symmetric matrix.
  if(t)
    S_bar.row(n+1)(1,n) = *t;
  else
    S_bar.row(n+1)(1,n) = 0.0;
  S_bar(n+1,n+1) = alpha;

  //
  // Update the factorization
  //
  if( force_refactorization ) {
    // Determine where to copy the original matrix to
    bool fact_in1 = false;
    if( S_indef_ ) {
      // S is already indefinite so let's copy the original into the storage
      // location other than the current one in case the newly nonsingular matrix
      // is singular or has the wrong inertia.
      fact_in1 = !fact_in1_;
    }
    else {
      // S is currently p.d. or n.d. so let's copy the new matrix
      // into the second storage location so as not to overwrite
      // the current cholesky factor in case the new matrix is singular
      // or has the wrong inertia.
      fact_in1 = false;
    }
    // Copy and factor the new matrix
    try {
      copy_and_factor_matrix(n+1,fact_in1);
    }
    catch( const DenseLinAlgLAPack::FactorizationException& excpt ) {
      std::ostringstream omsg;
      omsg
        << "MatrixSymAddDelBunchKaufman::augment_update(...): "
        << "Error, singular update but the original matrix was be maintianed:\n"
        << excpt.what();
      throw SingularUpdateException( omsg.str(), -1.0 );
    }
    // Compute the expected inertia
    Inertia expected_inertia = this->inertia();
    if( expected_inertia.zero_eigens == Inertia::UNKNOWN || add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN )
      expected_inertia = Inertia(); // Unknow inertia!
    else if( add_eigen_val == MSADU::EIGEN_VAL_NEG )
      ++expected_inertia.neg_eigens;
    else if( add_eigen_val == MSADU::EIGEN_VAL_POS )
      ++expected_inertia.pos_eigens;
    else
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should not happen!
    // Compute the actually inertia and validate that it is what is expected
    Inertia inertia;
    throw_exception = compute_assert_inertia(
      n+1,fact_in1,expected_inertia,"augment_update",pivot_tols
      ,&inertia,&omsg,&gamma);
    // Unset S_chol so that there is no chance of accedental modification.
    if(!S_indef_)
      S_chol_.init_setup(NULL);
    // Update the state variables
    ++S_size_;
    S_indef_      = true;
    fact_updated_ = true;
    fact_in1_     = fact_in1;
    inertia_      = inertia;
  }
  else {
    //
    // Don't update the factorization yet
    //
    if(!S_indef_) // We must set the inertia if it was definite!
      inertia_ = S_chol_.inertia();
    ++S_size_;
    S_indef_      = true;
    fact_updated_ = false;  // fact_in1_ is meaningless now
    // We need to keep track of the expected inertia!
    if( inertia_.zero_eigens == Inertia::UNKNOWN || add_eigen_val == MSADU::EIGEN_VAL_UNKNOWN )
      inertia_ = Inertia(); // Unknow inertia!
    else if( add_eigen_val == MSADU::EIGEN_VAL_NEG )
      ++inertia_.neg_eigens;
    else if( add_eigen_val == MSADU::EIGEN_VAL_POS )
      ++inertia_.pos_eigens;
    else
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should not happen!
  }
  if( throw_exception )
    throw WarnNearSingularUpdateException(omsg.str(),gamma);
}

void MatrixSymAddDelBunchKaufman::delete_update(
  size_type          jd
  ,bool              force_refactorization
  ,EEigenValType     drop_eigen_val
  ,PivotTolerances   pivot_tols
  )
{
  using BLAS_Cpp::upper;
  using BLAS_Cpp::lower;
  using DenseLinAlgPack::tri_ele;
  using DenseLinAlgPack::nonconst_tri_ele;
  typedef MatrixSymAddDelUpdateable  MSADU;

  assert_initialized();

  if( jd < 1 || S_size_ < jd )
    throw std::out_of_range(
      "MatrixSymAddDelBunchKaufman::delete_update(jd,...): "
      "Error, the indice jd must be 1 <= jd <= rows()" );

  bool                throw_exception = false; // If true then throw exception
  std::ostringstream  omsg;                    // Will be set if an exception has to be thrown.
  value_type          gamma;                   // ...

  if( !S_indef_ ) {
    //
    // The current matrix is p.d. or n.d. and so should the new matrix.
    //
    S_chol_.delete_update(jd,force_refactorization,drop_eigen_val,pivot_tols);
    --S_size_;
  }
  else {
    //
    // The current matrix is indefinite but the new matrix
    // could be p.d. or n.d. in which case we could switch
    // to the cholesky factorization.  If the user says the
    // new matrix will be p.d. or n.d. then we will take
    // his/her word for it and try the cholesky factorization
    // otherwise just recompute the indefinite factorization.
    //
    // ToDo: (2/2/01): We could also try some more fancy updating
    // procedures and try to get away with some potentially unstable
    // methods and resort to the full indefinite factorization if needed.
    // The sign of the dropped eigen value might tell us what we
    // might get away with?
    //
    // Update the factorization
    //
    Inertia inertia = S_chol_.inertia();
    if( (drop_eigen_val == MSADU::EIGEN_VAL_POS && inertia.pos_eigens == 1 )
      || (drop_eigen_val == MSADU::EIGEN_VAL_NEG && inertia.neg_eigens == 1 )
      || S_size_ == 2 )
    {
      // Lets take the users word for it and switch to a cholesky factorization.
      // To do this we will just let S_chol_ do it for us.  If the new matrix
      // turns out not to be what the user says, then we will resize the matrix
      // to zero and thrown an exception.
      try {
        // Delete row and column jd from M
        DMatrixSlice S = this->S(S_size_).gms();
        DenseLinAlgPack::delete_row_col( jd, &DenseLinAlgPack::nonconst_tri_ele(S,BLAS_Cpp::lower) ); 
        // Setup S_chol and factor the thing
        S_chol_.init_setup(&S_store1_(),Teuchos::null,0,true,true,true,false,0.0);
        S_chol_.initialize( this->S(S_size_-1), S_store1_.rows()-1
                  , force_refactorization, Inertia(), PivotTolerances() );
      }
      catch( const SingularUpdateException& excpt ) {
        S_size_ = 0;
        throw SingularUpdateException(
          std::string(
            "MatrixSymAddDelBunchKaufman::delete_update(...) : "
            "Error, the client incorrectly specified that the new "
            "matrix would be nonsingular and not indefinite:\n" )
          + excpt.what()
          , excpt.gamma
          );
      }
      catch( const WrongInertiaUpdateException& excpt ) {
        S_size_ = 0;
        throw WrongInertiaUpdateException(
          std::string(
            "MatrixSymAddDelBunchKaufman::delete_update(...) : "
            "Error, the client incorrectly specified that the new "
            "matrix would not be indefinite:\n" )
          + excpt.what()
          , excpt.gamma
          );
      }
      // If we get here the update succeeded and the new matrix is p.d. or n.d.
      --S_size_;
      S_indef_ = false;
    }
    else {
      // 
      // We have been given no indication that the new matrix is p.d. or n.d.
      // so we will assume it is indefinite and go from there.
      //
      DMatrixSlice  S  = this->S(S_size_).gms();
      if( force_refactorization ) {
        // 
        // Perform the refactorization carefully
        //
        const bool fact_in1 = !fact_in1_;
        // Copy the original into the unused storage location
        DMatrixSliceTriEle  DU = this->DU(S_size_,fact_in1);
        DenseLinAlgPack::assign(	&DU, tri_ele(S,lower) );
        // Delete row and column jd from the storage location for DU
        DenseLinAlgPack::delete_row_col( jd, &DU );
        // Now factor the matrix inplace
        try {
          factor_matrix(S_size_-1,fact_in1);
        }
        catch( const DenseLinAlgLAPack::FactorizationException& excpt ) {
          omsg
            << "MatrixSymAddDelBunchKaufman::delete_update(...): "
            << "Error, singular update but the original matrix was maintianed:\n"
            << excpt.what();
          throw SingularUpdateException( omsg.str(), -1.0 );
        }
        // Compute the expected inertia
        Inertia expected_inertia = this->inertia();
        if( expected_inertia.zero_eigens == Inertia::UNKNOWN || drop_eigen_val == MSADU::EIGEN_VAL_UNKNOWN )
          expected_inertia = Inertia(); // Unknow inertia!
        else if( drop_eigen_val == MSADU::EIGEN_VAL_NEG )
          --expected_inertia.neg_eigens;
        else if( drop_eigen_val == MSADU::EIGEN_VAL_POS )
          --expected_inertia.pos_eigens;
        else
          TEUCHOS_TEST_FOR_EXCEPT(true); // Should not happen!
        // Compute the exacted inertia and validate that it is what is expected
        Inertia inertia;
        throw_exception = compute_assert_inertia(
          S_size_-1,fact_in1,expected_inertia,"delete_update",pivot_tols
          ,&inertia,&omsg,&gamma);
        // If we get here the factorization worked out.
        --S_size_;
        S_indef_      = true;
        fact_updated_ = true;
        fact_in1_     = fact_in1;
        inertia_      = inertia;
      }
      // Delete the row and column jd from the original
      DenseLinAlgPack::delete_row_col( jd, &nonconst_tri_ele(S,lower) );
      if( !force_refactorization ) {
        //
        // The refactorization was not forced
        //
        --S_size_;
        S_indef_      = true;
        fact_updated_ = false;  // fact_in1_ is meaningless now
        // We need to keep track of the expected inertia!
        if( inertia_.zero_eigens == Inertia::UNKNOWN || drop_eigen_val == MSADU::EIGEN_VAL_UNKNOWN )
          inertia_ = Inertia(); // Unknow inertia!
        else if( drop_eigen_val == MSADU::EIGEN_VAL_NEG )
          --inertia_.neg_eigens;
        else if( drop_eigen_val == MSADU::EIGEN_VAL_POS )
          --inertia_.pos_eigens;
        else
          TEUCHOS_TEST_FOR_EXCEPT(true); // Should not happen!
      }
    }
  }
  if( throw_exception )
    throw WarnNearSingularUpdateException(omsg.str(),gamma);
}

// Overridden from MatrixSymOpNonsingSerial

size_type MatrixSymAddDelBunchKaufman::rows() const
{
  return S_size_;
}

std::ostream& MatrixSymAddDelBunchKaufman::output(std::ostream& out) const
{
  if( S_size_ ) {
    out << "Unfactored symmetric matrix stored as lower triangle (ignore upper nonzeros):\n"
      << S(S_size_).gms();
    if( S_indef_ ) {
      out << "Upper symmetric indefinite factor DU returned from sytrf(...) (ignore lower nonzeros):\n"
        << DU(S_size_,fact_in1_).gms();
      out << "Permutation array IPIV returned from sytrf(...):\n";
      for( size_type i = 0; i < S_size_; ++i )
        out << "  " << IPIV_[i];
      out << std::endl;
    }
    else {
      out << "Upper cholesky factor (ignore lower nonzeros):\n" << DU(S_size_,true).gms();
    }
  }
  else {
    out << "0 0\n";
  }
  return out;
}

void MatrixSymAddDelBunchKaufman::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const DVectorSlice& x, value_type b
  ) const
{
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(), S_size_, S_size_, M_trans, x.dim() );
  // Use the unfactored matrix since this is more accurate!
  DenseLinAlgPack::Vp_StMtV( y, a, S(S_size_), BLAS_Cpp::no_trans, x, b );
}

void MatrixSymAddDelBunchKaufman::V_InvMtV(
  DVectorSlice* y, BLAS_Cpp::Transp M_trans
  ,const DVectorSlice& x
  ) const
{
  const size_type n = S_size_;
  DenseLinAlgPack::Vp_MtV_assert_sizes( y->dim(), n, n, M_trans, x.dim() );
  if( S_indef_ ) {
    // Update the factorzation if needed
    if(!fact_updated_) {
      const bool fact_in1 = true;
      MatrixSymAddDelBunchKaufman
        *nc_this = const_cast<MatrixSymAddDelBunchKaufman*>(this);
      nc_this->copy_and_factor_matrix(S_size_,fact_in1); // May throw exceptions!
      nc_this->fact_updated_ = true;
      nc_this->fact_in1_     = fact_in1;
    }
    *y = x;
    DenseLinAlgLAPack::sytrs(
      DU(S_size_,fact_in1_), &const_cast<IPIV_t&>(IPIV_)[0]
      , &DMatrixSlice(y->raw_ptr(),n,n,n,1), &WORK_() );
  }
  else {
    AbstractLinAlgPack::V_InvMtV( y, S_chol_, M_trans, x );
  }
}

// Private

void MatrixSymAddDelBunchKaufman::assert_initialized() const
{
  if(!S_size_)
    throw std::logic_error(
      "MatrixSymAddDelBunchKaufman::assert_initialized() : "
      "Error, this matrix is not initialized yet" );
}

void MatrixSymAddDelBunchKaufman::resize_DU_store( bool in_store1 )
{
  if(!in_store1 && S_store2_.rows() < S_store1_.rows() )
    S_store2_.resize( S_store1_.rows(), S_store1_.cols() );
}

void MatrixSymAddDelBunchKaufman::copy_and_factor_matrix(
  size_type S_size, bool fact_in1 )
{
  DenseLinAlgPack::assign(
    &DU(S_size,fact_in1)
    , DenseLinAlgPack::tri_ele(S(S_size).gms(),BLAS_Cpp::lower) );
  factor_matrix(S_size,fact_in1);
}

void MatrixSymAddDelBunchKaufman::factor_matrix( size_type S_size, bool fact_in1 )
{
  // Resize the workspace first
  const size_type opt_block_size = 64; // This is an estimate (see sytrf(...))
  if( WORK_.dim() < S_store1_.rows() * opt_block_size )
    WORK_.resize(S_store1_.rows()*opt_block_size);
  if( IPIV_.size() < S_store1_.rows() )
    IPIV_.resize(S_store1_.rows());
  // Factor the matrix (will throw FactorizationException if singular)
  DenseLinAlgLAPack::sytrf( &DU(S_size,fact_in1), &IPIV_[0], &WORK_() );
}

bool MatrixSymAddDelBunchKaufman::compute_assert_inertia(
  size_type S_size, bool fact_in1, const Inertia& exp_inertia, const char func_name[]
  ,PivotTolerances pivot_tols, Inertia* comp_inertia, std::ostringstream* omsg, value_type* gamma )
{
  bool throw_exception = false;
  *gamma = 0.0;
  // Here we will compute the inertia given IPIV[] and D[] as described in the documentation
  // for dsytrf(...) (see the source code).
  const DMatrixSlice DU = this->DU(S_size,fact_in1).gms();
  const size_type      n = DU.rows();
  Inertia inertia(0,0,0);
  value_type max_diag = 0.0;
  value_type min_diag = std::numeric_limits<value_type>::max();
  for( size_type k = 1; k <= n; ) {
    const FortranTypes::f_int k_p = IPIV_[k-1];
    if( k_p > 0 ) {
      // D(k,k) is a 1x1 block.
      // Lets get the eigen value from the sign of D(k,k)
      const value_type D_k_k = DU(k,k), abs_D_k_k = std::fabs(D_k_k);
      if( D_k_k > 0.0 )
        ++inertia.pos_eigens;
      else
        ++inertia.neg_eigens;
      if(abs_D_k_k > max_diag) max_diag = abs_D_k_k;
      if(abs_D_k_k < min_diag) min_diag = abs_D_k_k;
      k++;
    }
    else {
      // D(k:k+1,k:k+1) is a 2x2 block.
      // This represents one positive eigen value and
      // on negative eigen value
      TEUCHOS_TEST_FOR_EXCEPT( !(  IPIV_[k] == k_p ) ); // This is what the documentation for xSYTRF(...) says!
      ++inertia.pos_eigens;
      ++inertia.neg_eigens;
      // To find the largest and smallest diagonals of U for L*U we must perform Gaussian
      // elimination on this 2x2 block:
      const value_type                                   // [ a   b ] = D(k:k+1,k:k+1) 
        a = DU(k,k), b = DU(k,k+1), c = DU(k+1,k+1),   // [ b   c ]
        abs_a = std::fabs(a), abs_b = std::fabs(b);  
      value_type pivot_1, pivot_2;
      if( abs_a > abs_b ) { // Pivot on a = D(k,k)
        pivot_1 = abs_a;              // [   1      ] * [ a   b ] = [ a      b     ]
        pivot_2 = std::fabs(c - b*b/a);  // [ -b/a  1  ]   [ b   c ]   [ 0  c - b*b/a ]
      }
      else {                // Pivot on b = D(k+1,k) = D(k,k+1)
        pivot_1 = abs_b;              // [   1      ] * [ b   c ] = [ b      c     ]
        pivot_2 = std::fabs(b - a*c/b);  // [ -a/b  1  ]   [ a   b ]   [ 0  b - a*c/b ]
      }
      if(pivot_1 > max_diag) max_diag = pivot_1;
      if(pivot_1 < min_diag) min_diag = pivot_1;
      if(pivot_2 > max_diag) max_diag = pivot_2;
      if(pivot_2 < min_diag) min_diag = pivot_2;
      k+=2;
    }
  }
  // gamma = min(|diag(i)|)/max(|diag(i)|)
  *gamma = min_diag / max_diag;
  // Now validate that the inertia is what is expected
    const bool
    wrong_inertia =
    ( exp_inertia.neg_eigens != Inertia::UNKNOWN
      && exp_inertia.neg_eigens != inertia.neg_eigens )
    || ( exp_inertia.pos_eigens != Inertia::UNKNOWN
       && exp_inertia.pos_eigens != inertia.pos_eigens ) ;
  PivotTolerances use_pivot_tols = S_chol_.pivot_tols();
  if( pivot_tols.warning_tol != PivotTolerances::UNKNOWN )
    use_pivot_tols.warning_tol = pivot_tols.warning_tol;
  if( pivot_tols.singular_tol != PivotTolerances::UNKNOWN )
    use_pivot_tols.singular_tol = pivot_tols.singular_tol;
  if( pivot_tols.wrong_inertia_tol != PivotTolerances::UNKNOWN )
    use_pivot_tols.wrong_inertia_tol = pivot_tols.wrong_inertia_tol;
  throw_exception = (
    *gamma == 0.0
    || !wrong_inertia && *gamma <= use_pivot_tols.warning_tol
    || !wrong_inertia && *gamma <= use_pivot_tols.singular_tol
    || wrong_inertia
    );
  // Create message and throw exceptions
  std::ostringstream onum_msg;
  if(throw_exception) {
    if(wrong_inertia) {
      onum_msg
        << "inertia = ("
        << inertia.neg_eigens << "," << inertia.zero_eigens << "," << inertia.pos_eigens
        << ") != expected_inertia = ("
        << exp_inertia.neg_eigens << "," << exp_inertia.zero_eigens << "," << exp_inertia.pos_eigens << ")"
        << std::endl;
    }
    onum_msg
      << "gamma = min(|diag(D)(k)|)/max(|diag(D)(k)|) = " <<  min_diag << "/" << max_diag
      << " = " << *gamma;
    *omsg
      << "MatrixSymAddDelBunchKaufman::"<<func_name<<"(...): ";
    if( *gamma == 0.0 || (!wrong_inertia && *gamma <= use_pivot_tols.singular_tol) ) {
      *omsg
        << "Singular update!\n" << onum_msg.str() << " <= singular_tol = "
        << use_pivot_tols.singular_tol;
      throw SingularUpdateException( omsg->str(), *gamma );
    }
    else if( wrong_inertia && *gamma <= use_pivot_tols.singular_tol ) {
      *omsg
        << "Singular update!\n" << onum_msg.str() << " <= wrong_inertia_tol = "
        << use_pivot_tols.wrong_inertia_tol;
      throw SingularUpdateException( omsg->str(), *gamma );
    }
    else if( wrong_inertia ) {
      *omsg
        << "Indefinite update!\n" << onum_msg.str() << " >= wrong_inertia_tol = "
        << use_pivot_tols.wrong_inertia_tol;
      throw WrongInertiaUpdateException( omsg->str(), *gamma );
    }
    else if( !wrong_inertia && *gamma <= use_pivot_tols.warning_tol ) {
      *omsg
        << "Warning, near singular update!\nsingular_tol = " << use_pivot_tols.singular_tol
        << " < " << onum_msg.str() << " <= warning_tol = "
        << use_pivot_tols.warning_tol;
      // Don't throw the exception till the very end!
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(true); // Only local programming error?
    }
  }
  // Return
  *comp_inertia = inertia;
  return throw_exception;
}

} // end namespace ConstrainedOptPack
