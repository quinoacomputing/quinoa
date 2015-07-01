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

#ifndef MATRIX_SYM_POS_DEF_BUNCH_KAUFMAN_H
#define MATRIX_SYM_POS_DEF_BUNCH_KAUFMAN_H

#include <vector>

#include "ConstrainedOptPack_MatrixSymAddDelUpdateableWithOpNonsingular.hpp"
#include "AbstractLinAlgPack_MatrixSymAddDelUpdateable.hpp"
#include "AbstractLinAlgPack_MatrixSymPosDefCholFactor.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsingSerial.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"

namespace ConstrainedOptPack {

/** \brief This class maintains the factorization of symmetric indefinite matrix
 * using a Bunch & Kaufman factorization.
 *
 * When the matix in question is positive definite or negative definite
 * then a cholesky factorization will be used for greater efficiency.
 * As much as possible the class trys to keep runtime costs down to
 * a minimum while still maintaining a stable factorization.
 *
 * ToDo: This matrix class could also handle the case of a single negative
 * or positive eigen value in an efficient manner as well.
 */
class MatrixSymAddDelBunchKaufman
  :public virtual MatrixSymOpNonsingSerial
  ,public virtual MatrixSymAddDelUpdateable
  ,public virtual MatrixSymAddDelUpdateableWithOpNonsingular
{
public:

  /// Initializes with 0x0 and pivot_tols == (0.0,0.0,0.0).
  MatrixSymAddDelBunchKaufman();

  /// Pivot tolerance used durring the cholesky factorization (it may be zero).
  void pivot_tols( PivotTolerances pivot_tols );
  /** \brief . */
  PivotTolerances	pivot_tols() const;

  /** @name Overridden from MatrixSymAddDelUpdateableWithOpNonsingular */
  //@{

  /** \brief . */
  const MatrixSymOpNonsing& op_interface() const;
  /** \brief . */
  MatrixSymAddDelUpdateable& update_interface();
  /** \brief . */
  const MatrixSymAddDelUpdateable& update_interface() const;

  //@}

  /** @name Overridden from MatrixSymAddDelUpdateable */
  //@{

  /** \brief . */
  void initialize(
    value_type         alpha
    ,size_type         max_size
    );
  /** \brief . */
  void initialize(
    const DMatrixSliceSym      &A
    ,size_type         max_size
    ,bool              force_factorization
    ,Inertia           inertia
    ,PivotTolerances   pivot_tols
    );
  /** \brief . */
  size_type max_size() const;
  /** \brief . */
  Inertia inertia() const;
  /** \brief . */
  void set_uninitialized();
  /** \brief . */
  void augment_update(
    const DVectorSlice  *t
    ,value_type        alpha
    ,bool              force_refactorization
    ,EEigenValType     add_eigen_val
    ,PivotTolerances   pivot_tols
    );
  /** \brief . */
  void delete_update(
    size_type          jd
    ,bool              force_refactorization
    ,EEigenValType     drop_eigen_val
    ,PivotTolerances   pivot_tols
    );

  //@}

  /** @name Overridden from MatrixSymOpNonsingSerial */
  //@{

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(
    DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const DVectorSlice& vs_rhs2, value_type beta
    ) const;
  /** \brief . */
  void V_InvMtV(
    DVectorSlice* vs_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const DVectorSlice& vs_rhs2
    )const;

private:

  // /////////////////////////////////////////////////////
  // Private types

  typedef std::vector<FortranTypes::f_int> IPIV_t;

  // /////////////////////////////////////////////////////
  // Private data members

  size_type S_size_;      // The size of the current symmetric matrix.  Size == 0 if flag for uninitialized.
  bool      S_indef_;     // True if the current matrix is indefinite.
  bool      fact_updated_;// True if the factorization for the current matrix is updated.  Only meaningful
                          // if S_indef_==true.
  bool      fact_in1_;    // If true then the current factorization is in S_store1_
                          // otherwise it is in S_store2_.  Only meaningful if S_indef_==true and
                          // fact_updated_==true
  MatrixSymAddDelUpdateable::Inertia
            inertia_;     // Inertial for the indefinite L*D*L' factorization.  If fact_updated_ == false
                          // then this will be UNKNOWN.  IF S_indef_==false then this is meaningless.
  DMatrix S_store1_;    // Storage for the factored matrix in the
                          // upper triangle as well as the original matrix
                          // in the lower triangle.  This uses the same storage scheme as
                          // in MatrixSymPosDefCholFactor.
  DMatrix S_store2_;    // Storage for the factorization also.  We need
                          // two storage locations for the L*D*L factorization
                          // in case an update is singular.  This will not
                          // be initialized for a p.d. or n.d. matrix.
  IPIV_t    IPIV_;        // Stores permutations computed by LAPACK
  mutable DVector
              WORK_;        // workspace
  MatrixSymPosDefCholFactor
          S_chol_;      // Used to factorize the matrix
                          // when it is p.d. or n.d.

  // /////////////////////////////////////////////////////
  // Private member funcitons.

  /** \brief Get view of DU.
   */
  DMatrixSliceTriEle DU(size_type S_size, bool fact_in1);
  /** \brief . */
  const DMatrixSliceTriEle DU(size_type S_size, bool fact_in1) const;
  /** \brief Get view of lower part of S.
   */
  DMatrixSliceSym S(size_type S_size);
  /** \brief . */
  const DMatrixSliceSym S(size_type S_size) const;
  /** \brief . */
  void assert_initialized() const;
  /** \brief . */
  void resize_DU_store( bool in_store1 );
  /** \brief Copy the original matrix into the new storage location and factorize it.
   *
   * Will throw DenseLinAlgLAPack::FactorizationException if singular.
   */
  void copy_and_factor_matrix( size_type S_size, bool fact_in1 );
  /** \brief Factor the current set matrix in-place (do not copy the original). 
   *
   * Will throw DenseLinAlgLAPack::FactorizationException if singular.
   */
  void factor_matrix( size_type S_size, bool fact_in1 );
  /** \brief Compute the new inertia and validate that it is what the client says it was.
   *
   * Will throw exceptions if the matrix is singular or has the wrong inertia.  If
   * the matrix is near singular then true will be returned, the update should
   * succeed but a warning exception should be thrown
   */
  bool compute_assert_inertia(
    size_type S_size, bool fact_in1
    ,const Inertia& expected_inertia, const char func_name[]
    ,PivotTolerances pivot_tols, Inertia* comp_inertia, std::ostringstream* err_msg, value_type* gamma );

  /// Not defined and not to be called.
  MatrixSymAddDelBunchKaufman( const MatrixSymAddDelBunchKaufman& );
  MatrixSymAddDelBunchKaufman& operator=( const MatrixSymAddDelBunchKaufman& );

};	// end class MatrixSymAddDelBunchKaufman

// ///////////////////////////
// Inline member functions

inline
DMatrixSliceTriEle MatrixSymAddDelBunchKaufman::DU(size_type S_size, bool fact_in1)
{
  resize_DU_store(fact_in1);
  return DenseLinAlgPack::nonconst_tri_ele(
    ( fact_in1 ? S_store1_ : S_store2_ )(1,S_size,2,S_size+1)
    ,BLAS_Cpp::upper );
}

inline
const DMatrixSliceTriEle MatrixSymAddDelBunchKaufman::DU(size_type S_size, bool fact_in1) const
{
  return DenseLinAlgPack::tri_ele(
    ( fact_in1 ? S_store1_ : S_store2_ )(1,S_size,2,S_size+1)
    ,BLAS_Cpp::upper);
}

inline
DMatrixSliceSym MatrixSymAddDelBunchKaufman::S(size_type S_size)
{
  return DenseLinAlgPack::nonconst_sym(
    S_store1_(2,S_size+1,1,S_size)
    , BLAS_Cpp::lower );
}

inline
const DMatrixSliceSym MatrixSymAddDelBunchKaufman::S(size_type S_size) const
{
  return DenseLinAlgPack::sym(
    S_store1_(2,S_size+1,1,S_size)
    , BLAS_Cpp::lower );
}

}	// namespace ConstrainedOptPack 

#endif	// MATRIX_SYM_POS_DEF_BUNCH_KAUFMAN_H
