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

#ifndef MATRIX_GEN_BANDED_H
#define MATRIX_GEN_BANDED_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Miref_count_ptr.h"
#include "MiReleaseResource.h"

namespace ConstrainedOptPack {
/** \brief Matrix subclass for general (possibly singular) banded matrices.
 * 
 * The banded matrix is stored by column in a simple flat rectangular matrix.
 * For example, for #m = 10, n = 8, kl = 3, ku = 2# the matrix #M# is stored in the
 * following format #MB# (same as for the BLAS routine xGBMV(...)):
 \begin{verbatim}

         M                                   MB
 [ x x x               ]
 [ x x x x             ]         [ o o x x x x x x x x ] \ ku = 2
 [ x x x x x           ]         [ o x x x x x x x x o ] /
 [ x x x x x x         ]    =>   [ x x x x x x x x o o ]
 [   x x x x x x       ]         [ x x x x x x x o o o ] \
 [     x x x x x x     ]         [ x x x x x x o o o o ] | kl = 3
 [       x x x x x x   ]         [ x x x x x o o o o o ] /
 [         x x x x x x ]           1 2 3 4 5 6 7 8 9 0
   1 2 3 4 5 6 7 8 9 0

 \end{verbatim}
 */
class MatrixGenBanded	: public MatrixOp
{
public:
  
  /** \brief . */
  typedef Teuchos::RCP<
    MemMngPack::ReleaseResource>  release_resource_ptr_t;

  // //////////////
    // Constructors

  /** \brief Construct and Initialize.
   *
   * This constructor just calls #this->initialize(...)#.
   */
  MatrixGenBanded(
    size_type                         m                       = 0
    ,size_type                        n                       = 0
    ,size_type                        kl                      = 0
    ,size_type                        ku                      = 0
    ,DMatrixSlice                   *MB                     = NULL
    ,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
    );

  // ///////////////////////////
    // Access representation

  /** \brief Initialize
   *
   * If called with all of the default arguments then #this# will become uninitialized.
   *
   * ToDo: Finish pre and post conditions!
   *
   * @param  m        [in] Determines the size of the banded matrix (m x n).  If
   *                  If #m == 0# then all of the following arguments should be left at
   * @param  n        [in] Determines the size of the banded matrix (m x n).
   * @param  kl       [in] Determines the lower band width of the matrix as defined by xGBMV(...).
   * @param  ku       [in] Determines the band width of the matrix as defined by xGBMV(...).
   * @param  MB       [in/state] If #MB != NULL# then this matrix (size (kl+ku+1) x n) is used to store
   *                  the original banded matrix #M# in the format of xGBMV(...).  This matrix must
   *                  be initialized on input.
   * @param  MB_release_resource_ptr
   *                  [in] Only significant if #MB != NULL#.  Points to a resource to
   *                  be released when #MB# is no longer needed.
   */
  void initialize(
    size_type                         m                       = 0
    ,size_type                        n                       = 0
    ,size_type                        kl                      = 0
    ,size_type                        ku                      = 0
    ,DMatrixSlice                   *MB                     = NULL
    ,const release_resource_ptr_t&    MB_release_resource_ptr = NULL
    );

  /** \brief . */
  size_type kl() const;
  /** \brief . */
  size_type ku() const;
  /** \brief Get view of MB.
   */
  DMatrixSlice& MB();
  /** \brief . */
  const DMatrixSlice& MB() const;

  // /////////////////////////////
  // Overridden from MatrixOp

  /** \brief . */
  size_type rows() const;
  /** \brief . */
  size_type cols() const;
  /** \brief . */
  size_type nz() const;
  /** \brief . */
  std::ostream& output(std::ostream& out) const;
  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const DVectorSlice& vs_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(DVectorSlice* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    , const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const DVectorSlice& vs_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(DVectorSlice* vs_lhs, value_type alpha
    , const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    , BLAS_Cpp::Transp M_rhs2_trans
    , const SpVectorSlice& sv_rhs3, value_type beta) const;

private:
  
  // /////////////////////////////
  // Private data members

  size_type                       m_;
  size_type                       n_;
  size_type                       kl_;
  size_type                       ku_;
  DMatrixSlice                  MB_;
  release_resource_ptr_t          MB_release_resource_ptr_;

  // /////////////////////////////
  // Private member functions

  void assert_initialized() const;

}; // end class MatrixGenBanded

// ///////////////////////////////////////////////////////
// Inline members for MatrixGenBanded

inline
size_type MatrixGenBanded::kl() const
{
  return kl_;
}

inline
size_type MatrixGenBanded::ku() const
{
  return ku_;
}

inline
DMatrixSlice& MatrixGenBanded::MB()
{
  return MB_;
}

inline
const DMatrixSlice& MatrixGenBanded::MB() const
{
  return MB_;
}

} // end namespace ConstrainedOptPack

#endif // MATRIX_GEN_BANDED_H
