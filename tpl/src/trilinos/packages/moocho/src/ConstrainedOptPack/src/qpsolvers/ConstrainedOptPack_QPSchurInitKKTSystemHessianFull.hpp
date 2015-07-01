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

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H

#include "ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp"

namespace ConstrainedOptPack {

/** \brief Implementation of initial KKT system for all variables initially free
 * and <tt>Ko = G</tt>.
 *
 * In this implementation, \c G must support the \c MatrixSymOpNonsing
 * interface.  Using this initial KKT system essentially make QPSchur use
 * a range space approach (Nocedal & Wright, 1999) for factorizing the KKT
 * system for the current active set.
 */
class QPSchurInitKKTSystemHessianFull
  : public QPSolverRelaxedQPSchur::InitKKTSystem 
{
public:

  /** @name Overridden from InitKKTSystem */
  //@{

  /** \brief Initialize the KKT system where all variables (except the relaxation variable)
   * are initially free and no constraints are in Ko.
   *
   * For this implementation:
   *
   * <tt>n_R = nd</tt>\\
   * <tt>i_x_free = emply (it is identity)</tt>\\
   * <tt>i_x_fixed[0] = nd+1</tt>\\
   * <tt>bnd_fixed[0] = LOWER</tt>\\
   * <tt>j_f_decomp[] = empty</tt>\\
   * <tt>b_X = etaL</tt>\\
   * <tt>Ko = G</tt>\\
   * <tt>fo = -g</tt>\\
   */
  void initialize_kkt_system(
    const Vector    &g
    ,const MatrixOp   &G
    ,value_type           etaL
    ,const Vector   *dL
    ,const Vector   *dU
    ,const MatrixOp   *F
    ,BLAS_Cpp::Transp     trans_F
    ,const Vector   *f
    ,const Vector   *d
    ,const Vector   *nu
    ,size_type            *n_R
    ,i_x_free_t           *i_x_free
    ,i_x_fixed_t          *i_x_fixed
    ,bnd_fixed_t          *bnd_fixed
    ,j_f_decomp_t         *j_f_decomp
    ,DVector               *b_X
    ,Ko_ptr_t             *Ko
    ,DVector               *fo
    ) const;

  //@}

}; // end class QPSchurInitKKTSystemHessianFull

} // end namesapce ConstrainedOptPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_FULL_H
