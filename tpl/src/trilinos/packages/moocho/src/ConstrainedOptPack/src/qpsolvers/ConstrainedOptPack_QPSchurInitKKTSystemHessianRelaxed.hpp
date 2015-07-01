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

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_RELAXED_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_RELAXED_H

#include "ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp"
#include "ConstrainedOptPack_QPSchurInitKKTSystemHessianFull.hpp"

namespace ConstrainedOptPack {

/** \brief Implementation of initial KKT system where all original variables
 * are free and all the relaxation variables are fixed.
 *
 * In this implementation, #G# should support the \Ref{MatrixSymHessianRelaxNonSing}
 * interface.   Otherwise, it will try the \Ref{MatrixSymWithOpFactorized} interface
 * using the base implementation of \Ref{QPSchurInitKKTSystemHessianFull}.
 */
class QPSchurInitKKTSystemHessianRelaxed
  : public QPSchurInitKKTSystemHessianFull
{
public:

  // ////////////////////////////////
  // Overridden from InitKKTSystem

  /** \brief Initialize the KKT system where the original variables are initiallly 
   * free and all the relaxation variables are fixed and their are no
   * constraints in Ko.
   *
   * The Hessian for the QP without the relaxation #G# is represented as
   * a \Ref{MatrixSymHessianRelaxNonSing} object and is:
     \begin{verbatim}
   G = [ G_orig     ]
       [          M ]
   \end{verbatim}
   * If #G# does not support the interface #MatrixSymHessianRelaxNonSing# then
   * the function #QPSchurInitKKTSystemHessianFull::initialize_kkt_system(...)#
   * will be called.
   *
   * Given the above parts of #G#, define: #[no,no] = size(G.G)# and
   * #[nr,nr] = size(G.M)#.  Then initial KKT system is defined as:
   *
   * #n_R = no#\\
   * #i_x_free.size() == 0# and #i_x_free is implicitly identity#\\
   * #i_x_fixed[l-1] = no + l, l = 1...nr#\\
   * #i_x_fixed[nr] = no+nr+1#\\
   * #bnd_fixed[l-1] = LOWER, l = 1...nr#\\
   * #bnd_fixed[nr] = LOWER#\\
   * #j_f_decomp[] = empty#\\
   * #b_X[l-1] = dL(no+l), l = 1...nr#\\
   * #b_X[nr] = etaL#\\
   * #Ko = G.G#\\
   * #fo = - g(1:no)#\\\
   */
  void initialize_kkt_system(
    const DVectorSlice&    g
    ,const MatrixOp&  G
    ,value_type           etaL
    ,const SpVectorSlice& dL
    ,const SpVectorSlice& dU
    ,const MatrixOp*  F
    ,BLAS_Cpp::Transp     trans_F
    ,const DVectorSlice*   f
    ,const DVectorSlice&   d
    ,const SpVectorSlice& nu
    ,size_type*           n_R
    ,i_x_free_t*          i_x_free
    ,i_x_fixed_t*         i_x_fixed
    ,bnd_fixed_t*         bnd_fixed
    ,j_f_decomp_t*        j_f_decomp
    ,DVector*              b_X
    ,Ko_ptr_t*            Ko
    ,DVector*              fo
    ) const;

private:
  QPSchurInitKKTSystemHessianFull  init_kkt_full_;

}; // end class QPSchurInitKKTSystemHessianRelaxed

} // end namesapce ConstrainedOptPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_RELAXED_H
