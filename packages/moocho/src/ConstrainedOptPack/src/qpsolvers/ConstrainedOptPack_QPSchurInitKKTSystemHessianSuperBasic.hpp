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

#ifndef QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_SUPER_BASIC_H
#define QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_SUPER_BASIC_H

#include "ConstrainedOptPack_QPSolverRelaxedQPSchur.hpp"
#include "ConstrainedOptPack_QPSchurInitKKTSystemHessianFull.hpp"

namespace ConstrainedOptPack {

/** \brief Implementation of initial KKT system for all variables initially fixed
 * and free where #Ko = B_RR#.
 *
 * In this implementation, #G# must support the \Ref{MatrixHessianSuperBasic}
 * interface.
 */
class QPSchurInitKKTSystemHessianSuperBasic
  : public QPSolverRelaxedQPSchur::InitKKTSystem 
{
public:

  // ////////////////////////////////
  // Overridden from InitKKTSystem

  /** \brief Initialize the KKT system where the variables are initiallly 
   * fixed and free and no constraints are in Ko.
   *
   * The Hessian for the QP without the relaxation #G# is represented as
   * a \Ref{MatrixHessianSuperBasic} object and is:
     *
   * #G = Q_R*B_RR*Q_R' + Q_R*op(B_RX)*Q_X' + Q_X*op(B_RX')*Q_R + Q_X*B_XX*Q_X'#
   *
   * If #G# does not support the interface #MatrixHessianSuperBasic# then
   * an exception will be thrown.
   *
   * Given the above parts of #G#, define: #[nd,nd_R] = size(G.Q_R)# and
   * #[nd,nd_X] = size(G.Q_X)#.  Then initial KKT system is defined as:
   *
   * #n_R = nd_R#\\
   * if #i_x_free.size() > 0#
   *    then #i_x_free[(G.Q_R.begin()+l-1)->col_j()-1] = (G.Q_R.begin()+l-1)->row_i(), l = 1...nd_R#\\
   * if #i_x_free.size() == 0# then #i_x_free is implicitly identity#\\
   * #i_x_fixed[(G.Q_X.begin()+l-1)->col_j()-1] = (G.Q_X.begin()+l-1)->row_i(), l = 1...nd_X#\\
   * #i_x_fixed[nd_X] = nd+1#\\
   * #bnd_fixed[l-1] = G.bnd_fixed[l-1], l = 1...nd_X#\\
   * #bnd_fixed[nd_X] = LOWER#\\
   * #j_f_decomp[] = empty#\\
   * #b_X[l-1] = { dL(i) if bnd_fixed[l-1] == LOWER or EQUALITY, dU(i) if bnd_fixed[l-1] == UPPER }#
   * #, l = 1...nd_X (where i = i_x_fixed[l-1])#\\
   * #b_X[nd_X] = etaL#\\
   * #Ko = G.B_RR#\\
   * #fo = - G.Q_R'*g - op(G.B_RX)*b_X(1:nd_X)#\\\
   *
   * Above, it is assumed that if #G.bnd_fixed[l-1] == EQUALITY#, that
   * #dL(G.i_x_fixed[l-1]) == dU(G.i_x_fixed[l-1]# but this may not be
   * inforced by this class.
   *
   * If the \Ref{MatrixHessianSuperBasic} interface is not suppored by #G# then a
   * \Ref{QPSchurInitKKTSystemHessianFull} strategy object is used to
   * try to initialize the KKT system.
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

}; // end class QPSchurInitKKTSystemHessianSuperBasic

} // end namesapce ConstrainedOptPack

#endif // QPSCHUR_INIT_KKT_SYSTEM_HESSIAN_SUPER_BASIC_H
