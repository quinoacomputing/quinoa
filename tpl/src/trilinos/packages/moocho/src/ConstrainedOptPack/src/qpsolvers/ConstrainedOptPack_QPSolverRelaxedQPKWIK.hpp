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

#ifndef QP_SOLVER_RELAXED_QPKWIK_H
#define QP_SOLVER_RELAXED_QPKWIK_H

#include <vector>

#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_F77_wrappers.h"

namespace ConstrainedOptPack {

/** \brief Solves Quadratic Programming (QP) problem using the primal-dual active-set
 * solver QPKWIK.
 *
 * In this implementation it is required that G support the \Ref{MatrixExtractInvCholFactor}
 * interface and is therefore quite restrictive on the likes of QPs it can solve.
 */
class QPSolverRelaxedQPKWIK : public QPSolverRelaxed
{
public:

  /** @name Initialization */
  //@{

  /// Set the maximum number of QP iterations as max_itr = max_qp_iter_frac * n.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac );

  /// Set the value of an infinite bound.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, infinite_bound );

  /** \brief . */
  QPSolverRelaxedQPKWIK(
      value_type        max_qp_iter_frac	= 10.0
      ,value_type       infinite_bound      = 1e+20
    );

  /** \brief . */
  ~QPSolverRelaxedQPKWIK();

  //@}

  /** @name Overridden from QPSolverRelaxed */
  //@{

  /** \brief . */
  QPSolverStats get_qp_stats() const;
  /** \brief . */
  void release_memory();

  //@}

protected:

  /** @name Overridden from QPSolverRelaxed */
  //@{

  /** \brief . */
  QPSolverStats::ESolutionType imp_solve_qp(
    std::ostream* out, EOutputLevel olevel, ERunTests test_what
    ,const Vector& g, const MatrixSymOp& G
    ,value_type etaL
    ,const Vector* dL, const Vector* dU
    ,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
    ,const Vector* eL, const Vector* eU
    ,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
    ,value_type* obj_d
    ,value_type* eta, VectorMutable* d
    ,VectorMutable* nu
    ,VectorMutable* mu, VectorMutable* Ed
    ,VectorMutable* lambda, VectorMutable* Fd
    );

  //@}

private:

  // //////////////////////////////////////////////////////////////
  // Private types

  /** \brief . */
  typedef FortranTypes::f_int f_int;
  /** \brief . */
  typedef std::vector<f_int>  IBND_t;
  /** \brief . */
  typedef std::vector<f_int>  IACTSTORE_t;
  /** \brief . */
  typedef std::vector<f_int>  IACT_t;
  /** \brief . */
  typedef std::vector<f_int>  ISTATE_t;

  // //////////////////////////////////////////////////////////////
  // Private Data Members.

  QPSolverStats   qp_stats_;

  // Inverse mapping for IBND_INV(j) == k <-> IBND(k) == j
  IBND_t          IBND_INV_;

  // Parameters to QPKWIK

  /** \brief . */
  f_int      N_;
  /** \brief . */
  f_int      M1_;
  /** \brief . */
  f_int      M2_;
  /** \brief . */
  f_int      M3_;
  /** \brief . */
  DVector          GRAD_;
  /** \brief . */
  DMatrix       UINV_AUG_;
  /** \brief . */
  f_int      LDUINV_AUG_;
  /** \brief . */
  IBND_t          IBND_;
  /** \brief . */
  DVector          BL_;
  /** \brief . */
  DVector          BU_;
  /** \brief . */
  DMatrix       A_;
  /** \brief . */
  f_int		LDA_;
  /** \brief . */
  DVector          YPY_;
  /** \brief . */
  f_int      IYPY_;
  /** \brief . */
  f_int      WARM_;
  /** \brief . */
  value_type      NUMPARAM_[3];
  /** \brief . */
  f_int      MAX_ITER_;

  // Input / Output

  /** \brief . */
  DVector          X_;
  /** \brief . */
  f_int      NACTSTORE_;
  /** \brief . */
  IACTSTORE_t     IACTSTORE_;
  /** \brief . */
  f_int      INF_;
  
  // Output

  /** \brief . */
  f_int      NACT_;
  /** \brief . */
  IACT_t          IACT_;
  /** \brief . */
  DVector          UR_;
  /** \brief . */
  value_type      EXTRA_;
  /** \brief . */
  f_int      ITER_;
  /** \brief . */
  f_int      NUM_ADDS_;
  /** \brief . */
  f_int      NUM_DROPS_;
  
  // Internal state

  /** \brief . */
  ISTATE_t        ISTATE_;

  // Workspace

  /** \brief . */
  f_int      LRW_;
  /** \brief . */
  DVector          RW_;

}; // end class QPSolverRelaxedQPKWIK

} // end namespace ConstrainedOptimizationPackTypes

#endif // QP_SOLVER_RELAXED_QPKWIK_H
