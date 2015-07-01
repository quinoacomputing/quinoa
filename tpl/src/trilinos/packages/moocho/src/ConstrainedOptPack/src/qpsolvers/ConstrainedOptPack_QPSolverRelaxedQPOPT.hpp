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

#ifdef CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT

#ifndef QP_SOLVER_RELAXED_QPOPT_H
#define QP_SOLVER_RELAXED_QPOPT_H

#include "ConstrainedOptPack_QPSolverRelaxedQPOPTSOL.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief QPSolver subclass that uses QPOPT.
 *
 * ToDo: Finish documentation.
 */  
class QPSolverRelaxedQPOPT : public QPSolverRelaxedQPOPTSOL
{
public:

  /** \brief . */
  typedef QPSolverRelaxedQPOPTSOL inherited;

  /** \brief Set the maximum number of QP iterations as max_itr = max_qp_iter_frac * n.
    */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac );

  /** \brief . */
  QPSolverRelaxedQPOPT(
    value_type max_qp_iter_frac	= 10.0
    );

  /** \brief . */
  ~QPSolverRelaxedQPOPT();

  // /////////////////////////////////
  // Overridden from QPSolverRelaxed

  /** \brief . */
  void release_memory();

protected:

  // /////////////////////////////////////////////////////////////
  // Overridden from QPSolverRelaxedQPOPTSOL

  /** \brief . */
  f_int liwork(f_int N, f_int NCLIN) const;
  /** \brief . */
  f_int lrwork(f_int N, f_int NCLIN) const;
  /** \brief . */
  EInform call_qp_solver(bool warm_start);

private:

  // ////////////////////////////
  // Private types

  /** \brief . */
  enum EQPOPTInform {
    STRONG_LOCAL_MIN      = 0,
    WEAK_LOCAL_MIN        = 1,
    UNBOUNDED             = 2,
    INFEASIBLE            = 3,
    ITMAX_EXCEEDED        = 4,
    MAX_DOF_TOO_SMALL     = 5,
    INVALID_INPUT         = 6,
    PROB_TYPE_NOT_REGOG   = 7
  };

  // ////////////////////////////
  // Private data members

  // extra QPOPT control and input parameters.

  // control

  f_int        ITMAX_;
  f_dbl_prec   BIGBND_;
  f_dbl_prec   FEATOL_;

  // input/output

  f_int           LDA_;
  f_int           LDH_;
  f_dbl_prec*     H_;
  f_int           INFORM_;

  // ////////////////////////////
  // Private member functions

  // not defined and not to be called.
  QPSolverRelaxedQPOPT(const QPSolverRelaxedQPOPT&);
  QPSolverRelaxedQPOPT& operator=(const QPSolverRelaxedQPOPT&);

};	// end class QPSolverRelaxedQPOPT

}	// end namespace ConstrainedOptPack

#endif // QP_SOLVER_RELAXED_QPOPT_H

#endif // CONSTRAINED_OPTIMIZATION_PACK_USE_QPOPT
