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

#ifndef EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
#define EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H

#include "MoochoPack_EvalNewPointTailoredApproach_Step.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Implements "orthogonal" decompostion for "Tailored Appraoch".
 *
 * Computes:
 \verbatim
 py = inv(I + D*D') * py
 Y  = [ I; -D' ]
 Uy = ???
 \endverbatim
 */
class EvalNewPointTailoredApproachOrthogonal_Step
  : public EvalNewPointTailoredApproach_Step
{
public:

  /** \brief . */
  EvalNewPointTailoredApproachOrthogonal_Step(
    const deriv_tester_ptr_t                &deriv_tester
    ,const bounds_tester_ptr_t              &bounds_tester
    ,EFDDerivTesting                        fd_deriv_testing = FD_DEFAULT
    );

protected:

  /** @name Overridden from EvalNewPointTailoredApproach_Step */
  //@{

  /** \brief . */
  void uninitialize_Y_Uy(
    MatrixOp         *Y
    ,MatrixOp        *Uy
    );
  /** \brief . */
  void calc_py_Y_Uy(
    const NLPDirect       &nlp
    ,const D_ptr_t        &D
    ,VectorMutable        *py
    ,MatrixOp             *Y
    ,MatrixOp             *Uy
    ,EJournalOutputLevel  olevel
    ,std::ostream         &out
    );
  /** \brief . */
  void recalc_py(
    const MatrixOp           &D
    ,VectorMutable           *py
    ,EJournalOutputLevel     olevel
    ,std::ostream            &out
    );
  /** \brief . */
  void print_calc_py_Y_Uy(
    std::ostream& out, const std::string& leading_str
    ) const;

  //@}

private:

  // ///////////////////////////////
  // Private types

  /** \brief . */
  typedef Teuchos::RCP<MatrixSymOpNonsing>  S_ptr_t;

  // ///////////////////////////////
  // Private data members

  S_ptr_t   S_ptr_;

  // //////////////////////////////
  // Private member functions

  // not defined and not to be called
  EvalNewPointTailoredApproachOrthogonal_Step();

};	// end class EvalNewPointTailoredApproachOrthogonal_Step

}	// end namespace MoochoPack 

#endif	// EVAL_NEW_POINT_TAILORED_APPROACH_ORTHOGONAL_STEP_H
