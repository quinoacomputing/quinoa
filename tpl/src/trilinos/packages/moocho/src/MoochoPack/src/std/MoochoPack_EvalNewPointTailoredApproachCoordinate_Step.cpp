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

#include "MoochoPack_EvalNewPointTailoredApproachCoordinate_Step.hpp"
#include "ConstrainedOptPack_MatrixIdentConcatStd.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_MatrixZero.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

EvalNewPointTailoredApproachCoordinate_Step::EvalNewPointTailoredApproachCoordinate_Step(
  const deriv_tester_ptr_t& 	  deriv_tester
  ,const bounds_tester_ptr_t&	  bounds_tester
  ,EFDDerivTesting              fd_deriv_testing
  )
  :EvalNewPointTailoredApproach_Step(deriv_tester,bounds_tester,fd_deriv_testing)
{}

// protected

void EvalNewPointTailoredApproachCoordinate_Step::uninitialize_Y_Uy(
  MatrixOp         *Y
  ,MatrixOp        *Uy
  )
{
  // Nothing to free
}

void EvalNewPointTailoredApproachCoordinate_Step::calc_py_Y_Uy(
  const NLPDirect       &nlp
  ,const D_ptr_t        &D
  ,VectorMutable        *py
  ,MatrixOp             *Y
  ,MatrixOp             *Uy
  ,EJournalOutputLevel  olevel
  ,std::ostream         &out
  )
{
  namespace rcp = MemMngPack;
  using Teuchos::dyn_cast;

  MatrixIdentConcatStd
    &cY = dyn_cast<MatrixIdentConcatStd>(*Y);
  //
  // Y = [      I     ] space_xD  
  //     [    Zero    ] space_xI
  //        space_xD
  //
  VectorSpace::space_ptr_t
    space_x  = nlp.space_x(),
    space_xD = space_x->sub_space(nlp.var_dep())->clone(),
    space_xI = space_x->sub_space(nlp.var_indep())->clone();
  cY.initialize(
    space_x                                                // space_cols
    ,space_xD                                              // space_rows
    ,MatrixIdentConcatStd::BOTTOM                          // top_or_bottom
    ,1.0                                                   // alpha
    ,Teuchos::rcp(
      new MatrixZero(
        space_xI    // space_cols
        ,space_xD   // space_rows
        ) )                                            // D_ptr
    ,BLAS_Cpp::no_trans                                    // D_trans
    );
  // py is not altered here!
}

void EvalNewPointTailoredApproachCoordinate_Step::recalc_py(
  const MatrixOp          &D
  ,VectorMutable          *py
  ,EJournalOutputLevel    olevel
  ,std::ostream           &out
  )
{
  // py is not altered here!
}

void EvalNewPointTailoredApproachCoordinate_Step::print_calc_py_Y_Uy(
  std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Coordinate decomposition\n"
    << L << "py_k = py_k\n"
    << L << "Y = [ I ; 0 ] <: R^(n x m) [0 represented using MatrixZero]\n"
    << L << "Uy = Gc(var_dep,con_undecomp)\'\n"
    ;
}

}	// end namespace MoochoPack 
