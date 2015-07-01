#if 0

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

#include <sstream>
#include <typeinfo>
#include <iomanip>

#include "MoochoPack_ReducedHessianExactStd_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "AbstractLinAlgPack_MatrixSymDenseInitialize.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "ConstrainedOptPack/src/VectorWithNorms.h"
#include "NLPInterfacePack_NLPSecondOrder.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DMatrixAsTriSym.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "Midynamic_cast_verbose.h"

namespace MoochoPack {

bool ReducedHessianExactStd_Step::do_step(
    Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  , poss_type assoc_step_poss)
{
  using Teuchos::dyn_cast;
  using DenseLinAlgPack::nonconst_sym;
  using AbstractLinAlgPack::Mp_StMtMtM;
  typedef AbstractLinAlgPack::MatrixSymDenseInitialize	MatrixSymDenseInitialize;
  typedef AbstractLinAlgPack::MatrixSymOp			MatrixSymOp;
  using ConstrainedOptPack::NLPSecondOrder;

  NLPAlgo	&algo	= rsqp_algo(_algo);
  NLPAlgoState	&s		= algo.rsqp_state();
  NLPSecondOrder
#ifdef _WINDOWS
        &nlp	= dynamic_cast<NLPSecondOrder&>(algo.nlp());
#else
        &nlp	= dyn_cast<NLPSecondOrder>(algo.nlp());
#endif
  MatrixSymOp
    *HL_sym_op = dynamic_cast<MatrixSymOp*>(&s.HL().get_k(0));

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
  }

  // problem size
  size_type	n		= nlp.n(),
        r		= nlp.r(),
        nind	= n - r;

  // Compute HL first (You may want to move this into its own step later?)

  if( !s.lambda().updated_k(-1) ) {
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "Initializing lambda_km1 = nlp.get_lambda_init ... \n";
    }
    nlp.get_init_lagrange_mult( &s.lambda().set_k(-1).v(), NULL );
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
      out << "||lambda_km1||inf = " << s.lambda().get_k(-1).norm_inf() << std::endl;
    }
    if( static_cast<int>(olevel) >= static_cast<int>(PRINT_VECTORS) ) {
      out << "lambda_km1 = \n" << s.lambda().get_k(-1)();
    }
  }

  nlp.set_HL(	HL_sym_op );
  nlp.calc_HL( s.x().get_k(0)(), s.lambda().get_k(-1)(), false );

  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ITERATION_QUANTITIES) ) {
    s.HL().get_k(0).output( out << "\nHL_k = \n" );
  }

  // If rHL has already been updated for this iteration then just leave it.
  if( !s.rHL().updated_k(0) ) {

    if( !HL_sym_op ) {
      std::ostringstream omsg;
      omsg
        << "ReducedHessianExactStd_Step::do_step(...) : Error, "
        << "The matrix HL with the concrete type "
        << typeName(s.HL().get_k(0)) << " does not support the "
        << "MatrixSymOp iterface";
      throw std::logic_error( omsg.str() );
    }		

    MatrixSymDenseInitialize
      *rHL_sym_init = dynamic_cast<MatrixSymDenseInitialize*>(&s.rHL().set_k(0));
    if( !rHL_sym_init ) {
      std::ostringstream omsg;
      omsg
        << "ReducedHessianExactStd_Step::do_step(...) : Error, "
        << "The matrix rHL with the concrete type "
        << typeName(s.rHL().get_k(0)) << " does not support the "
        << "MatrixSymDenseInitialize iterface";
      throw std::logic_error( omsg.str() );
    }		

    // Compute the dense reduced Hessian
    DMatrix rHL_sym_store(nind,nind);
    DMatrixSliceSym rHL_sym(rHL_sym_store(),BLAS_Cpp::lower);
    Mp_StMtMtM( &rHL_sym, 1.0, MatrixSymOp::DUMMY_ARG, *HL_sym_op
          , s.Z().get_k(0), BLAS_Cpp::no_trans, 0.0 );

    if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
      out << "\nLower triangular partion of dense reduced Hessian (ignore nonzeros above diagonal):\nrHL_dense = \n" << rHL_sym_store(); 
    }
  
    // Set the reduced Hessain
    rHL_sym_init->initialize( rHL_sym );

    if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
      s.rHL().get_k(0).output( out << "\nrHL_k = \n" );
    }
  }

  return true;
}

void ReducedHessianExactStd_Step::print_step(
    const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type
  , poss_type assoc_step_poss, std::ostream& out, const std::string& L ) const
{
  out
    << L << "*** Calculate the exact reduced Hessian of the Lagrangian\n"
    << L << "if lambda_km1 is not updated then\n"
    << L << "    lambda_km1 = nlp.get_lambda_init\n"
    << L << "end\n"
    << L << "HL_k = HL(x_k,lambda_km1) <: R^(n+m) -> R^(n x n)\n"
    << L << "if rHL_k is not updated then\n"
    << L << "    rHL_dense = Z_k' * HL_k * Z_k  (MatrixSymOp interface for HL_k)\n"
    << L << "    rHL_k = rHL_dense (MatrixSymDenseInitialize interface for rHL_k)\n"
    << L << "end\n";
}

}	// end namespace MoochoPack 

#endif // 0
