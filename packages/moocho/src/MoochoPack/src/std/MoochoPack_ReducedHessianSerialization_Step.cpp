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

#include <fstream>

#include "MoochoPack_ReducedHessianSerialization_Step.hpp"
#include "MoochoPack_moocho_algo_conversion.hpp"
#include "MoochoPack_Exceptions.hpp"
#include "IterationPack_print_algorithm_step.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixSymInitDiag.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "SerializationPack_Serializable.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace MoochoPack {

ReducedHessianSerialization_Step::ReducedHessianSerialization_Step(
  const std::string    &reduced_hessian_input_file_name
  ,const std::string   &reduced_hessian_output_file_name
  )
  :reduced_hessian_input_file_name_(reduced_hessian_input_file_name)
  ,reduced_hessian_output_file_name_(reduced_hessian_output_file_name)
{}

// Overridden from AlgorithmStep

bool ReducedHessianSerialization_Step::do_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  using Teuchos::dyn_cast;
  using SerializationPack::Serializable;

  NLPAlgo       &algo = rsqp_algo(_algo);
  NLPAlgoState  &s    = algo.rsqp_state();
  
  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  EJournalOutputLevel ns_olevel = algo.algo_cntr().null_space_journal_output_level();
  std::ostream& out = algo.track().journal_out();

  // print step header.
  if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
    using IterationPack::print_algorithm_step;
    print_algorithm_step( _algo, step_poss, type, assoc_step_poss, out );
  }

  IterQuantityAccess<MatrixSymOp>  &rHL_iq = s.rHL();

  if( !rHL_iq.updated_k(0) && reduced_hessian_input_file_name().length() ) {
    int k_last_offset = rHL_iq.last_updated();
    if( k_last_offset == IterQuantity::NONE_UPDATED ) {
      if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out
          << "\nNo previous matrix rHL was found!\n"
          << "\nReading in the matrix rHL_k from the file \""<<reduced_hessian_input_file_name()<<"\" ...\n";
      }
      MatrixSymOp	&rHL_k = rHL_iq.set_k(0);
      Serializable &rHL_serializable = dyn_cast<Serializable>(rHL_k);
      std::ifstream reduced_hessian_input_file(reduced_hessian_input_file_name().c_str());
      TEUCHOS_TEST_FOR_EXCEPTION(
        !reduced_hessian_input_file, std::logic_error
        ,"ReducedHessianSerialization_Step::do_step(...): Error, the file \""<<reduced_hessian_input_file_name()<<"\""
        " could not be opened or contains no input!"
        );
      rHL_serializable.unserialize(reduced_hessian_input_file);
      if( (int)ns_olevel >= (int)PRINT_ALGORITHM_STEPS ) {
        out << "\nrHL_k.rows() = " << rHL_k.rows() << std::endl;
        out << "\nrHL_k.cols() = " << rHL_k.cols() << std::endl;
        if(algo.algo_cntr().calc_matrix_norms())
          out << "\n||rHL_k||inf    = " << rHL_k.calc_norm(MatrixOp::MAT_NORM_INF).value << std::endl;
        if(algo.algo_cntr().calc_conditioning()) {
          const MatrixSymOpNonsing
            *rHL_ns_k = dynamic_cast<const MatrixSymOpNonsing*>(&rHL_k);
          if(rHL_ns_k)
            out << "\ncond_inf(rHL_k) = " << rHL_ns_k->calc_cond_num(MatrixOp::MAT_NORM_INF).value << std::endl;
        }
      }
      if( (int)ns_olevel >= (int)PRINT_ITERATION_QUANTITIES )
        out << "\nrHL_k = \n" << rHL_k;
      // Validate the space
      const MatrixOp &Z_k = s.Z().get_k(0);
      const VectorSpace &null_space = Z_k.space_rows();
      TEUCHOS_TEST_FOR_EXCEPTION(
        !null_space.is_compatible(rHL_k.space_cols()) || !null_space.is_compatible(rHL_k.space_rows())
        ,std::runtime_error
        ,"ReducedHessianSerialization_Step::do_step(...): Error, the read-in reduced Hessian of dimension "
        << rHL_k.rows() << " x " << rHL_k.cols() << " is not compatible with the null space of dimension "
        "Z_k.cols() = " << Z_k.cols() << "!"
        );
    }
  }
  
  return true;
  
}

void ReducedHessianSerialization_Step::finalize_step(
  Algorithm& _algo, poss_type step_poss, IterationPack::EDoStepType type
  ,poss_type assoc_step_poss
  )
{
  using Teuchos::dyn_cast;
  using SerializationPack::Serializable;

  const NLPAlgo       &algo = rsqp_algo(_algo);
  const NLPAlgoState  &s    = algo.rsqp_state();

  EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
  std::ostream& out = algo.track().journal_out();

  const IterQuantityAccess<MatrixSymOp>  &rHL_iq = s.rHL();
  
  int k_last_offset = rHL_iq.last_updated();

  if( k_last_offset != IterQuantity::NONE_UPDATED && reduced_hessian_output_file_name().length() ) {
    if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
      out
        << "\nSerializing the matrix rHL_k("<<k_last_offset<<") to the file "
        << "\""<<reduced_hessian_output_file_name()<<"\" ...\n";
    }
    const MatrixSymOp	&rHL_k = rHL_iq.get_k(k_last_offset);
    const Serializable &rHL_serializable = dyn_cast<const Serializable>(rHL_k);
    std::ofstream reduced_hessian_output_file(reduced_hessian_output_file_name().c_str());
    TEUCHOS_TEST_FOR_EXCEPTION(
      !reduced_hessian_output_file, std::logic_error
      ,"ReducedHessianSerialization_Step::finalize_step(...): Error, the file \""<<reduced_hessian_output_file_name()<<"\""
      " could not be opened!"
      );
    rHL_serializable.serialize(reduced_hessian_output_file);
  }
}

void ReducedHessianSerialization_Step::print_step(
  const Algorithm& algo, poss_type step_poss, IterationPack::EDoStepType type, poss_type assoc_step_poss
  ,std::ostream& out, const std::string& L
  ) const
{
  out
    << L << "*** Read in the reduced Hessian of the Lagrangian rHL from a file.\n"
    << L << "if (rHL_k is not updated and reduced_hessian_input_file_name != \"\") then\n"
    << L << "  k_last_offset = last iteration rHL was updated for\n"
    << L << "  if k_last_offset==NONE_UPDATED then\n"
    << L << "    rHL_serializable = dyn_cast<Serializable>(rHL_k) *** Throws exception if fails!\n"
    << L << "    Unserialize into rHL_serializable from the file \""<<reduced_hessian_input_file_name()<<"\"\n"
    << L << "  else\n"
    << L << "    *** There is some reduced Hessian that exists for some past iteration so\n"
    << L << "    *** we will let some other step object initialize itQ\n"
    << L << "  end\n"
    << L << "end\n"
    << L << "*** Note: On finalization, this step object will serialize rHL_k to the file:\n"
    << L << "***   \""<<reduced_hessian_output_file_name()<<"\""
    ;
}

}	// end namespace MoochoPack 
