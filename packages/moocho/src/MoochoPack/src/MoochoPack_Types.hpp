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

#ifndef REDUCED_SPACE_SQP_PACK_TYPES_H
#define REDUCED_SPACE_SQP_PACK_TYPES_H

#include "ConstrainedOptPack_Types.hpp"
#include "IterationPack_Types.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerbosityLevel.hpp"

namespace MoochoPack {

using Teuchos::RCP;

// using types from ConstrainedOptPack
#include "ConstrainedOptPack_PublicTypes.ud"

// using types from IterationPack
#include "IterationPack_PublicTypes.ud"

/** \brief enum for journal output. */
enum EJournalOutputLevel {
  PRINT_NOTHING = 0
  ,PRINT_BASIC_ALGORITHM_INFO = 1
  ,PRINT_ALGORITHM_STEPS = 2
  ,PRINT_ACTIVE_SET = 3
  ,PRINT_VECTORS = 4
  ,PRINT_ITERATION_QUANTITIES = 5
};

/** \brief Conver to Teuchos::EVerbosityLevel. */
inline Teuchos::EVerbosityLevel convertToVerbLevel( const EJournalOutputLevel output_level )
{
  switch(output_level) {
    case PRINT_NOTHING:
    case PRINT_BASIC_ALGORITHM_INFO:
      return Teuchos::VERB_NONE;
    case PRINT_ACTIVE_SET:
    case PRINT_ALGORITHM_STEPS:
      return Teuchos::VERB_LOW;
    case PRINT_VECTORS:
      return Teuchos::VERB_HIGH;
    case PRINT_ITERATION_QUANTITIES:
      return Teuchos::VERB_EXTREME;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should never get here!
  }
  return Teuchos::VERB_NONE; // Should never be called!
}

// public interface classes

class NLPAlgoState;
class NLPSolverClientInterface;
class NLPAlgoClientInterface;
class NLPAlgoConfig;

//

class NLPAlgo;

}	// end namespace MoochoPack 

#endif // REDUCED_SPACE_SQP_PACK_TYPES_H
