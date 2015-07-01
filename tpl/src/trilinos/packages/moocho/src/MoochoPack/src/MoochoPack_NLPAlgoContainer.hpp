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

#ifndef RSQP_ALGO_CONTAINER_H
#define RSQP_ALGO_CONTAINER_H

#include "MoochoPack_NLPAlgoClientInterface.hpp"
#include "MoochoPack_NLPAlgoInterface.hpp"
#include "MoochoPack_NLPAlgoConfig.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"

namespace MoochoPack {

/** \brief Implementation for NLPAlgo solver.
 *
 * Acts as a container for NLPAlgo.  This class is hidden from clients
 * by not exposing it to them in header files.
 */
class NLPAlgoContainer : public NLPAlgoClientInterface {
public:

  /** @name Constructors / initializers */
  //@{

  /// Members for <<std comp>> of the algorithm object algo.
  STANDARD_COMPOSITION_MEMBERS( NLPAlgoInterface, algo );

  /// Construct a container with no configuration object set.
  NLPAlgoContainer()
  {}

  //@}

  /** @name Overridden from NLPAlgoClientInterface */
  //@{

  /** \brief . */
  void set_config(const config_ptr_t& config);
  /** \brief . */
  config_ptr_t& get_config();
  /** \brief . */
  const config_ptr_t& get_config() const;
  /** \brief . */
  NLPAlgoConfig& config();
  /** \brief . */
  const NLPAlgoConfig& config() const;

  //@}

  /** @name Overridden from NLPSolverClientInterface */
  //@{

  /** \brief . */
  EFindMinReturn find_min();
  /** \brief . */
  void configure_algorithm(std::ostream* trase_out);
  /** \brief . */
  void print_algorithm(std::ostream& out) const;
  /** \brief . */
  void set_algo_timing( bool algo_timing );
  /** \brief . */
  bool algo_timing() const;
  /** \brief . */
  void print_algorithm_times( std::ostream& out ) const;

  //@}

private:

  config_ptr_t			config_;

  // Assert that the object has been set up properly and throw exception if it has not
  void assert_valid_setup() const;

  // Not defined and not to be called
  NLPAlgoContainer(const NLPAlgoContainer&);
  NLPAlgoContainer& operator=(const NLPAlgoContainer&);

};	// end class NLPAlgoContainer

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_CONTAINER_H
