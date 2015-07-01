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

#ifndef RSQP_ALGO_CLIENT_INTERFACE_H
#define RSQP_ALGO_CLIENT_INTERFACE_H

#include "MoochoPack_NLPSolverClientInterface.hpp"

namespace MoochoPack {

/** \brief Interface that smart clients use to set the algorithm configuration
 * object that defines the rSQP algorithm to be used to solve the NLP.
 *
 * ToDo: Finish documentation!
 */
class NLPAlgoClientInterface : public NLPSolverClientInterface {
public:

  /** @name Public Types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<NLPAlgoConfig>	config_ptr_t;

  //@}

  /** @name «std comp» members for config. */
  //@{

  /** \brief . */
  virtual void set_config(const config_ptr_t& config) = 0;
  /** \brief . */
  virtual config_ptr_t& get_config() = 0;
  /** \brief . */
  virtual const config_ptr_t& get_config() const = 0;
  /** \brief . */
  virtual NLPAlgoConfig& config() = 0;
  /** \brief . */
  virtual const NLPAlgoConfig& config() const = 0;

  //@}
  
  /** \brief Causes the algorithm to be configured.
   *
   * Causes the \c config object to configure the algorithm
   * to be ready to solve an NLP or print the algorithm.
   *
   * May be called after the \c nlp, \c track and \c config objects
   * are set.
   *
   * Must be  called before \c print_algorithm() or \c find_min() are called.
    */
  virtual void configure_algorithm(std::ostream* trase_out = 0) = 0;

  /// Print the configured algorithm
  virtual void print_algorithm(std::ostream& out) const = 0;

private:

#ifdef DOXYGEN_COMPILE // Strictly for doxygen diagrams
  /** \brief . */
  NLPAlgoConfig    *config;
#endif

};	// end class NLPAlgoClientInterface

}	// end namespace MoochoPack

#endif	// RSQP_ALGO_CLIENT_INTERFACE_H
