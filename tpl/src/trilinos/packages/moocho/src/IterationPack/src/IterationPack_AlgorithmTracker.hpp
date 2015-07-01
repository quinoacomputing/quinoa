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

#ifndef ALGORITHM_TRACK_H
#define ALGORITHM_TRACK_H

#include <iosfwd>

#include "IterationPack_Types.hpp"
#include "Teuchos_RCP.hpp"

namespace IterationPack {

/** \brief Used to ouput iteration results and other information.
  *
  * This interface can be implemented by outside clients of an iterative
  * algorithm to monitor or "track" the progress of the algorithm.
  *
  * ToDo: Write more documentation!
  */
class AlgorithmTracker {
public:

  /** \brief . */
  virtual ~AlgorithmTracker() {}

  /** @name Public types */
  //@{

  /** \brief . */
  typedef Teuchos::RCP<std::ostream>    ostream_ptr_t;

  //@}

  /** @name Constructors */
  //@{

  /** \brief Construct with an output stream for journal_out.
   *
   * Preconditions:<ul>
   * <li> <tt>journal_out.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>).
   * </ul>
   */
  AlgorithmTracker(const ostream_ptr_t& journal_out);
  
  //@}
  
  /** @name Algorithm iteration state notification */
  //@{
  
  /** \brief Reinitialize the track object right before it is used.
   *
   * The default implementation does nothing.
   */
  virtual void initialize();

  /** \brief Output information about an iteration just completed.
    *
    * The default just does nothing.
    */
  virtual void output_iteration(const Algorithm& algo) const;

  /** \brief Output information about a just completed algorithm.
    *
    * The default just does nothing.
    */
  virtual void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

  //@}

  /** @name Journal file access */
  //@{

  /** \brief Set a smart pointer to the journal file.
    */
  virtual void set_journal_out(const ostream_ptr_t& journal_out);

  /** \brief Get the smart pointer to the journal file.
   */
  const ostream_ptr_t& get_journal_out() const;

  /** \brief Return a reference to a <tt>std::ostream</tt> to be used to output debug information 
    * and the like.
    */
  virtual std::ostream& journal_out() const;

  //@}

private:

#ifndef DOXYGEN_COMPILE
  ostream_ptr_t   journal_out_;
#endif

  // not defined and not to be called
  AlgorithmTracker();
  
};	// end class AlgorithmTracker

}	// end namespace IterationPack 

#endif // ALGORITHM_TRACK_H
