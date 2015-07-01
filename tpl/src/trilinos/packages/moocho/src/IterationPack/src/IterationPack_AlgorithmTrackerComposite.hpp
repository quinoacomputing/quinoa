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

#ifndef ALGORITHM_TRACK_COMPOSITE_H
#define ALGORITHM_TRACK_COMPOSITE_H

#include <list>

#include "IterationPack_AlgorithmTracker.hpp"
#include "Teuchos_RCP.hpp"

namespace IterationPack {

/** \brief This class acts a composite container for other \c AlgorithmTracker objects.
 *
 * This class exposes a <tt>std::list<AlgorithmTracker*></tt> object and lets the client
 * manipulate the list.  It is up to the client to maintain this list.
 *
 * See the "Composite" pattern in "Design Patterns", Gama et al, 1995.
 */
class AlgorithmTrackerComposite : public AlgorithmTracker {
public:

  /** \brief . */
  typedef Teuchos::RCP<AlgorithmTracker>      track_ptr_t;
  /** \brief . */
  typedef std::list<track_ptr_t>                                    track_list_t;
  /** \brief . */
  AlgorithmTrackerComposite(const ostream_ptr_t& journal_out);
  /// Give access to the list of \c AlgorithmTracker object pointers.
  track_list_t& tracks();
  /** \brief . */
  const track_list_t& tracks() const;

  /**  @name Overridden from AlgorithmTracker */
  //@{

  /** \brief . */
  void initialize();
  /** \brief . */
  void output_iteration(const Algorithm& algo) const;
  /** \brief . */
  void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

  //@}

private:

#ifdef DOXYGEN_COMPILE
  AlgorithmTracker  *tracks;
#else
  track_list_t    tracks_;
#endif

};	// end class AlgorithmTrackerComposite

// ///////////////////////////////////
// Inline members

inline
AlgorithmTrackerComposite::track_list_t&
AlgorithmTrackerComposite::tracks()
{ 
  return tracks_;
}

inline
const AlgorithmTrackerComposite::track_list_t&
AlgorithmTrackerComposite::tracks() const
{ 
  return tracks_;
}

}	// end namespace IterationPack 

#endif	// ALGORITHM_TRACK_COMPOSITE_H
