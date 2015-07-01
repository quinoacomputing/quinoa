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

#ifndef ACT_SET_STATS_H
#define ACT_SET_STATS_H

#include "MoochoPack_Types.hpp"

namespace MoochoPack {

/** \brief Class for storing statistics about the changes in the active set
  * of an SQP algorithm
  */
class ActSetStats {
public:

  // Public types

  /// Set to this value if a statistic is not known.
  enum { NOT_KNOWN = -1 };

  // Public interface

  /// Construct all unknowns
  ActSetStats()
    : num_active_(NOT_KNOWN), num_adds_(NOT_KNOWN), num_drops_(NOT_KNOWN)
    , num_active_indep_(NOT_KNOWN), num_adds_indep_(NOT_KNOWN), num_drops_indep_(NOT_KNOWN)
  {}

  /// Initialize the statistics
  void set_stats(
    int num_active, int num_adds, int num_drops
    ,int num_active_indep, int num_adds_indep, int num_drops_indep
    )
  {
    num_active_        = num_active;
    num_adds_          = num_adds;
    num_drops_         = num_drops;
    num_active_indep_  = num_active_indep;
    num_adds_indep_    = num_adds_indep;
    num_drops_indep_   = num_drops_indep;
  }

  /** \brief . */
  int num_active() const
  {
    return num_active_;
  }
  /** \brief . */
  int	num_adds() const
  {
    return num_adds_;
  }
  /** \brief . */
  int	num_drops() const
  {
    return num_drops_;
  }

  /** \brief . */
  int num_active_indep() const
  {
    return num_active_indep_;
  }
  /** \brief . */
  int	num_adds_indep() const
  {
    return num_adds_indep_;
  }
  /** \brief . */
  int	num_drops_indep() const
  {
    return num_drops_indep_;
  }

private:
  int num_active_;
  int	num_adds_;
  int	num_drops_;
  int num_active_indep_;
  int	num_adds_indep_;
  int	num_drops_indep_;

};	// end class ActSetStats

}	// end namespace MoochoPack

#endif	// ACT_SET_STATS_H
