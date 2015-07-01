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

#ifndef MOOCHO_TRACKER_XML_SUMMARY_HPP
#define MOOCHO_TRACKER_XML_SUMMARY_HPP

#include "IterationPack_Algorithm.hpp"
#include "IterationPack_AlgorithmTracker.hpp"
#include "MoochoPack_Types.hpp"

namespace MoochoPack {


using IterationPack::Algorithm;
using IterationPack::EAlgoReturn;

/** \brief This class outputs an XML summary file of the algorithm 
 *   results and performance
 */
class MoochoTrackerXMLSummary
  : public IterationPack::AlgorithmTracker
{
public:

  /// Construct with an output stream
  MoochoTrackerXMLSummary(
    const Teuchos::RCP<std::ostream> &journal_out
    ,const std::string xml_filename
    ,const std::string problem_name
    ,const std::string algorithm_description
    );

  /// Set the output stream for summary outputting
  //void set_output_stream(const ostream_ptr_t& o);

  /// Get the output stream for summary outputting.
  //const ostream_ptr_t& get_output_stream() const;

  /// Output a basic file (with failed status)
  //   that will be overwritten if there is no
  //   exception
  void output_pre_file() const;

  /** @name Overridden from AlgorithmTracker */
  //@{

  /** \brief . */
  void output_iteration(const Algorithm& algo) const;
  /** \brief . */
  void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;
  
  //@}

protected:

  /// Print the header to the output
  void open_problem_element( std::ostream& out, const Algorithm& algo) const;
  void close_problem_element( std::ostream& out) const;

private:

  mutable value_type obj_value_;
  mutable value_type c_norm_value_;

  std::string xml_filename_;
  std::string problem_name_;
  std::string algorithm_description_;

  // Not defined and not to be called
  MoochoTrackerXMLSummary();

};	// end class MoochoTrackerXMLSummary

}	// end namespace MoochoPack 

#endif	// RSQP_TRACK_SUMMARY_STD_H
