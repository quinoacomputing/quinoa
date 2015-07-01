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

#include <iomanip>

#include "IterationPack_AlgorithmTrackTesting.hpp"
#include "IterationPack_Algorithm.hpp"

namespace IterationPack {

void AlgorithmTrackTesting::output_iteration(const Algorithm& algo) const {
  std::ostream &o = journal_out();
  o
    << "\ntrack.output_iteration(algo) called for the iteration k = "
    << algo.state().k() << std::endl;
  o
    << "\nStep times for the current iteration (in seconds):\n";
  const int n = algo.num_steps();
  std::vector<double> step_times(n+1);
  algo.get_step_times_k(0,&step_times[0]);
  o << "  step_id:time = ";
  for( int l = 0; l < n+1; ++l )
    o << " " << (l+1) << ":" << step_times[l];
  o << "\n  total time = " << step_times[n] << std::endl;
  
}

void AlgorithmTrackTesting::output_final(const Algorithm& algo, EAlgoReturn algo_return) const {
  char algo_return_name[6][50] =
    {
      "TERMINATE_TRUE"
      ,"TERMINATE_FALSE"
      ,"MAX_ITER_EXCEEDED"
      ,"MAX_RUN_TIME_EXCEEDED"
      ,"INTERRUPTED_TERMINATE_TRUE"
      ,"INTERRUPTED_TERMINATE_FALSE"
    };
  std::ostream &o = journal_out();
  o << "\ntrack.output_final(algo,algo_return) called for the iteration k = "
    << algo.state().k() << " and algo_return = " << algo_return_name[algo_return]
    << std::endl;
  o << "Timing (in seconds) statistics for step 0 : ";
  double total, average, min, max, percent;
  algo.get_final_step_stats(0,&total,&average,&min,&max,&percent);
  o << "total = " << total << ", average = " << average << ", min = " << min
    << ", max = " << max << ", percent = " << percent << std::endl;
}

}	// end namespace IterationPack 
