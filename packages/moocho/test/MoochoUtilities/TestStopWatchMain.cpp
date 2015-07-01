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

#include "StopWatchPack_stopwatch.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

int main(int argc, char* argv[])
{

  using std::cout;
  using std::endl;
  using StopWatchPack::stopwatch;
  
  bool success = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  try {
    
    stopwatch timer;
    
    // Finding min resolution.
    double min_resolution = 0.0; // in seconds
    int total_num_calls   = 0;
    {
      cout	<< "\n*** Measuring miminum resolution.\n";
      timer.start();
      double last_time = timer.read();
      const int max_num_samples = 20;
      int num_samples = 0;
      int num_calls = 0;
      while( num_samples < max_num_samples ) {
        double time = timer.read();
        num_calls++;
        if( time - last_time > 0.0 ) {
          cout	<< "time_diff = " << time - last_time
                << ", num_calls = " << num_calls << endl;
          min_resolution += time - last_time;
          ++total_num_calls;
          last_time = time;
          num_calls = 0;
          num_samples++;
        }
      }
      min_resolution /= total_num_calls;
    }
    
    std::cerr << "Minimum stopwatch resolution = " << min_resolution << " sec\n";
    
    // Finding increasing resolution.
    {
      cout	<< "\n*** Measuring increasing resolution.\n";
      timer.start();
      double start_time = timer.read(), last_time = start_time;
      const int max_num_samples = 20;
      int num_samples = 0;
      int num_calls = 0;
      while( num_samples < max_num_samples ) {
        double time = timer.read();
        num_calls++;
        if( time - last_time > 0.0 ) {
          cout	<< "time = " << time - start_time
                << ", num_calls = " << num_calls << endl;
          last_time = time;
          num_calls = 0;
          num_samples++;
        }
      }
    }

  } // end try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);
  
  if(success)
    cout << "\nEnd Result: TEST PASSED" << std::endl;
  
  return 0;
  
}
