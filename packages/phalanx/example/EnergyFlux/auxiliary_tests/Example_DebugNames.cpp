// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Traits.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  GlobalMPISession mpi_session(&argc, &argv);

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of debug naming testing
    // *********************************************************************
    {
      cout << "\nTesting Debug naming scheme.\n";
      // Evaluation Types
      cout << typeAsString<MyTraits::Residual>() << endl;
      cout << typeAsString<MyTraits::Jacobian>() << endl;
      // Data Types
      cout << typeAsString<double>() << endl;
      cout << typeAsString<MyTraits::FadType>() << endl;
      cout << typeAsString< MyVector<double> >() << endl;
      cout << typeAsString< MyVector<MyTraits::FadType> >() << endl;
      cout << typeAsString< MyTensor<double> >() << endl;
      cout << typeAsString< MyTensor<MyTraits::FadType> >() << endl;
    }
    
    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    cout << "\nRun has completed successfully!\n" << endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const exception& e) {
    cout << "************************************************" << endl;
    cout << "************************************************" << endl;
    cout << "Exception Caught!" << endl;
    cout << "Error message is below\n " << e.what() << endl;
    cout << "************************************************" << endl;
  }
  catch (...) {
    cout << "************************************************" << endl;
    cout << "************************************************" << endl;
    cout << "Unknown Exception Caught!" << endl;
    cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
