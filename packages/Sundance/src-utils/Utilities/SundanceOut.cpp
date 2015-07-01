/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "Teuchos_Array.hpp"

using namespace Sundance;
using namespace Teuchos;
using namespace std;


namespace Sundance
{

void writeTable(std::ostream& os, const Tabs& tab,
  const Array<double>& a, int cols)
{
  int rows = a.size() / cols;

  for (int i=0; i<rows; i++)
  {
    os << tab << setw(10) << i << ":";
    for (int j=0; j<cols; j++) 
      os << setw(12) << setprecision(6) << a[i*cols+j];
    os << std::endl;
  }
  int n = a.size() - rows * cols;
  if (n==0) return ;
  os << tab << setw(10) << rows << ":" ;
  for (int j=0; j<n; j++) 
    os << setw(12) << setprecision(6) << a[rows*cols+j];
  os << std::endl;
}


void writeTable(std::ostream& os, const Tabs& tab,
  const Array<int>& a, int cols)
{
  int rows = a.size() / cols;

  for (int i=0; i<rows; i++)
  {
    os << tab << setw(10) << i << ":";
    for (int j=0; j<cols; j++) 
      os << setw(10) << a[i*cols+j];
    os << std::endl;
  }
  int n = a.size() - rows * cols;
  if (n==0) return ;
  os << tab << setw(10) << rows << ":" ;
  for (int j=0; j<n; j++) 
    os << setw(10) << a[rows*cols+j];
  os << std::endl;
}


}







