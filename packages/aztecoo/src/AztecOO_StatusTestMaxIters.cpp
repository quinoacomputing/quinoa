
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "AztecOO_StatusTestMaxIters.h"

AztecOO_StatusTestMaxIters::AztecOO_StatusTestMaxIters(int MaxIters) 
  : AztecOO_StatusTest() {

  if (MaxIters < 1) MaxIters_ = 1;
    
  MaxIters_ = MaxIters;
  status_ = Unchecked;
}

AztecOO_StatusType
AztecOO_StatusTestMaxIters::CheckStatus(int CurrentIter, 
                                        Epetra_MultiVector * CurrentResVector, 
                                        double CurrentResNormEst,
                                        bool SolutionUpdated)
{
  (void)CurrentResVector;
  (void)CurrentResNormEst;
  (void)SolutionUpdated;
  status_ = Unconverged;
  Niters_ = CurrentIter;
  if (Niters_ >= MaxIters_)
    status_ = Failed;
  return status_;
}


std::ostream& AztecOO_StatusTestMaxIters::Print(std::ostream& stream, int indent) const {

  for (int j = 0; j < indent; j ++)
    stream << ' ';
  PrintStatus(stream, status_);
  stream << "Number of Iterations = ";
  stream << Niters_;
  stream << ((Niters_<MaxIters_) ? " < " : ((Niters_==MaxIters_) ? " = " : " > "));
  stream << MaxIters_;
  stream << std::endl;
 return stream;
}
