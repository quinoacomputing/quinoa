/*! @HEADER */
/*
************************************************************************

                CTrilinos:  C interface to Trilinos
                Copyright (2009) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

3. Neither the name of the Corporation nor the names of the
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact M. Nicole Lemaster (mnlemas@sandia.gov)

************************************************************************
*/
/*! @HEADER */

#ifdef HAVE_CTRILINOS_EXPERIMENTAL
// ----------   Includes   ----------
#include <iostream>
#include "Epetra_CrsMatrix.h"
#include "CNOX_Interface.hpp"

#include "EpetraExt_RowMatrixOut.h"

//-----------------------------------------------------------------------------
CNOX_Interface::CNOX_Interface(int* nelems, double* statevector,
         const LOCA::ParameterVector& pVector_, const Epetra_Comm& comm_,
         void* blackbox_res_, void* blackbox_prec_,
         void (*residualFunction_)(double *, double *, int, void *),
         void (*precFunction_)(double *, double *, int, double*, void *)) :
  N(*nelems),
  comm(comm_),
  pVector(pVector_),
  blackbox_res(blackbox_res_),
  blackbox_prec(blackbox_prec_),
  residualFunction(residualFunction_),
  precFunction(precFunction_)
{

  globalMap = Teuchos::rcp(new Epetra_Map(-1, N, 0, comm));

  solution = Teuchos::rcp(new Epetra_Vector(Copy, *globalMap, statevector));
}

CNOX_Interface::~CNOX_Interface()
{
}

bool CNOX_Interface::computeF(const Epetra_Vector& x, Epetra_Vector& F, FillType flag)
{
  F.PutScalar(0.0);
  residualFunction(x.Values(), F.Values(), N, blackbox_res);

  return true;
}

void CNOX_Interface::setParameters(const LOCA::ParameterVector& params)
{
  pVector = params;
}

void CNOX_Interface::printSolution(const Epetra_Vector& x, double conParam)
{
    cout << setprecision(5) << "Solution at: " << conParam << "  is:  "
    << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " "
    << x[4] << " " << x[5] << " ... " << x[N-1] << endl;
}

Teuchos::RCP<Epetra_Vector> CNOX_Interface::getVector() const
{ return solution;}

// Preconditioner is 2 steps. Only computePreconditioner is given the state,
//  which we store in solution.
bool CNOX_Interface::computePreconditioner(const Epetra_Vector& x,
                                              Epetra_Operator& Prec,
                                              Teuchos::ParameterList* p)
{
  (*solution) = x;
  return true;
}

// ApplyInverse is the action of the preconditioner. The state was already sent above.
int CNOX_Interface::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
   precFunction(Y(0)->Values(), X(0)->Values(), N, solution->Values(), blackbox_prec);

   return 0;
}
//-----------------------------------------------------------------------------
#endif //HAVE_CTRILINOS_EXPERIMENTAL
