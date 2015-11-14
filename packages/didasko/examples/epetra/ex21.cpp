// @HEADER
// ***********************************************************************
//
//                      Didasko Tutorial Package
//                 Copyright (2005) Sandia Corporation
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
// Questions about Didasko? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ***********************************************************************
// @HEADER

// Use of Epetra_Operator.
// This code must be run with one process

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"

// ==================== //
// TriDiagonal Operator //
// -------------------- //

class TriDiagonalOperator : public Epetra_Operator
{

  public:

    // constructor
    TriDiagonalOperator( double diag_minus_one,
        double diag,
        double diag_plus_one,
        Epetra_Map & Map) :
      Map_(Map),
      diag_minus_one_(diag_minus_one),
      diag_(diag),
      diag_plus_one_(diag_plus_one)
  {}

    // application of the tridiagonal operator
    int Apply( const Epetra_MultiVector & X,
        Epetra_MultiVector & Y ) const
    {
      int Length = X.MyLength();

      // maybe some error checks on MultiVector Lenghts
      // for the future...

      for( int vec=0 ; vec<X.NumVectors() ; ++vec ) {

        // one-dimensional problems here
        if( Length == 1 ) {
          Y[vec][0] = diag_ * X[vec][0];
          break;
        }

        // more general case (Lenght >= 2)

        // first row
        Y[vec][0] = diag_ * X[vec][0] + diag_plus_one_ * X[vec][1];

        // intermediate rows
        for( int i=1 ; i<Length-1 ; ++i ) {
          Y[vec][i] = diag_ * X[vec][i] + diag_plus_one_ * X[vec][i+1]
            + diag_minus_one_ * X[vec][i-1];
        }
        // final row
        Y[vec][Length-1] = diag_ * X[vec][Length-1]
          + diag_minus_one_ * X[vec][Length-2];
      }

      return true;
    }

    // other function
    int SetUseTranspose( bool UseTranspose)
    {
      return(-1); // not implemented
    }

    int ApplyInverse( const Epetra_MultiVector & X,
        Epetra_MultiVector & Y ) const
    {
      return(-1); // not implemented
    }

    double NormInf() const
    {
      return(abs(diag_) + abs(diag_minus_one_) + abs(diag_plus_one_));
    }

    const char * Label () const
    {
      return("TriDiagonalOperator");
    }

    bool UseTranspose() const
    {
      return(false);
    }

    bool HasNormInf () const
    {
      return(true);
    }


    const Epetra_Comm & Comm() const
    {
      return(Map_.Comm());
    }

    const Epetra_Map & OperatorDomainMap() const
    {
      return(Map_);
    }

    const Epetra_Map & OperatorRangeMap() const
    {
      return(Map_);
    }


  private:

    Epetra_Map Map_;
    double diag_minus_one_;   // value in the sub-diagonal
    double diag_;             // value in the diagonal
    double diag_plus_one_;    // value in the super-diagonal

};

// =========== //
// main driver //
// ----------- //

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if( Comm.NumProc() != 1 ) {
    if( Comm.MyPID() == 0 ) {
      cerr << "This is mono-process example\n"
        << "Please run with one processo only\n";
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  // global dimension of the problem, could be any positive number
  int NumGlobalElements( 5 );

  // linear decomposition (for simplicity, could be general)
  Epetra_Map Map(NumGlobalElements,0,Comm );

  // define two vectors based on Map
  Epetra_Vector x(Map);
  Epetra_Vector y(Map);
  x.PutScalar(1.0);

  // define a linear operator, as previously defined in class
  // TriDiagonalOperator

  TriDiagonalOperator TriDiagOp(-1.0,2.0,-1.0,Map);

  TriDiagOp.Apply(x,y);

  cout << x;
  cout << y;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure Didasko with:\n"
      "--enable-epetra");

  return 0;
}
#endif

