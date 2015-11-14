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

// Use of Epetra_Operator
// This code should be run with at least two processes

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
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"

// auxiliary function to local an index in an integer vector

int find( int key, int vector[], int Length )
{
  for( int i=0 ; i<Length ; ++i ) {
    if( vector[i] == key ) return i;
  }
  return -1;
}

// ==================== //
// TriDiagonal Operator //
// -------------------- //

// NOTE: Will not work with IndexBase != 0
class TriDiagonalOperator : public Epetra_Operator
{

  public:

    // constructor
    TriDiagonalOperator( double diag_minus_one,
        double diag,
        double diag_plus_one,
        const Epetra_Map & Map) :
      Map_( Map ),
      diag_minus_one_(diag_minus_one),
      diag_(diag),
      diag_plus_one_(diag_plus_one)
  {
    // build the importer
    // Each local node will need the node+1 and node-1
    // (except for global node 0 and global node NumGlobalElemenets-1
    NumMyElements_ = Map_.NumMyElements();
    NumGlobalElements_ = Map_.NumGlobalElements();
    int* MyGlobalElements = new int[NumMyElements_];
    Map_.MyGlobalElements(MyGlobalElements);

    // count the nodes required from other processors
    // (this will be an upper bound of the # of required nodes
    // because I may count twice an external node)
    int count=0;
    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int globalIndex = MyGlobalElements[i];
      // no -1 node for the first node of the grid
      if( globalIndex>0 )
        if( Map.LID(globalIndex-1) == -1 ) ++count;
      // now +1 node for the last node of the grid
      if( globalIndex<NumGlobalElements_-1 )
        if( Map.LID(globalIndex+1) == -1 ) ++count;
      ++count;
    }

    // now allocate space for local nodes and external nodes
    // (an external node is a node required for the matrix-vector
    // product, but owned by another process)
    int Length = count;
    std::vector<int> ListOfNodes(Length);

    count=0;
    for( int i=0 ; i<NumMyElements_ ; ++i ) {
      int globalIndex = MyGlobalElements[i];
      // no -1 node for the first node of the grid
      if( globalIndex>0 ) {
        if( Map.LID(globalIndex-1) == -1 )
          if( find( globalIndex-1, ListOfNodes.data(), Length) == -1 ) {
            ListOfNodes[count] = globalIndex-1;
            ++count;
          }
      }
      // now +1 node for the last node of the grid
      if( globalIndex<NumGlobalElements_-1 ) {
        if( Map.LID(globalIndex+1) == -1 ) {
          if( find( globalIndex+1, ListOfNodes.data(), Length) == -1 ) {
            ListOfNodes[count] = globalIndex+1;
            ++count;
          }
        }
      }
      ListOfNodes[count] = globalIndex;
      ++count;
    }
    /*
       cout << "count = " << count << endl;
       for( int i=0 ; i<count ; i++ ) {
       cout << "ListOfNodes[" << i << "] = " << ListOfNodes[i] << endl;
       }
       */
    // create a Map defined using ListOfNodes

    ImportMap_ = new Epetra_Map(-1,count,ListOfNodes.data(),0,Map_.Comm());

    Importer_ = new  Epetra_Import(*ImportMap_,Map_);

    delete[] MyGlobalElements;

    return;

  }

    // application of the tridiagonal operator
    int Apply( const Epetra_MultiVector & X,
        Epetra_MultiVector & Y ) const
    {

      cout << X;

      // maybe some error checks on MultiVector Lenghts
      // for the future...

      Epetra_MultiVector Xext((*ImportMap_),X.NumVectors());

      // this will contain local nodes and the required extenal nodes
      Xext.Import(X,*Importer_,Insert);

      for( int i=0 ; i<X.MyLength() ; ++i ) {

        int globalRow = Map_.GID(i);
        int iMinusOne = (*ImportMap_).LID(globalRow-1);
        int iPlusOne = (*ImportMap_).LID(globalRow+1);

        printf("%d %d %d\n", globalRow, iMinusOne, iPlusOne);

        for( int vec=0 ; vec<X.NumVectors() ; ++vec ) {

          Y[vec][i] = diag_ * X[vec][i];

          if( iMinusOne != -1 )
            Y[vec][i] += diag_minus_one_ * Xext[vec][iMinusOne];

          if( iPlusOne != -1 )
            Y[vec][i] += diag_plus_one_ * Xext[vec][iPlusOne];

        }
      }

      return true;
    }

    // other function
    int SetUseTranspose( bool UseTranspose)
    {
      return(0);
    }

    int ApplyInverse( const Epetra_MultiVector & X,
        Epetra_MultiVector & Y ) const
    {
      return 0;
    }

    double NormInf() const
    {
      return( abs(diag_) + abs(diag_minus_one_) + abs(diag_plus_one_) );
    }

    const char * Label () const
    {
      return "TriDiagonalOperator";
    }

    bool UseTranspose() const
    {
      return false;
    }

    bool HasNormInf () const
    {
      return true;
    }


    const Epetra_Comm & Comm() const
    {
      return( Map_.Comm() );
    }

    const Epetra_Map & OperatorDomainMap() const
    {
      return( Map_ );
    }

    const Epetra_Map & OperatorRangeMap() const
    {
      return( Map_ );
    }


  private:

    Epetra_Map Map_;
    int NumMyElements_;
    int NumGlobalElements_;
    double diag_minus_one_;   // value in the sub-diagonal
    double diag_;             // value in the diagonal
    double diag_plus_one_;    // value in the super-diagonal
    Epetra_Import *Importer_;
    Epetra_Map *ImportMap_;

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

  // global dimension of the problem, could be any positive number
  int NumGlobalElements( 5 );

  // linear decomposition (for simplicity, could be general)
  Epetra_Map Map(NumGlobalElements,0,Comm );

  // define two vectors based on Map
  Epetra_Vector x(Map);
  Epetra_Vector y(Map);
  int NumMyElements = Map.NumMyElements();
  Epetra_IntSerialDenseVector MyGlobalElements(NumMyElements);
  Map.MyGlobalElements( MyGlobalElements.Values() );

  // x is a linear function, Laplace applied to it
  // should be zero except for the boundary nodes
  for( int i=0 ; i<NumMyElements ; ++i )
    x[i] = 1.0*MyGlobalElements[i];

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

