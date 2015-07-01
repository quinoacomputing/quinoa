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

#include "euclid_Helpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>
#include <fstream>
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Ifpack_Hypre.h"
#include "Epetra_MpiComm.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <sstream>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;

Epetra_CrsMatrix::Epetra_CrsMatrix* newCrsMatrix(int N){
  int ierr = 0;
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  
  RCP<Epetra_Map>            Map;
  // pointer to the matrix to be created
  Epetra_CrsMatrix*      Matrix;/*
  // container for parameters
  Teuchos::ParameterList GaleriList;
  int nx = N * Comm.NumProc();
  int ny = N * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);

  Map = rcp(Galeri::CreateMap("Cartesian2D", Comm, GaleriList));
  Matrix   = Galeri::CreateCrsMatrix("Laplace2D", Map.get(), GaleriList);
 */ 
  Map = rcp(new Epetra_Map(N*Comm.NumProc(), N, 0, Comm));
  Matrix = new Epetra_CrsMatrix(Copy, *Map, 3);

  int ilower = N*Comm.MyPID(); int iupper = N*(Comm.MyPID()+1)-1;
  for(int i = ilower; i <= iupper; i++){
    if(i != 0 && i != N*Comm.NumProc()-1){
      Teuchos::Array<int> indices; 
      Teuchos::Array<double> values;
      values.resize(3); indices.resize(3);
      indices[0] = i-1;
      indices[1] = i;
      indices[2] = i+1;
      values[0] = -1.0;
      values[1] = 4.0;
      values[2] = -1.0;
      Matrix->InsertGlobalValues(i, 3, &values[0], &indices[0]);
    } else if(i != N*Comm.NumProc()-1){
      Teuchos::Array<int> indices;
      Teuchos::Array<double> values;
      values.resize(2); indices.resize(2);
      indices[0] = 0;
      indices[1] = 1;
      values[0] = 4.0;
      values[1] = -1.0;
      Matrix->InsertGlobalValues(0, 2, &values[0], &indices[0]);
    } else{
      Teuchos::Array<int> indices;
      Teuchos::Array<double> values;
      values.resize(2); indices.resize(2);
      indices[0] = i-1;
      indices[1] = i;
      values[0] = -1.0;
      values[1] = 4.0;
      Matrix->InsertGlobalValues(i, 2, &values[0], &indices[0]);
    }
  }

  ierr += Matrix->FillComplete();
  if(ierr != 0){
    printf("\n\nMatrix was not created properly!\n\n");
  }
  return Matrix;
}


bool EquivalentVectors(Epetra_MultiVector &X, Epetra_MultiVector &Y, const double tol){
  bool retVal = true;
  string fileName;
  stringstream out;
  out << "Failure." << X.Comm().MyPID();
  fileName = out.str();
  ofstream myfile;
  
  int num_vectors = X.NumVectors();
  if(Y.NumVectors() != num_vectors){
    printf("Multivectors do not have same number of vectors.\n");
    return false;
  }
  
  for(int j = 0; j < num_vectors; j++){
    if(X.MyLength() != Y.MyLength()){
      printf("Vectors are not same size on local processor.\n");
      return false;
    }
    for(int i = 0; i < Y.MyLength(); i++){
      if(fabs(X[j][i] - Y[j][i]) > tol){
        if(retVal == true) {
          myfile.open(fileName.c_str());
          myfile << setw(8) << "Vector" << setw(8) << "Entry" << setw(12) << "X" << setw(12) << "Y\n";
        }
        myfile << setw(8) << j << setw(8) << i << setw(12) << X[j][i] << setw(12) << Y[j][i] << endl;    
        retVal = false;      
      }
    }
  }
  if(retVal == false) myfile.close();
  Teuchos::Array<int> vals; vals.resize(X.Comm().NumProc());
  int my_vals[1]; my_vals[0] = (int)retVal;
  X.Comm().GatherAll(my_vals, &vals[0], 1);
  for(int i = 0; i < X.Comm().NumProc(); i++){
    if(vals[i] == false){
      retVal = false;
    }
  }
  if(!retVal && !Y.Comm().MyPID()) printf("See Failure file for details\n");
  return retVal;
}

