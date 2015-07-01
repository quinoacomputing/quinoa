/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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


#ifndef PLAYA_MATRIXMARKETREADER_HPP
#define PLAYA_MATRIXMARKETREADER_HPP

#include <iostream>
#include <iomanip>
#include <vector>
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_RowMatrixOut.h"
#include "PlayaLinearOperatorImpl.hpp"
#include "PlayaEpetraMatrix.hpp"
#include "PlayaEpetraVectorSpace.hpp"
#include "Epetra_SerialComm.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

using namespace EpetraExt;


namespace Playa
{
  /**
   * This function reads a real matrix from a Matrix Market file and creates
   * an Epetra_CrsMatrix wrapped in a Playa LinearOperator object.
   *
   * The hard work of reading the file is done by an EpetraExt function. The
   * rest is a matter of setting up RCPs and creating a VectorSpace for the
   * domain and range of the operator.
   */
  LinearOperator<double> readMatrixMarketFile(const string& fileName)
  {
    /* Call the epetraext function that reads the matrix and forms
     * an Epetra matrix object. */
#ifdef HAVE_MPI
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif
    Epetra_CrsMatrix *A;
    int ierr = MatrixMarketFileToCrsMatrix(fileName.c_str(),comm,A,false,true);

    TEUCHOS_TEST_FOR_EXCEPTION(ierr<0, std::runtime_error,
			       "MatrixMarket reader failed with ierr=" << ierr);

    /* Capture the matrix pointer into an RCP */
    RCP<Epetra_CrsMatrix> Aptr = rcp(A);

    /* Get the map for the matrix's domain space */
    const Epetra_Map& dom = A->DomainMap();

    /* Form a vector space. Class EpetraVectorSpace is a subtype of class
     * Playa::VectorSpaceBase<double>. */
    EpetraVectorSpace* vs = new EpetraVectorSpace(rcp(new Epetra_Map(dom)));
    /* Capture into an RCP */
    RCP<const VectorSpaceBase<double> > vsb = rcp(vs);
    /* Capture into a vector space handle */
    VectorSpace<double> vecSpace(vsb);	

    /* Form the LinearOperator object. Note: the domain and range spaces are
     * identical here, which limits us to square matrices. This should be
     * fixed at some point. */
    RCP<LinearOperatorBase<double> > LOAb = rcp(new EpetraMatrix(Aptr,vecSpace,vecSpace));
    LinearOperator<double> LOA(LOAb);
    return LOA; 
  }

  void writeMatrixMarketFile(const string& filename,
			     const LinearOperator<double>& A,
			     const string& comment="")
  {
    RCP<const Epetra_CrsMatrix> crs = EpetraMatrix::getConcretePtr(A);

    int ierr = EpetraExt::RowMatrixToMatrixMarketFile(filename.c_str(),
						      *crs,
						      0,
						      comment.c_str());
    TEUCHOS_TEST_FOR_EXCEPT(ierr<0);
  }
}

#endif
