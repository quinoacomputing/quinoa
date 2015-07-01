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


#include "PlayaEpetraVectorType.hpp"
#include "PlayaEpetraVectorSpace.hpp"
#include "PlayaEpetraGhostImporter.hpp"
#include "PlayaEpetraMatrixFactory.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"
#include "PlayaOut.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_RefCountPtr.hpp"
#include "PlayaEpetraMatrix.hpp"

using namespace Playa;
using namespace Teuchos;

EpetraVectorType::EpetraVectorType()
{;}


RCP<const VectorSpaceBase<double> > 
EpetraVectorType::createSpace(int /*dimension*/,
  int nLocal,
  const int* localIndices,
  const MPIComm& comm) const
{
#ifdef HAVE_MPI
  Epetra_MpiComm epComm(comm.getComm());
#else
  Epetra_SerialComm epComm;
#endif

  TEUCHOS_TEST_FOR_EXCEPTION(nLocal < 0, std::runtime_error, "negative vector size n=" << nLocal);

	RCP<Epetra_Map> map = rcp(new Epetra_Map(-1, nLocal,
      (int*) localIndices,
      0, epComm));

	return rcp(new EpetraVectorSpace(map));
}

RCP<GhostImporter<double> > 
EpetraVectorType::createGhostImporter(const VectorSpace<double>& space,
                                      int nGhost,
                                      const int* ghostIndices) const
{
  const EpetraVectorSpace* p 
    = dynamic_cast<const EpetraVectorSpace*>(space.ptr().get());

  TEUCHOS_TEST_FOR_EXCEPTION(p==0, std::runtime_error,
                     "non-epetra vector space [" << space.description() << "] given as "
                     "argument to EpetraVectorType::createGhostImporter()");

  return rcp(new EpetraGhostImporter(p->epetraMap(), nGhost, ghostIndices));
  
}

RCP<MatrixFactory<double> >
EpetraVectorType::createMatrixFactory(const VectorSpace<double>& domain,
                                      const VectorSpace<double>& range) const
{
  RCP<const EpetraVectorSpace> pd 
    = rcp_dynamic_cast<const EpetraVectorSpace>(domain.ptr());

  RCP<const EpetraVectorSpace> pr 
    = rcp_dynamic_cast<const EpetraVectorSpace>(range.ptr());


  TEUCHOS_TEST_FOR_EXCEPTION(pd.get()==0, std::runtime_error, 
                     "incompatible domain space given to "
                     "EpetraVectorType::createMatrix()");

  TEUCHOS_TEST_FOR_EXCEPTION(pr.get()==0, std::runtime_error, 
                     "incompatible range space given to "
                     "EpetraVectorType::createMatrix()");

  //  RCP<SingleScalarTypeOp<double> > A = rcp(new EpetraMatrix(pd, pr));

  return rcp(new EpetraMatrixFactory(pd, pr));
}






