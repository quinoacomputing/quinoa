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

#include "PlayaEpetraVector.hpp"
#include "PlayaEpetraVectorSpace.hpp"
#include "Teuchos_Assert.hpp"

#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "PlayaVectorImpl.hpp"
#include "PlayaLinearOperatorImpl.hpp"
#endif

using Teuchos::RCP;
using namespace Playa;

EpetraVector::EpetraVector(const VectorSpace<double>& vs)
  : SingleChunkVector<double>(), 
    vecSpace_(vs),
    epetraVec_(), 
    numLocalElements_(vs.numLocalElements())
{
  const EpetraVectorSpace* epvs 
    = dynamic_cast<const EpetraVectorSpace*>(vs.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(epvs==0, std::runtime_error,
    "could not cast vector space to EpetraVectorSpace in "
    "EpetraVector ctor");

  epetraVec_ = rcp(new Epetra_Vector(*(epvs->epetraMap())));
}



EpetraVector
::EpetraVector(const VectorSpace<double>& vs,
  const RCP<Epetra_Vector>& vec)
  : SingleChunkVector<double>(), 
    vecSpace_(vs), 
    epetraVec_(vec), 
    numLocalElements_(vs.numLocalElements())
{
  const EpetraVectorSpace* epvs 
    = dynamic_cast<const EpetraVectorSpace*>(vs.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(epvs==0, std::runtime_error,
    "could not cast vector space to EpetraVectorSpace in "
    "EpetraVector ctor");
}





double& EpetraVector::operator[](int globalIndex) 
{
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  return (*epetraVec())[myMap.LID(globalIndex)];
}

void EpetraVector::setElement(int index, const double& value)
{
  int loc_index[1] = { index };
  epetraVec()->ReplaceGlobalValues(1, const_cast<double*>(&value), 
    loc_index);
}

void EpetraVector::addToElement(int index, const double& value)
{
//  cout << "adding (" << index << ", " << value << ")" << std::endl;
  int loc_index[1] = { index };
  epetraVec()->SumIntoGlobalValues(1, const_cast<double*>(&value), 
    loc_index);
}

const double& EpetraVector::getElement(int index) const 
{
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  return (*epetraVec())[myMap.LID(index)];
}

void EpetraVector::getElements(const int* globalIndices, int numElems,
  Teuchos::Array<double>& elems) const
{
  elems.resize(numElems);
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  RCP<const Epetra_Vector> epv = epetraVec();

  for (int i=0; i<numElems; i++)
  {
    elems[i] = (*epv)[myMap.LID(globalIndices[i])];
  }
}

void EpetraVector::setElements(int numElems, const int* globalIndices,
  const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  int ierr = vec->ReplaceGlobalValues(numElems, globalIndices, values);
  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, "ReplaceGlobalValues returned "
    "ierr=" << ierr << " in EpetraVector::setElements()");
}

void EpetraVector::addToElements(int numElems, const int* globalIndices,
  const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  int ierr = vec->SumIntoGlobalValues(numElems, globalIndices, values);
  TEUCHOS_TEST_FOR_EXCEPTION(ierr < 0, std::runtime_error, "SumIntoGlobalValues returned "
    "ierr=" << ierr << " in EpetraVector::addToElements()");
}

void EpetraVector::finalizeAssembly()
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  vec->GlobalAssemble();
}


void EpetraVector::print(std::ostream& os) const 
{
  epetraVec()->Print(os);
}


const Epetra_Vector& EpetraVector::getConcrete(const Playa::Vector<double>& tsfVec)
{
  const EpetraVector* epv 
    = dynamic_cast<const EpetraVector*>(tsfVec.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return *(epv->epetraVec());
}

Epetra_Vector& EpetraVector::getConcrete(Playa::Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return *(epv->epetraVec());
}


Epetra_Vector* EpetraVector::getConcretePtr(Playa::Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEUCHOS_TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
    "EpetraVector::getConcrete called on a vector that "
    "could not be cast to an EpetraVector");
  return epv->epetraVec().get();
}


void EpetraVector::update(const double& alpha, 
  const VectorBase<double>* other, const double& gamma)
{
  const EpetraVector* epv 
    = dynamic_cast<const EpetraVector*>(other);
  epetraVec_->Update(alpha, *(epv->epetraVec_), gamma);
}



void EpetraVector::update(
  const double& alpha, const VectorBase<double>* x,
  const double& beta, const VectorBase<double>* y,
  const double& gamma)
{
  const EpetraVector* epx 
    = dynamic_cast<const EpetraVector*>(x);
  const EpetraVector* epy
    = dynamic_cast<const EpetraVector*>(y);

  epetraVec_->Update(alpha, *(epx->epetraVec_), 
    beta, *(epy->epetraVec_), gamma);
}


double EpetraVector::dot(const VectorBase<double>* other) const 
{
  double rtn = 0.0;
  const EpetraVector* epv 
    = dynamic_cast<const EpetraVector*>(other);
  epetraVec_->Dot(*(epv->epetraVec_), &rtn);

  return rtn;
}

double EpetraVector::norm2() const 
{
  double rtn = 0.0;
  epetraVec_->Norm2(&rtn);

  return rtn;
}
