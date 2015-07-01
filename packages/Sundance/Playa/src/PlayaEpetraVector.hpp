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

#ifndef PLAYA_EPETRAVECTOR_HPP
#define PLAYA_EPETRAVECTOR_HPP

#include "PlayaDefs.hpp"
#include "PlayaPrintable.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaLoadableVector.hpp"
#include "PlayaSingleChunkVector.hpp"
#include "Epetra_FEVector.h"
#include "Epetra_Vector.h"
#include "PlayaEpetraVectorSpace.hpp"


namespace Playa
{
using Teuchos::RCP;
/**
 * Playa::VectorBase wrapper for Epetra_Vector
 */
class EpetraVector : public SingleChunkVector<double>,
                     public LoadableVector<double>,
                     public Printable,
                     public Teuchos::Describable
{
public:

  /** Construct with an Epetra vector space. */
  EpetraVector(const VectorSpace<double>& vs);

  /** Construct with an Epetra vector space
      and an existing Epetra vector. */
  EpetraVector(const VectorSpace<double>& vs,
    const RCP<Epetra_Vector>& vec);


  /** \name VectorBase interface */
  //@{
  /** Return the space in which this vector lives */
  RCP< const VectorSpaceBase<double> > space() const {return vecSpace_.ptr();}
  //@}

  /** \name IndexableVector interface */
  //@{
  /** read the element at the given global index */
  virtual const double& operator[](int globalIndex) const 
    {return getElement(globalIndex);}

  /** writable access to the element at the given global index */
  virtual double& operator[](int globalIndex) ;
  //@}

  /** \name LoadableVector interface */
  //@{
  /** set a single element */
  void setElement(int globalIndex, const double& value);

  /** add to a single element */
  void addToElement(int globalIndex, const double& value);

  /** set a group of elements */
  void setElements(int numElems, const int* globalIndices, 
    const double* values);


  /** add to a group of elements */
  void addToElements(int numElems, const int* globalIndices, 
    const double* values);

  /** */
  void finalizeAssembly();
  //@}

  /** \name AccessibleVector interface */
  //@{
  /** */
  const double& getElement(int globalIndex) const ;

  /** */
  void getElements(const int* globalIndices, int numElems,
    Teuchos::Array<double>& elems) const ;
  //@}

  /** \name Printable interface */
  //@{
  /** Print to a stream */
  void print(std::ostream& os) const ;
  //@}
      

  /** */
  const RCP<Epetra_Vector>& epetraVec() const 
    {return epetraVec_;}

  /** */
  RCP<Epetra_Vector>& epetraVec() {return epetraVec_;}

  /** Get a read-only Epetra_Vector */
  static const Epetra_Vector& getConcrete(const Playa::Vector<double>& tsfVec);
  /** Get a read-write Epetra_Vector */
  static Epetra_Vector& getConcrete(Playa::Vector<double>& tsfVec);
  /** Get a read-write Epetra_Vector pointer */
  static Epetra_Vector* getConcretePtr(Playa::Vector<double>& tsfVec);

  
  /** */
  virtual void update(const double& alpha, const VectorBase<double>* other,
    const double& gamma);
    

  /** */
  virtual void update(
    const double& alpha, const VectorBase<double>* x,
    const double& beta, const VectorBase<double>* y,
    const double& gamma) ;

  /** */
  virtual double dot(const VectorBase<double>* other) const ;

  /** */
  virtual double norm2() const ;
protected:    

  /** \name Single chunk data access interface */
  //@{
  /** */
  virtual const double* dataPtr() const {return &(epetraVec_->operator[](0));}
  /** */
  virtual double* dataPtr() {return &(epetraVec_->operator[](0));}

  /** Size of the (single) chunk of data values */
  virtual int chunkSize() const {return numLocalElements_;}
  //@}

private:

  VectorSpace<double> vecSpace_;

  RCP<Epetra_Vector> epetraVec_;

  int numLocalElements_;
  
};
  
}


#endif
