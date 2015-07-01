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

#ifndef PLAYA_SERIAL_VECTOR_HPP
#define PLAYA_SERIAL_VECTOR_HPP

#include "PlayaDefs.hpp"
#include "PlayaSingleChunkVector.hpp"
#include "PlayaLoadableVector.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "Teuchos_Array.hpp"

namespace Playa
{
using namespace Teuchos;

template <class Scalar> class Vector;

/**
 * Playa implementation of a serial vector, implementing the LoadableVector
 * interface allowing an application to access elements. 
 * If created in SPMD, this will be replicated on
 * all processors.
 */
class SerialVector : public SingleChunkVector<double>,
                     public LoadableVector<double>,
                     public Describable
{
public:

  /** Construct with a vector space. */
  SerialVector(const VectorSpace<double>& vs);

  /** \name VectorBase interface */
  //@{
  /** Access to the space in which this vector lives */
  RCP<const VectorSpaceBase<double> > space() const {return vecSpace_.ptr();}
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

  /** \name Diagnostics */
  //@{
  /** */
  std::string description() const ;
  //@}



  /** \name Access through global indices */
  //@{
  /** get the batch of elements at the given global indices */
  void getElements(const int* globalIndices, int numElems,
    Array<double>& elems) const ;
  //@}

  /** */
  static const SerialVector* getConcrete(const Vector<double>& x);
  /** */
  static SerialVector* getConcrete(Vector<double>& x);

      
 /** \name Single chunk data access interface */
  //@{
  /** */
  virtual const double* dataPtr() const {return &(data_[0]);}
  /** */
  virtual double* dataPtr() {return &(data_[0]);}

  /** Size of the (single) chunk of data values */
  virtual int chunkSize() const {return dim_;}
  //@}



  
private:

  VectorSpace<double> vecSpace_;

  Array<double> data_;

  int dim_;
};
  
}


#endif
