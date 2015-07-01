/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
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

#ifndef SUNDANCE_ARRAYOFTUPLES_H
#define SUNDANCE_ARRAYOFTUPLES_H

#include "Teuchos_Array.hpp"

namespace Sundance
{
  using namespace Teuchos;

  /** 
   * Class ArrayOfTuples packs an heterogeneous array of tuples into 
   * a single 1D array. 
   */
  template<class T>
    class ArrayOfTuples
    {
    public:
      /** Empty ctor */
      ArrayOfTuples();

      /** Constructor specifying the size of each tuple, but not the 
       * number of tuples */
      ArrayOfTuples(int tupleSize);

      /** Constructor specifying both the size and number of the tuples */
      ArrayOfTuples(int numTuples, int tupleSize);

      /** Returns the number of tuples */
      int length() const {return numTuples_;}

      /** Returns the size of the tuples */
      int tupleSize() const {return tupleSize_;}

      /** Change the number of tuples */
      void resize(int newSize) {
        numTuples_ = newSize;
        data_.resize(numTuples_*tupleSize_);
      }

      /** Change the number and size of the tuples */
      void resize(int newSize, int newTupleSize) 
        {
        numTuples_ = newSize;
        tupleSize_ = newTupleSize;
        data_.resize(numTuples_*tupleSize_);
        }

      /** Reserve memory for a number of tuples */
      void reserve(int newSize) {data_.reserve(newSize*tupleSize_);}

      /** Specify the size of the tuples */
      void setTupleSize(int tupleSize) {tupleSize_ = tupleSize;}

      /** Get the j-th entry in the i-th tuple */
      const T& value(int i, int j) const {return data_[i*tupleSize_+j];}

      /** Get the j-th entry in the i-th tuple */
      T& value(int i, int j) {return data_[i*tupleSize_+j];}

      /** Append a new tuple to the array */
      void append(const Array<T>& x);

      /** Append a new tuple to the array */
      void append(const T* x, int n);

    private:

      int numTuples_;
      int tupleSize_;

      Array<T> data_;

    };

  template<class T> inline ArrayOfTuples<T>::ArrayOfTuples()
    : numTuples_(0), tupleSize_(0), data_()
    {;}

  template<class T> inline ArrayOfTuples<T>::ArrayOfTuples(int tupleSize)
    : numTuples_(0), tupleSize_(tupleSize), data_()
    {;}

  template<class T> inline ArrayOfTuples<T>::ArrayOfTuples(int numTuples, int tupleSize)
    : numTuples_(numTuples), tupleSize_(tupleSize), data_(numTuples*tupleSize)
    {;}

  template<class T> inline void ArrayOfTuples<T>::append(const Array<T>& x)
    {
      for (int i=0; i<x.length(); i++)
        {
          data_.append(x[i]);
        }
      numTuples_++;
    }

  template<class T> inline void ArrayOfTuples<T>::append(const T* x, int n)
    {
      for (int i=0; i<n; i++)
        {
          data_.append(x[i]);
        }
      numTuples_++;
    }


}

#endif
