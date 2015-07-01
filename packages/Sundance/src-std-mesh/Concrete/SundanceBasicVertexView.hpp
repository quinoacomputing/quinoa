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

#ifndef SUNDANCE_BASICVERTEXVIEW_H
#define SUNDANCE_BASICVERTEXVIEW_H

#include "SundanceDefs.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
using namespace Teuchos;

/**
 * VertexView is a read-only "view" of a cell's vertices, where the 
 * vertices are stored contiguously in a large master array. By working 
 * with views, we can greatly reduce the number of temporary arrays 
 * created during hashtable searches for existing vertex arrays.
 */
class VertexView
{
public:
  /** empty ctor, needed for storing VertexViews in Teuchos hashtables */
  VertexView() : base_(0), offset_(0), length_(0) {;}
  /** Construct a view into an array
   * \param bese pointer to the start of the master data array. By 
   * using double indirection, the master array can be resized or 
   * relocated and VertexViews can remain valid. 
   * \param  offset the index of the vertex subarray being viewed. 
   * \param length the number of vertices included in this view.
   */
  VertexView(int** base, int offset, int length)
    : base_(base), offset_(offset), length_(length) {;}

  /**
   * Return a hash code for the vertex set. 
   */
  int hashCode() const ;

  /** 
   * Test equality between two vertex sets. 
   * Two vertex sets are equal when their vertices are identical.
   */
  bool operator==(const VertexView& other) const ;

  /**
   * Write to a std::string
   */
  std::string toString() const ;


private:
  int** base_;
  int offset_;
  int length_;
};

/*
 * Two vertex sets are equal when their vertices are identical.
 * 
 */
inline bool VertexView::operator==(const VertexView& other) const
{
  /* For efficiency's sake, skip the test for equal lengths because we
   * can assume the caller is only comparing equal length vertex views */
  int* p = *base_ + offset_*length_;
  int* op = *(other.base_) + other.offset_*length_;

  for (int i=0; i<length_; i++)
  {
    if (p[i] != op[i]) return false;
  }
  return true;
}


}

namespace Teuchos
{
/** \relates VertexView */
inline int hashCode(const Sundance::VertexView& v) 
{return v.hashCode();}

/** \relates VertexView */
inline std::string toString(const Sundance::VertexView& v) 
{return v.toString();}
}

#endif
