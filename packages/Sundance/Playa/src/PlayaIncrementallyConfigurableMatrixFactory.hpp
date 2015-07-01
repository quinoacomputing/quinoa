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

#ifndef PLAYA_INCREMENTALLYCONFIGURABLEMATRIXFACTORY_HPP
#define PLAYA_INCREMENTALLYCONFIGURABLEMATRIXFACTORY_HPP

#include "PlayaDefs.hpp"

namespace Playa
{
  /** 
   * Class IncrementallyConfigurableMatrixFactory provides an abstract 
   * interface for row-at-a-time configuration of matrix factories.
   */
  class IncrementallyConfigurableMatrixFactory
  {
  public:
    /** Virtual dtor */
    virtual ~IncrementallyConfigurableMatrixFactory(){;}

    /** Initialize a set of nonzero elements in the matrix's graph.
     * @param globalRowIndex the global index of the row to which these
     * elements belong.
     * @param nElemsToInsert the number of elements being inserted in this
     * step
     * @param globalColumnIndices array of column indices. Must 
     * be nElemsToInsert in length. 
     */
    virtual void initializeNonzerosInRow(int globalRowIndex,
                                         int nElemsToInsert,
                                         const int* globalColumnIndices) = 0 ;

    /** 
     * Initialize nonzeros in a batch of rows. 
     */
    virtual void initializeNonzeroBatch(int numRows, 
                                        int rowBlockSize,
                                        const int* globalRowIndices,
                                        int numColumnsPerRow,
                                        const int* globalColumnIndices,
                                        const int* skipRow);

    /** Finalize values of the matrix. This is a hook for any
     * implementation-dependent steps that must be done after
     * loading of elements. */
    virtual void finalize() = 0 ;

  private:
    
    
  };

  /* Default implementation of initializeElementBatch */
  inline void IncrementallyConfigurableMatrixFactory
  ::initializeNonzeroBatch(int numRows, 
                           int rowBlockSize,
                           const int* globalRowIndices,
                           int numColumnsPerRow,
                           const int* globalColumnIndices,
                           const int* skipRow)
  {
    int numRowBlocks = numRows/rowBlockSize;
    int row = 0;

    for (int rb=0; rb<numRowBlocks; rb++)
      {
        const int* cols = globalColumnIndices + rb*numColumnsPerRow;
        for (int r=0; r<rowBlockSize; r++, row++)
          {
            if (skipRow[row]) continue;
            initializeNonzerosInRow(globalRowIndices[row], 
                                    numColumnsPerRow, cols);
          }
      }
  }
}

#endif
