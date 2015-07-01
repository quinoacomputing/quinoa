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

#ifndef SUNDANCE_MATRIXVECTORASSEMBLYKERNEL_H
#define SUNDANCE_MATRIXVECTORASSEMBLYKERNEL_H

#include "SundanceDefs.hpp"
#include "SundanceVectorFillingAssemblyKernel.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * MatrixVectorAssemblyKernel does assembly of a matrix and vector
 */
class MatrixVectorAssemblyKernel : public VectorFillingAssemblyKernel
{
public:
  /** */
  MatrixVectorAssemblyKernel(
    const Array<RCP<DOFMapBase> >& rowMap,
    const Array<RCP<Array<int> > >& isBCRow,
    const Array<int>& lowestLocalRow,
    const Array<RCP<DOFMapBase> >& colMap,
    const Array<RCP<Array<int> > >& isBCCol,
    const Array<int>& lowestLocalCol,
    LinearOperator<double> A,
    Array<Vector<double> > b,
    bool partitionBCs,
    int verb)
    : VectorFillingAssemblyKernel(rowMap, isBCRow, lowestLocalRow, 
      b, partitionBCs, verb),
      mat_(rowMap.size()),
      cmb_(colMap, isBCCol, lowestLocalCol, partitionBCs, verb)
    {
      init(rowMap, colMap, A, partitionBCs);
    }

  /** */
  void prepareForWorkSet(
    const Array<Set<int> >& requiredTests,
    const Array<Set<int> >& requiredUnks,
    RCP<StdFwkEvalMediator> mediator) ;

  /** */
  void fill(bool isBC,
    const IntegralGroup& group,
    const RCP<Array<double> >& localValues) ;

protected:

  /** */
  void init(
  const Array<RCP<DOFMapBase> >& rowMap,
  const Array<RCP<DOFMapBase> >& colMap,
  LinearOperator<double> A,
  bool partitionBCs);

  /** */
  void writeLSMs(int blockRow, int blockCol,
    bool useCofacetCells,
    int numTestNodes, 
    int nTestFuncs, 
    int testFuncIndex, 
    const Array<int>& rowDof,
    int numUnkNodes, 
    int nUnkFuncs, 
    int unkFuncIndex, 
    const Array<int>& colDof,
    const Array<double>& localValues) const ;

  /** */
  void insertLocalMatrixBatch(
    bool isBCRqc,
    bool useCofacetCells,
    const Array<int>& testID, 
    const Array<int>& testBlock, 
    const Array<int>& unkID,
    const Array<int>& unkBlock,
    const Array<double>& localValues) const ;

protected:
  const MapBundle& rmb() const {return mapBundle();}
  const MapBundle& cmb() const {return cmb_;}

private:
  LinearOperator<double> A_;
  Array<Array<LoadableMatrix<double>* > > mat_;
  mutable MapBundle cmb_;
};

}



#endif
