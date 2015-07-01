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

#ifndef SUNDANCE_VECTORFILLINGASSEMBLYKERNEL_H
#define SUNDANCE_VECTORFILLINGASSEMBLYKERNEL_H

#include "SundanceDefs.hpp"
#include "SundanceAssemblyKernelBase.hpp"
#include "PlayaLoadableVector.hpp"
#include "SundanceMapBundle.hpp"

namespace Sundance
{
using namespace Teuchos;


/**
 * VectorFillingAssemblyKernel provides a common implementation of
 * multivector-filling capabilities needed by matrix-vector assembly, 
 * multivector assembly, and functional gradient evaluation. 
 */
class VectorFillingAssemblyKernel : public AssemblyKernelBase
{
public:
  /**
   * Ctor takes several arguments:
   * \param dofMap is an array of DOFMap ptrs, one for each block 
   *
   * \param isBCIndex is an array of ptrs to arrays of ints (bools). The value 
   * (*isBCIndex[b])[d] indicates whether dof #d in block #b is or is not 
   * an essential BC dof. 
   *
   * \param lowestLocalIndex stores the lowest locally-owned DOF index for each 
   * block 
   * 
   * \param b multivector to be filled
   *
   * \param partitionBC whether dirichlet BCs are stored in a separate block
   *
   * \param verb verbosity level
   */
  VectorFillingAssemblyKernel(
    const Array<RCP<DOFMapBase> >& dofMap,
    const Array<RCP<Array<int> > >& isBCIndex,
    const Array<int>& lowestLocalIndex,
    Array<Vector<double> >& b,
    bool partitionBCs,
    int verb
    );

  /** */
  virtual ~VectorFillingAssemblyKernel(){;}

protected:

  /** 
   * Insert local vector values into a global vector
   */
  void insertLocalVectorBatch(bool isBCRqc,
    bool useCofacetCells,
    const Array<int>& funcID,  
    const Array<int>& funcBlock,
    const Array<int>& mvIndices,  
    const Array<double>& localValues) const ;

  /** 
   * Build fast lookup tables of DOFs for local cells.
   */
  void buildLocalDOFMaps(
    const RCP<StdFwkEvalMediator>& mediator,
    IntegrationCellSpecifier intCellSpec,
    const Array<Set<int> >& requiredFuncs)  ;

  

protected:
  const MapBundle& mapBundle() const {return mapBundle_;}

private:
  /** Multivector to be filled */
  Array<Vector<double> > b_;   
  /** Loadable view of the multivector to be filled */
  Array<Array<RCP<LoadableVector<double> > > > vec_;
  /** The map bundle collects together several DOF-mapping data structures */
  mutable MapBundle mapBundle_;

};


}



#endif
