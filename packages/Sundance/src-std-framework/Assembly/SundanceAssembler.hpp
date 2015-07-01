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

#ifndef SUNDANCE_ASSEMBLER_H
#define SUNDANCE_ASSEMBLER_H

#include "SundanceDefs.hpp"
#include "PlayaLoadableVector.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorType.hpp"
#include "Teuchos_HashSet.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaIncrementallyConfigurableMatrixFactory.hpp"
#include "PlayaCollectivelyConfigurableMatrixFactory.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceMesh.hpp"
#include "SundanceEvalContext.hpp"
#include "SundanceIntegrationCellSpecifier.hpp"
#include "SundanceComputationType.hpp"


namespace Sundance
{
using namespace Teuchos;
using namespace Playa;

class EquationSet;
class EvaluatableExpr;
class EvalManager;
class EvalVector;
class DiscreteSpace;
class DiscreteFunction;
class CellFilter;
class DOFMapBase;
class IntegralGroup;
class StdFwkEvalMediator;
class AssemblyKernelBase;

typedef std::set<int> ColSetType;

/** 
 * 
 */
class Assembler 
{
public:
  /** */
  Assembler(
    const Mesh& mesh, 
    const RCP<EquationSet>& eqn,
    const Array<VectorType<double> >& rowVectorType,
    const Array<VectorType<double> >& colVectorType,
    bool partitionBCs);


  /** */
  Assembler(
    const Mesh& mesh, 
    const RCP<EquationSet>& eqn);
      
  /** */
  const Array<RCP<DOFMapBase> >& rowMap() const 
    {return rowMap_;}

  /** */
  const Array<RCP<DOFMapBase> >& colMap() const 
    {return colMap_;}

  /** */
  const Array<RCP<DiscreteSpace> >& solutionSpace() const 
    {return externalColSpace_;}

  /** */
  const Array<RCP<DiscreteSpace> >& rowSpace() const 
    {return externalRowSpace_;}

  /** */
  VectorSpace<double> solnVecSpace() const ;

  /** */
  VectorSpace<double> rowVecSpace() const ;

  /** */
  const Array<RCP<Set<int> > >& bcRows() const { return bcRows_; }

  /** Allocate, but do not fill, the matrix */
  Playa::LinearOperator<double> allocateMatrix() const ;

  /** */
  void assemble(Playa::LinearOperator<double>& A,
    Array<Vector<double> >& b) const ;

  /** */
  void assembleSensitivities(Playa::LinearOperator<double>& A,
    Array<Vector<double> >& b) const ;


  /** */
  void assemble(Array<Vector<double> >& b) const ;

  /** */
  void evaluate(double& value,
    Array<Vector<double> >& gradient) const ;

  /** */
  void evaluate(double& value) const ;

  /** */
  static int& workSetSize() ;

      
  /** */
  void getGraph(int br, int bc,
    Array<int>& graphData,
    Array<int>& rowPtrs,
    Array<int>& nnzPerRow) const ;
      
  /** */
  void incrementalGetGraph(int br, int bc, 
    IncrementallyConfigurableMatrixFactory* mf) const ;

  /** flushes all configuration , so that it enforces the reassemble of the matrix*/
  void flushConfiguration() const;

  /** */
  static int& numAssembleCalls() {static int rtn=0; return rtn;}

  /** */
  static bool& matrixEliminatesRepeatedCols() {static bool x = false; return x;}

  /** */
  const RCP<EquationSet>& eqnSet() const 
    {return eqn_;}

  /** */
  int maxWatchFlagSetting(const std::string& param) const ;

  /** */
  static Time& assemblyTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("assembly"); 
      return *rtn;
    }

  /** */
  static Time& configTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("matrix config"); 
      return *rtn;
    }
  
  /** */
  static Time& fillTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("matrix/vector fill"); 
      return *rtn;
    }
  

private:

  /** */
  void init(const Mesh& mesh, const RCP<EquationSet>& eqn);

  /** */
  bool detectInternalBdry(int cellDim, const CellFilter& filter) const ;

  /** */
  void displayEvaluationResults(
    const EvalContext& context, 
    const EvaluatableExpr* evalExpr, 
    const Array<double>& constantCoeffs, 
    const Array<RCP<EvalVector> >& vectorCoeffs) const ;

  /** */
  void assemblyLoop(const ComputationType& compType,
    RCP<AssemblyKernelBase> kernel) const ;

  /** */
  bool matNeedsConfiguration() const;

  /** */
  void configureMatrix(LinearOperator<double>& A,
    Array<Vector<double> >& b) const ;

  /** */
  void configureVector(Array<Vector<double> >& b) const ;

  /** */
  void configureMatrixBlock(int br, int bc, 
    LinearOperator<double>& A) const ;

  /** */
  void configureVectorBlock(int br, Vector<double>& b) const ;

  /** */
  Array<Array<int> > findNonzeroBlocks() const ;

  /** */
  IntegrationCellSpecifier whetherToUseCofacets(
    const Array<RCP<IntegralGroup> >& groups,
    const EvaluatableExpr* ee,
    bool isMaximalCell,
    int verb) const ;

  
  
  /** */
  static int defaultWorkSetSize() {static int rtn=100; return rtn;}

  bool partitionBCs_;

  mutable int numConfiguredColumns_;

  Mesh mesh_;

  RCP<EquationSet> eqn_;

  Array<RCP<DOFMapBase> > rowMap_;

  Array<RCP<DOFMapBase> > colMap_;

  Array<RCP<DiscreteSpace> > externalRowSpace_;

  Array<RCP<DiscreteSpace> > externalColSpace_;

  Array<RCP<DiscreteSpace> > privateRowSpace_;

  Array<RCP<DiscreteSpace> > privateColSpace_;

  Array<RCP<Set<int> > > bcRows_;

  Array<RegionQuadCombo> rqc_;

  Map<ComputationType, Array<EvalContext> > contexts_;

  Map<ComputationType, Array<int> > skipRqc_;

  Array<int> isBCRqc_;

  Array<int> isInternalBdry_;

  Map<ComputationType, Array<Array<RCP<IntegralGroup> > > > groups_;

  Array<RCP<StdFwkEvalMediator> > mediators_;

  Map<ComputationType, Array<const EvaluatableExpr*> > evalExprs_;

  RCP<EvalManager> evalMgr_;

  Array<RCP<Array<int> > > isBCRow_;

  Array<RCP<Array<int> > > isBCCol_;

  Array<RCP<std::set<int> > > remoteBCCols_;

  Array<int> lowestRow_;

  Array<int> lowestCol_;

  Array<VectorType<double> > rowVecType_;

  Array<VectorType<double> > colVecType_;

  Map<int, int> testIDToBlockMap_;

  Map<int, int> unkIDToBlockMap_;

  Map<int, int> fixedParamIDToVectorNumber_;

  Map<ComputationType, Array<IntegrationCellSpecifier> > rqcRequiresMaximalCofacets_;

  /** Cached reference to the previously assembled matrix
   *  A null value signals that matrix must be assembled */
  mutable LinearOperator<double> cachedAssembledMatrix_;
};

}



#endif
