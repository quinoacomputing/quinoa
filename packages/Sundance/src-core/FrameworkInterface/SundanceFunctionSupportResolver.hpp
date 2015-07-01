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

#ifndef SUNDANCE_FUNCTION_SUPPORT_RESOLVER_H
#define SUNDANCE_FUNCTION_SUPPORT_RESOLVER_H

#include "SundanceDefs.hpp"
#include "SundanceExpr.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceOrderedHandle.hpp"
#include "SundanceCellFilterStub.hpp"


namespace Sundance
{

class SumOfIntegrals;
class SumOfBCs;
class CommonFuncDataStub;

/** */
class FunctionSupportResolver
{
public:
  /** */
  FunctionSupportResolver(
    const Expr& eqns,
    const Expr& bcs,
    const Array<Expr>& vars,
    const Array<Expr>& unks,
    const Expr& unkParams,
    const Expr& params,
    const Array<Expr>& fixedFields,
    bool isVariational);
  
  
  /** \name Getting information about functions */
  //@{
  /** Returns the number of variational function blocks */
  int numVarBlocks() const {return varFuncs_.size();}

  /** Returns the number of unknown function blocks */
  int numUnkBlocks() const {return unkFuncs_.size();}

  /** Returns the number of unknown parameters */
  int numUnkParams() const {return unkParams_.size();}

  /** Returns the number of fixed parameters */
  int numFixedParams() const {return fixedParams_.size();}

  /** Returns the number of variational functions in this block */
  int numVars(int block) const {return varFuncs_[block].size();}

  /** Returns the number of unk functions in this block */
  int numUnks(int block) const {return unkFuncs_[block].size();}

  /** Returns the number of variational function IDs in this block.
   * This will differ from the number of variational functions in cases
   * where a vector field uses a single vector-valued basis rather
   * than scalar bases for each component. */
  int numVarIDs(int block) const {return varFuncData_[block].size();}

  /** Returns the number of unk function IDs in this block.
   * This will differ from the number of unknown functions in cases
   * where a vector field uses a single vector-valued basis rather
   * than scalar bases for each component. */ 
  int numUnkIDs(int block) const {return unkFuncData_[block].size();}

  /** Returns the data for the i-th variational function in block b */
  RCP<const CommonFuncDataStub> varFuncData(int b, int i) const {return varFuncData_[b][i];}

  /** Returns the data for the i-th unknown function in block b */
  RCP<const CommonFuncDataStub> unkFuncData(int b, int i) const {return unkFuncData_[b][i];}

  /** Returns the i-th unknown parameter */
  const Expr& unkParam(int i) const {return unkParams_[i];}

  /** Determine whether a given func ID is listed as a 
   * variational function in this equation set */
  bool hasVarID(int fid) const 
    {return varIDToBlockMap_.containsKey(fid);}

  /** Determine whether a given func ID is listed as a unk function 
   * in this equation set */
  bool hasUnkID(int fid) const 
    {return unkIDToBlockMap_.containsKey(fid);}

  /** Determine whether a given func ID is listed as a unk parameter 
   * in this equation set */
  bool hasUnkParamID(int fid) const 
    {return unkParamIDToReducedUnkParamIDMap_.containsKey(fid);}

  /** Determine whether a given func ID is listed as a fixed parameter 
   * in this equation set */
  bool hasFixedParamID(int fid) const 
    {return fixedParamIDToReducedFixedParamIDMap_.containsKey(fid);}

  /** get the block number for the variational function having the
   * specified unreduced funcID */
  int blockForVarID(int varID) const ;

  /** get the block number for the unknown function having the
   * specified unreduced funcID */
  int blockForUnkID(int unkID) const ;
  //@}


  /** \name Finding the functions that appear on regions */
  //@{
  /** Returns the variational functions that appear explicitly
   * on the d-th region */
  const Set<int>& varsOnRegion(int d) const 
    {return varsOnRegions_.get(regions_[d]);}

  /** Returns the unknown functions that appear explicitly on the
   * d-th region. */
  const Set<int>& unksOnRegion(int d) const 
    {return unksOnRegions_.get(regions_[d]);}

  /** Returns the variational functions that 
   * appear in BCs on the d-th region.
   * We can use this information to tag certain rows as BC rows */
  const Set<int>& bcVarsOnRegion(int d) const 
    {return bcVarsOnRegions_.get(regions_[d]);}

  /** Returns the unknown functions that appear in BCs on the d-th region.
   * We can use this information to tag certain columns as BC
   * columns in the event we're doing symmetrized BC application */
  const Set<int>& bcUnksOnRegion(int d) const 
    {return bcUnksOnRegions_.get(regions_[d]);}

  /** Returns the reduced variational functions that appear explicitly
   * on the d-th region */
  const Array<Set<int> >& reducedVarsOnRegion(const OrderedHandle<CellFilterStub>& r) const 
    {return reducedVarsOnRegions_[indexForRegion(r)];}

  /** Returns the reduced unknown functions that appear explicitly on the
   * d-th region. */
  const Array<Set<int> >& reducedUnksOnRegion(const OrderedHandle<CellFilterStub>& r) const 
    {return reducedUnksOnRegions_[indexForRegion(r)];}
  //@}

      


  /** \name Transforming between unreduced and reduced function IDs */
  //@{
  /** get the reduced ID for the variational function having the
   * specified unreduced funcID */
  int reducedVarID(int varID) const ;

  /** get the reduced ID for the unknown 
   * function having the given funcID */
  int reducedUnkID(int unkID) const ;

  /** get the reduced ID for the unk parameter
   * having the given funcID */
  int reducedUnkParamID(int unkID) const ;

  /** get the reduced ID for the fixed parameter
   * having the given funcID */
  int reducedFixedParamID(int unkID) const ;

  /** get the unreduced funcID for a variational function
   * as specified by a reduced ID and block index */
  int unreducedVarID(int block, int reducedVarID) const 
    {return unreducedVarID_[block][reducedVarID];}

  /** get the unreduced funcID for an unknown function
   * as specified by a reduced ID and block index */
  int unreducedUnkID(int block, int reducedUnkID) const 
    {return unreducedUnkID_[block][reducedUnkID];}

  /** get the unreduced funcID for an unknown parameter
   * as specified by a reduced ID  */
  int unreducedUnkParamID(int reducedUnkParamID) const 
    {return unreducedUnkParamID_[reducedUnkParamID];}

  /** get the unreduced funcID for a fixed parameter
   * as specified by a reduced ID */
  int unreducedFixedParamID(int reducedFixedParamID) const 
    {return unreducedFixedParamID_[reducedFixedParamID];}

  /** Return the map from fixed param ID to reduced fixed param ID */
  const Map<int, int>& fixedParamIDToReducedFixedParamIDMap() const
    {return fixedParamIDToReducedFixedParamIDMap_;}
  
  //@}


  /** \name Finding integration regions for the equation set */
  //@{
  /** Returns the number of regions on which pieces of the equation
   * or BCs are defined. */
  int numRegions() const {return regions_.size();}
      
  /** Returns the d-th region for this equation set */
  const RCP<CellFilterStub>& region(int d) const 
    {return regions_[d].ptr();}

  /** Returns the index of the given region */
  int indexForRegion(const OrderedHandle<CellFilterStub>& region) const ;

  /** Whether a region has BCs */
  bool isBCRegion(int d) const 
    {return bcVarsOnRegions_.containsKey(regions_[d]);}


  /** Return the set of regions on which the specified 
   * test func appears. */
  const Set<OrderedHandle<CellFilterStub> >& 
  regionsForTestFunc(int unreducedTestID) const ;
      
  /** Return the set of regions on which the specified 
   * unknown func appears */
  const Set<OrderedHandle<CellFilterStub> >& 
  regionsForUnkFunc(int unreducedUnkID) const ;
  //@}


  /**
   * Flatten a spectral expression into a list of its coefficients
   */
  Expr flattenSpectral(const Expr& input) const ;
  /**
   * Flatten a spectral expression into a list of its coefficients
   */
  Array<Expr> flattenSpectral(const Array<Expr>& input) const ;

  /** Whether essential BCs are present */
  bool hasBCs() const ;

  /** Access to integrals */
  const SumOfIntegrals* integralSum() const {return integralSum_;}

  /** Access to BCs */
  const SumOfBCs* bcSum() const {return bcSum_;}

  /** */
  const Array<Expr>& unks() const {return unkFuncs_;}

  /** */
  const Array<Expr>& vars() const {return varFuncs_;}

  /** */
  const Array<Expr>& fixedFields() const {return fixedFields_;}

  /** */
  const Expr& fixedParams() const {return fixedParams_;}

  /** */
  const Expr& unkParams() const {return unkParams_;}

  /** */
  const Set<int>& varFuncSet() const {return varFuncSet_;}
  /** */
  const Set<int>& unkFuncSet() const {return unkFuncSet_;}
  /** */
  const Set<int>& unkParamSet() const {return unkParamSet_;}
  /** */
  const Set<int>& fixedParamSet() const {return fixedParamSet_;}
  
      
private:

  /** */
  Expr eqns_;

  /** */
  Expr bcs_;

  /** */
  const SumOfIntegrals* integralSum_;

  /** */
  const SumOfBCs* bcSum_;

  /** */
  Set<int> varFuncSet_;

  /** */
  Set<int> unkFuncSet_;

  /** */
  Set<int> unkParamSet_;

  /** */
  Set<int> fixedParamSet_;

  /** */
  Array<OrderedHandle<CellFilterStub> > regions_;

  /** */
  Map<OrderedHandle<CellFilterStub>, int> regionToIndexMap_;

  /** */
  Map<OrderedHandle<CellFilterStub>, Set<int> > varsOnRegions_;

  /** */
  Map<OrderedHandle<CellFilterStub>, Set<int> > unksOnRegions_;

  /** */
  Array<Array<Set<int> > > reducedVarsOnRegions_;

  /** */
  Array<Array<Set<int> > > reducedUnksOnRegions_;

  /** */
  Map<OrderedHandle<CellFilterStub>, Set<int> > bcVarsOnRegions_;

  /** */
  Map<OrderedHandle<CellFilterStub>, Set<int> > bcUnksOnRegions_;

  /** */
  Map<int, Set<OrderedHandle<CellFilterStub> > > testToRegionsMap_;

  /** */
  Map<int, Set<OrderedHandle<CellFilterStub> > > unkToRegionsMap_;

  /** var function data for this equation set */
  Array<Array<RCP<const CommonFuncDataStub> > > varFuncData_;

  /** unknown function data for this equation set */
  Array<Array<RCP<const CommonFuncDataStub> > > unkFuncData_;

  /** var functions for this equation set */
  Array<Expr> varFuncs_;

  /** unknown functions for this equation set */
  Array<Expr> unkFuncs_;

  /** fixed functions for this equation set */
  Array<Expr> fixedFields_;

  /** The point in function space about which the equations
   * are linearized */
  Array<Expr> unkLinearizationPts_;

  /** unknown parameters for this equation set */
  Expr unkParams_;

  /** fixed parameters for this equation set */
  Expr fixedParams_;

  /** map from variational function funcID to that function's
   * position in list of var functions */
  Array<Map<int, int> > varIDToReducedIDMap_;

  /** map from unknown function funcID to that function's
   * position in list of unk functions */
  Array<Map<int, int> > unkIDToReducedIDMap_;

  /** map from unknown param funcID to that param's
   * position in list of unk params */
  Map<int, int> unkParamIDToReducedUnkParamIDMap_;

  /** map from fixed param funcID to that param's
   * position in list of fixed params */
  Map<int, int> fixedParamIDToReducedFixedParamIDMap_;

  /** map from variational function funcID to that function's
   * position in list of var blocks */
  Map<int, int> varIDToBlockMap_;

  /** map from unknown function funcID to that function's
   * position in list of unk blocks */
  Map<int, int> unkIDToBlockMap_;

  /** Map from (block, unreduced var ID) to reduced ID */
  Array<Array<int> > reducedVarID_;

  /** Map from (block, unreduced unk ID) to reduced ID */
  Array<Array<int> > reducedUnkID_;

  /** Map from unreduced unk ID to reduced ID */
  Array<int> reducedUnkParamID_;

  /** Map from unreduced fixed ID to reduced ID */
  Array<int> reducedFixedParamID_;

  /** Map from (block, reduced varID) to unreduced varID */
  Array<Array<int> > unreducedVarID_;

  /** Map from (block, reduced unkID) to unreduced unkID */
  Array<Array<int> > unreducedUnkID_;

  /** Map from reduced unkParamID to unreduced unkParamID */
  Array<int> unreducedUnkParamID_;

  /** Map from reduced fixedParamID to unreduced fixedParamID */
  Array<int> unreducedFixedParamID_;

  /** Flag indicating whether this equation set is nonlinear */
  bool isNonlinear_;
      
  /** Flag indicating whether this equation set is 
   * a variational problem */
  bool isVariationalProblem_;
};

/** */
RCP<const CommonFuncDataStub> getSharedFunctionData(const FuncElementBase* f);

}
 
#endif
