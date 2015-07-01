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

#ifndef SUNDANCE_ASSEMBLYKERNELBASE_H
#define SUNDANCE_ASSEMBLYKERNELBASE_H

#include "SundanceDefs.hpp"
#include "SundanceDOFMapBase.hpp"
#include "SundanceEquationSet.hpp"
#include "SundanceDiscreteSpace.hpp"
#include "SundanceDiscreteFunction.hpp"
#include "SundanceIntegralGroup.hpp"
#include "SundanceElementIntegral.hpp"
#include "SundanceGrouperBase.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceStdFwkEvalMediator.hpp"
#include "SundanceEvaluatableExpr.hpp"
#include "PlayaLoadableVector.hpp"
#include "PlayaLoadableMatrix.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaVectorDecl.hpp"
#include "PlayaVectorType.hpp"
#include "Teuchos_HashSet.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaIncrementallyConfigurableMatrixFactory.hpp"
#include "PlayaCollectivelyConfigurableMatrixFactory.hpp"

namespace Sundance
{
using namespace Teuchos;

class StdFwkEvalMediator;
class IntegralGroup;

/** 
 * AssemblyKernelBase abstracts the operations that must be done in an
 * assembly loop. Regardless of whether the assembly loop is doing
 * matrix/vector fill, vector fill, functional/gradient evaluation,
 * or functional evaluation, the assembly loop will involve
 * <ul>
 * <li> A preprocessing step before the main assembly loop
 * <li> A preprocessing step before each work set is started
 * <li> A fill step after each integral group is done
 * <li> A postprocessing step after the main assembly loop is done
 * </ul>
 * The first of these is done by the subclass constructor. The others
 * are done using the pure virtual functions of this class.
 *
 * It is assumed that any data structures to be filled -- such as a matrix,
 * a vector, or simply a number -- are stored internally in the assembly
 * kernel subclass, and that they persist between preprocessing and fill 
 * calls.
 */
class AssemblyKernelBase
{
public:
  /** */
  AssemblyKernelBase(int verb) : verb_(verb) {;}

  /** */
  virtual ~AssemblyKernelBase(){;}

  /**  
   * Do preprocessing steps needed before integrating the current
   * work set. 
   *
   * The default implementation does nothing. 
   */
  virtual void prepareForWorkSet(
    const Array<Set<int> >& requiredTests,
    const Array<Set<int> >& requiredUnks,
    RCP<StdFwkEvalMediator> mediator) {;}

  /** 
   * Adds the results of the current integral
   * group into the assembly results. 
   * \param isBC whether the current group is a replace-style 
   * boundary condition
   * \param group the current integral group
   * \param localValues the results of integrating the current integral group
   */
  virtual void fill(bool isBC,
    const IntegralGroup& group,
    const RCP<Array<double> >& localValues) = 0 ;  

  /** 
   * Hook to do any finalization steps after the main assembly loop, 
   * for example, doing an all-reduce on locally computed functional values. 
   * The default implementation does nothing. */
  virtual void postLoopFinalization() {;}

  /** verbosity level */
  int verb() const {return verb_;}

  /** set verbosity level.
   * (This function needs to be virtual because certain subclasses need specialized
   * implementations that propagate verbosity to children 
  */
  virtual void setVerb(int verb) {verb_=verb;}

private:
  int verb_;
};



}




#endif
