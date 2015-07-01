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

#include "SundanceOut.hpp"
#include "PlayaTabs.hpp"
#include "SundanceVectorAssemblyKernel.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_TimeMonitor.hpp"


using namespace Sundance;
using namespace Sundance;
using namespace Sundance;
using namespace Teuchos;
using namespace Playa;
using std::setw;
using std::endl;
      
VectorAssemblyKernel::VectorAssemblyKernel(
  const Array<RCP<DOFMapBase> >& dofMap,
  const Array<RCP<Array<int> > >& isBCIndex,
  const Array<int>& lowestLocalIndex,
  Array<Vector<double> >& b,
  bool partitionBCs,
  int verb
  )
  : VectorFillingAssemblyKernel(dofMap, isBCIndex, lowestLocalIndex,
    b, partitionBCs, verb)
{}


void VectorAssemblyKernel::fill(
  bool isBC, 
  const IntegralGroup& group,
  const RCP<Array<double> >& localValues) 
{
  Tabs tab0;
  SUNDANCE_MSG1(verb(), tab0 << "in VectorAssemblyKernel::fill()");

  TEUCHOS_TEST_FOR_EXCEPT(!group.isOneForm());

  bool useCofacets = group.usesMaximalCofacets();

  if (group.isOneForm())
  {
    insertLocalVectorBatch(isBC, useCofacets, 
      group.testID(), group.testBlock(), group.mvIndices(),
      *localValues);
  }

  SUNDANCE_MSG1(verb(), tab0 << "done VectorAssemblyKernel::fill()");
}
  

void VectorAssemblyKernel:: prepareForWorkSet(
  const Array<Set<int> >& requiredTests,
  const Array<Set<int> >& /* requiredUnks */,
  RCP<StdFwkEvalMediator> mediator)
{
  Tabs tab0;
  SUNDANCE_MSG1(verb(), tab0 
    << "in VectorAssemblyKernel::prepareForWorkSet()");
  IntegrationCellSpecifier intCellSpec = mediator->integrationCellSpec();
  buildLocalDOFMaps(mediator, intCellSpec, requiredTests);
  SUNDANCE_MSG1(verb(), tab0 << "done VectorAssemblyKernel::prepareForWorkSet()");
}
