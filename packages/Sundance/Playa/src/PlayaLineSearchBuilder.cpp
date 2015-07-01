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


#include "PlayaLineSearchBuilder.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PlayaTabs.hpp"
#include "PlayaOut.hpp"
#include "PlayaSimpleBacktracking.hpp"

namespace Playa
{

RCP<LineSearchBase> 
LineSearchBuilder::createLineSearch(const ParameterList& params,
  int verb)
{
  Tabs tab(0);
  PLAYA_MSG1(verb, tab << "LineSearchBuilder::createLineSearch()");
  Tabs tab1;
  PLAYA_MSG2(verb, tab1 << "params=" << params);
  
  TEUCHOS_TEST_FOR_EXCEPTION(params.name() != "Line Search",
    std::runtime_error, 
    "LineSearchBuilder::getLineSearch() expected parameter list named "
    "\"Line Search\", got name [" << params.name() << "]");

  const std::string& lsType = getParameter<string>(params, "Type");

  RCP<LineSearchBase> ls;

  if (lsType=="Simple Backtracking")
  {
    PLAYA_MSG2(verb, tab1 << "found Simple Backtracking LS");
    ls = rcp(new SimpleBacktracking(params));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(ls.get()==0, 
    std::runtime_error, 
    "LineSearchBuilder::getLineSearch() could not construct a valid line "
    "search object from parameter list " << params);
    
  return ls;
}

}
