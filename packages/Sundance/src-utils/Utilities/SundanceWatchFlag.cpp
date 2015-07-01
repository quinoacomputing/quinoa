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

#include "SundanceWatchFlag.hpp"
#include "SundanceOut.hpp"

using namespace Sundance;


WatchFlag::WatchFlag(const std::string& name,
  const ParameterList& params)
  : name_(name), params_(rcp(new ParameterList(params)))
{
  if (name_.size() > 0) isActiveMap().put(name, true);
  else isActiveMap().put(name, false);
}

void WatchFlag::activate() 
{
  isActiveMap()[name()] = true;
}

void WatchFlag::deactivate() 
{
  isActiveMap()[name()] = false;
}

bool WatchFlag::isActive() const 
{
  return isActiveMap().get(name());
}

XMLObject WatchFlag::toXML() const 
{
  XMLObject xml("WatchFlag");
  xml.addAttribute("name", name());
  return xml;
}

int WatchFlag::param(const std::string& name) const 
{
  return params_->get<int>(name);
}


void WatchFlag::setParam(const std::string& name, int val)
{
  params_->set<int>(name, val);
}




RCP<ParameterList> WatchFlag::defaultParams()
{
  static RCP<ParameterList> rtn=rcp(new ParameterList());
  static bool first=true;
  if (first)
  {
    rtn->set<int>("evaluation", 0);
    rtn->set<int>("evaluator setup", 0);
    rtn->set<int>("discrete function evaluation", 0);
    rtn->set<int>("symbolic preprocessing", 0);
    rtn->set<int>("equation set setup", 0);
    rtn->set<int>("assembler setup", 0);
    rtn->set<int>("assembly loop", 0);
    rtn->set<int>("dof map setup", 0);
    rtn->set<int>("dof map access", 0);
    rtn->set<int>("eval mediator", 0);
    rtn->set<int>("integration setup", 0);
    rtn->set<int>("integration", 0);
    rtn->set<int>("integral transformation", 0);    
    rtn->set<int>("fill", 0);
    rtn->set<int>("matrix config", 0);
    rtn->set<int>("vector config", 0);
    rtn->set<int>("solve control", 0);
    first = false;
  }
  return rtn;
}
