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

#include "PlayaTabs.hpp"
#include "SundanceOut.hpp"
#include "SundanceRegionQuadCombo.hpp"

using namespace Sundance;
using namespace Teuchos;


RegionQuadCombo::RegionQuadCombo()
  : id_(-1), domain_(), quad_(), paramCurve_() , watch_()
{;}

RegionQuadCombo::RegionQuadCombo(const RCP<CellFilterStub>& domain,
  const RCP<QuadratureFamilyStub>& quad,
  const ParametrizedCurve& paramCurve ,
  const WatchFlag& watch)
  : id_(getID(domain, quad,watch)), domain_(domain), quad_(quad),
  paramCurve_(paramCurve) , watch_(watch)
{;}


int RegionQuadCombo::getID(const RCP<CellFilterStub>& domain,
  const RCP<QuadratureFamilyStub>& quad,
  const WatchFlag& watch)
{
  RegTriple p(domain, quad, watch);

  if (!domainAndQuadToIDMap().containsKey(p))
    {
      int id = topID();
      domainAndQuadToIDMap().put(p, id);
    }
  return domainAndQuadToIDMap().get(p);
}

string RegionQuadCombo::toString() const
{
  TeuchosOStringStream oss;
  Tabs tabs;

  oss << "Integration Region" << std::endl;
  {
    oss << tabs << "cell filter=" << domain_->description() << std::endl;
    oss << tabs << "quadrature rule=" << quad_->description() << std::endl;
    oss << tabs << "watchpoint=[" << watch().name() << "]" << std::endl;
  }
  return oss.str();
}

Sundance::Map<RegTriple, int>& RegionQuadCombo::domainAndQuadToIDMap()
{
  static Sundance::Map<RegTriple, int> rtn = Sundance::Map<RegTriple, int>();
  return rtn;
}





