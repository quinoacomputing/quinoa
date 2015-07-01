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

#ifndef SUNDANCE_REGIONQUADCOMBO_H
#define SUNDANCE_REGIONQUADCOMBO_H


#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_Utils.hpp"
#include "SundanceWatchFlag.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceOrderedTuple.hpp"
#include "SundanceOrderedHandle.hpp"
#include "SundanceParametrizedCurve.hpp"


namespace Sundance
{
using namespace Teuchos;

/** */
typedef OrderedTriple<OrderedHandle<CellFilterStub>,
                      OrderedHandle<QuadratureFamilyStub>,
                      WatchFlag> RegTriple;
/** 
 * Expressions may appear in more than one subregions of a problem,
 * for instance in an internal domain and also on a boundary. On
 * those different subregions, a given expression might be subject
 * to different sets of functional derivatives; thus, different
 * evaluation regions might have different sparsity patterns. 
 * It is therefore necessary to build and store
 * sparsity information on a region-by-region basis. 
 *
 * Class RegionQuadCombo is used as an identifier for regions. The
 * only thing it needs to do is to be useable as a key in a STL map.
 */
class RegionQuadCombo 
{
public:
  /** */
  RegionQuadCombo();

  /** */
  RegionQuadCombo(const RCP<CellFilterStub>& domain,
    const RCP<QuadratureFamilyStub>& quad,
    const ParametrizedCurve& paramCurve = ParametrizedCurve::returnDummyCurve(),
    const WatchFlag& watch = WatchFlag());


  /** */
  inline bool operator==(const RegionQuadCombo& other) const
    {return id_==other.id_;}

  /** */
  std::string toString() const ;

  /** */
  bool operator<(const RegionQuadCombo& other) const
    {return id_ < other.id_;}

  /** */
  const RCP<CellFilterStub>& domain() const {return domain_;}

  /** */
  const RCP<QuadratureFamilyStub>& quad() const 
    {return quad_;}

  /** */
  const WatchFlag& watch() const 
    {return watch_;}

  /** */
   const ParametrizedCurve& paramCurve() const
     {return paramCurve_;}

private:

  /** */
  int id_;

  /** */
  RCP<CellFilterStub> domain_;

  /** */
  RCP<QuadratureFamilyStub> quad_;

  /** Such RQC might have one curve */
  ParametrizedCurve paramCurve_;

  /** */
  WatchFlag watch_;
          
  /** */
  static int getID(const RCP<CellFilterStub>& domain,
    const RCP<QuadratureFamilyStub>& quad,
    const WatchFlag& watch);

  /** */
  static int topID() {static int rtn=0; return rtn++;}

  /** */
  static Map<RegTriple, int>& domainAndQuadToIDMap() ;
};

}

namespace std
{
/** \relates Sundance::RegionQuadCombo*/
inline ostream& operator<<(std::ostream& os, 
  const Sundance::RegionQuadCombo& c)
{
  os << c.toString();
  return os;
}
}

namespace Teuchos
{
/** \relates Sundance::RegionQuadCombo */
inline std::string toString(const Sundance::RegionQuadCombo& h)
{return h.toString();}

}

#endif
