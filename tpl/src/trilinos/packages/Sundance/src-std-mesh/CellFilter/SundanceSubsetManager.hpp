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

#ifndef SUNDANCE_SUBSETMANAGER_H
#define SUNDANCE_SUBSETMANAGER_H

#include "SundanceDefs.hpp"
#include "SundanceMap.hpp"
#include "SundanceSet.hpp"
#include "SundanceCellFilter.hpp"
#include "SundanceCellFilterBase.hpp"

namespace Sundance
{
/** 
 * SubsetManager keeps maps from cell filters to their subsets. This is
 * a device to avoid circular references between subset and superset.
 */

class SubsetManager
{
public:
  /** */
  static void registerSubset(const CellFilter& filter, 
    const CellFilter& subset)
    {
      if (!subsetMap().containsKey(filter))
        subsetMap().put(filter, Sundance::Set<CellFilter>());
      subsetMap().get(filter).put(subset);
    }

  /** */
  static void registerDisjoint(const CellFilter& filter, 
    const CellFilter& subset)
    {
      if (!disjointMap().containsKey(filter))
        disjointMap().put(filter, Sundance::Set<CellFilter>());
      disjointMap().get(filter).put(subset);
    }

  /** */
  static void registerLabeledSubset(const CellFilter& filter, 
    const Set<int>& label, const CellFilter& subset)
    {
      if (!labeledMap().containsKey(filter))
        labeledMap().put(filter, Sundance::Map<Set<int>, CellFilter>());
      labeledMap().get(filter).put(label, subset);
    }


  /** */
  static const Sundance::Set<CellFilter>& getSubsets(const CellFilter& filter)
    {
      if (!subsetMap().containsKey(filter))
        subsetMap().put(filter, Sundance::Set<CellFilter>());
      return subsetMap().get(filter);
    }

  /** */
  static const Sundance::Set<CellFilter>& getDisjoints(const CellFilter& filter)
    {
      if (!disjointMap().containsKey(filter))
        disjointMap().put(filter, Sundance::Set<CellFilter>());
      return disjointMap().get(filter);
    }

  /** */
  static const Sundance::Map<Set<int>, CellFilter>& 
  getLabeledSubsets(const CellFilter& filter)
    {
      if (!labeledMap().containsKey(filter))
        labeledMap().put(filter, Sundance::Map<Set<int>, CellFilter>());
      return labeledMap().get(filter);
    }

private:
  
  /** */
  static Sundance::Map<CellFilter, Sundance::Set<CellFilter> >& subsetMap()
    {
      static Sundance::Map<CellFilter, Sundance::Set<CellFilter> > rtn;
      return rtn;
    }
  
  /** */
  static Sundance::Map<CellFilter, Sundance::Set<CellFilter> >& disjointMap()
    {
      static Sundance::Map<CellFilter, Sundance::Set<CellFilter> > rtn;
      return rtn;
    }
  
  /** */
  static Sundance::Map<CellFilter, Sundance::Map<Set<int>, CellFilter> >& labeledMap()
    {
      static Sundance::Map<CellFilter, Sundance::Map<Set<int>, CellFilter> > rtn;
      return rtn;
    }
};

}



#endif
