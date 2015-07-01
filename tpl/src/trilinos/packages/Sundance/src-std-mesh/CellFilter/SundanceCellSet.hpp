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

#ifndef SUNDANCE_CELLSET_H
#define SUNDANCE_CELLSET_H

#include "SundanceDefs.hpp"
#include "SundanceCellSetBase.hpp"
#include "SundanceCellPredicate.hpp"
#include "PlayaHandle.hpp"


namespace Sundance
{
using namespace Teuchos;
  
/** 
 * CellSet is, you guessed it, a set of cells in a mesh. Cells are 
 * represented by their LID relative to the mesh. 
 * 
 * 
 * @see CellFilter, CellIterator
 **/
class CellSet : public Playa::Handle<CellSetBase>
{
public:
  /* handle boilerplate */
  HANDLE_CTORS(CellSet, CellSetBase);

  /** Construct from an explicit set of cells */
  CellSet(const Mesh& mesh, int cellDim,
    const CellType& cellType,
    const Set<int>& cellLIDs);
      

  /** The ID number of the mesh in which these cells exist */
  int meshID() const {return ptr()->meshID();}
      
  /** The mesh in which these cells exist */
  const Mesh& mesh() const {return ptr()->mesh();}

  /** Indicate whether the cells in this set are null cells */
  bool isNull() const {return ptr().get()==0 || ptr()->dimension() < 0;}

  /** The dimension of the cells contained in this set */
  int dimension() const {return ptr()->dimension();}

  /** The type of the cells contained in this set */
  const CellType& cellType() const {return ptr()->cellType();}

  /** An iterator pointing to the beginning of the set */
  CellIterator begin() const {return ptr()->begin();}

  /** An iterator pointing to the end of the set */
  CellIterator end() const {return ptr()->end();}

  /** Return a cell set that is the union of this set and another set */
  CellSet setUnion(const CellSet& other) const ;

  /** Return a cell set that is the intersection
   *  of this set and another set */
  CellSet setIntersection(const CellSet& other) const ;

  /** Return a cell set that is the difference
   *  of this set and another set */
  CellSet setDifference(const CellSet& other) const ;

  /** */
  CellSet subset(const RCP<CellPredicate>& test) const ;


  /** Determine whether all cells in this set are
   * facets of cells in the other set */
  bool areFacetsOf(const CellSet& other) const ;

  /** */
  bool operator<(const CellSet& other) const ;

  /** */
  int numCells() const ;

private:
  void checkCompatibility(const std::string& op, const CellSet& other) const ;
};


}

namespace std
{
inline ostream& operator<<(std::ostream& os, 
  const Sundance::CellSet& c)
{
  c.print(os);
  return os;
}
}



#endif
