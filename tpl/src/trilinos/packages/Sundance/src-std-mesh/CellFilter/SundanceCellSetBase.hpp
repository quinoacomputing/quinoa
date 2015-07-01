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

#ifndef SUNDANCE_CELLSETBASE_H
#define SUNDANCE_CELLSETBASE_H

#include "SundanceDefs.hpp"
#include "SundanceSet.hpp"
#include "SundanceMap.hpp"
#include "SundanceCellType.hpp"
#include "SundanceCellIterator.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "SundanceMesh.hpp"
#include "SundanceNoncopyable.hpp"
#include "PlayaHandleable.hpp"
#include "PlayaHandle.hpp"

namespace Sundance
{
using namespace Teuchos;

/** 
 * CellSetBase is the base class for cell sets. There are two cell
 * set subtypes: ExplicitCellSet and ImplicitCellSet.
 *
 * @see CellFilter
 **/
class CellSetBase : public ObjectWithClassVerbosity<CellSetBase>,
                    public Playa::Printable,
                    public Noncopyable,
                    public Playa::Handleable<CellSetBase>
{
public:
  /** Construct, initializing to an empty set */
  CellSetBase(const Mesh& mesh, int cellDim,
    const CellType& cellType);

  /** Return an iterator pointing to the first element in the set */
  virtual CellIterator begin() const = 0 ;

  /** Return an iterator containing the past-the-end value */
  virtual CellIterator end() const = 0 ;

  /** Return the type of cells in this set */
  const CellType& cellType(const CellType& cellType) const
    {return cellType_;}

  /** The ID number of the mesh in which these cells exist */
  int meshID() const {return mesh_.id();}

  /** The dimension of the cells contained in this set */
  int dimension() const {return dim_;}
      
  /** The mesh in which these cells exist */
  const Mesh& mesh() const {return mesh_;}

  /** The type of the cells contained in this set */
  const CellType& cellType() const {return cellType_;}

  /** */
  bool lessThan(const CellSetBase* other) const ;

  /** */
  virtual bool internalLessThan(const CellSetBase* other) const = 0 ;

private:

  /** the mesh in which the set exists */
  Mesh mesh_;

  /** the type of cell in the set */
  CellType cellType_;

  /** the dimension of the cells in the set */
  int dim_;
      
};
}



#endif
