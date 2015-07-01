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


#ifndef SUNDANCE_CELLPREDICATEBASE_H
#define SUNDANCE_CELLPREDICATEBASE_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "PlayaHandleable.hpp"
#include "Teuchos_XMLObject.hpp"
#include "PlayaPrintable.hpp"
#include "Teuchos_Describable.hpp"
#include "SundanceObjectWithVerbosity.hpp"
#include <typeinfo>

namespace Sundance
{
using namespace Teuchos;
  
/** 
 * CellPredicateBase is the base class for predicate objects
 * that test cells
 * against some condition. A simulation developer needing
 * some specialized method for identifying cells might implement
 * a custom cell predicate by extending this function. However,
 * the most common cases, selection by cell label or cell position,
 * have already been implemented 
 * in LabelCellPredicate and PositionalCellPredicate.
 */
class CellPredicateBase 
  : public Playa::Handleable<CellPredicateBase>,
    public Noncopyable,
    public ObjectWithClassVerbosity<CellPredicateBase>
{
public:
  /** Empty ctor */
  CellPredicateBase();

  /** virtual dtor */
  virtual ~CellPredicateBase();
      
  /** Test the predicate on a batch of cells */
  virtual void testBatch(const Array<int>& cellLID,
    Array<int>& results) const = 0 ;

      
  /** Set the current mesh and dimension on which cells are to be tested */
  virtual void setMesh(const Mesh& mesh, int cellDim) const 
    {mesh_ = mesh; cellDim_ = cellDim;}

  /** Write to XML */
  virtual XMLObject toXML() const = 0 ;

  /** */
  virtual bool lessThan(const CellPredicateBase* other) const = 0 ;

  /** */
  virtual std::string description() const = 0 ;


  /** */
  virtual std::string typeName() const {return typeid(*this).name();}
protected:

  /** */
  const Mesh& mesh() const {return mesh_;}


  /** */
  int cellDim() const {return cellDim_;}

private:
  mutable Mesh mesh_;

  mutable int cellDim_;



};

}

#endif
