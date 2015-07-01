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

#ifndef SUNDANCE_EXODUSWRITER_H
#define SUNDANCE_EXODUSWRITER_H


#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"
#include "SundanceCellType.hpp"

namespace Sundance
{
/**
 * ExodusWriter writes a mesh or fields to an ExodusII file
 */
class ExodusWriter : public FieldWriterBase
{
public:
  /** */
  ExodusWriter(const std::string& filename) 
    : FieldWriterBase(filename) {;}
    
  /** virtual dtor */
  virtual ~ExodusWriter(){;}

  /** */
  virtual void write() const ;

  /** Return a ref count pointer to self */
  virtual RCP<FieldWriterBase> getRcp() {return rcp(this);}

  /** */
  void writeParallelInfo(const std::string& filename) const ;


private:    
  /** */
  void getCharpp(const Array<std::string>& s, Array<const char*>& p) const ;

  /** */
  void findNodeSets(
    Array<CellFilter>& nodesetFilters,
    Array<int>& omnipresentFuncs,
    Array<RCP<Array<int> > >& funcsForNodeset,
    Array<RCP<Array<int> > >& nodesForNodeset,
    Array<int>& nsID,
    Array<int>& nNodesPerSet,
    Array<int>& nsNodePtr,
    RCP<Array<int> > allNodes
    ) const ;

  /** */
  void findBlocks(
    Array<CellFilter>& blockFilters,
    Array<int>& omnipresentFuncs,
    Array<RCP<Array<int> > >& funcsForBlock,
    Array<RCP<Array<int> > >& elemsForBlock,
    Array<int>& elemIDs,
    Array<int>& nElemsPerBlock,
    Array<int>& blockElemPtr,
    RCP<Array<int> > allElems
    ) const ;



  /** */
  void offset(Array<int>& x) const ;

  /** */
  std::string elemType(const CellType& type) const ;

  /** */
  void writeMesh(int exoID, 
    const Array<CellFilter>& nodesetFilters,
    const Array<int>& nsID,
    const Array<int>& nNodesPerSet,
    const Array<int>& nsNodePtr,
    const RCP<Array<int> >& allNodes) const ;

  /** */
  void writeFields(int exoID, 
    const Array<CellFilter>& nodesetFilters,
    const Array<int>& omnipresentNodalFuncs,
    const Array<int>& omnipresentElemFuncs,
    const Array<RCP<Array<int> > >& funcsForNodeset,
    const Array<RCP<Array<int> > >& nodesForNodeset,
    const Array<int>& nsID) const ;
    
    
};


/** 
 * ExodusWriterFactory produces an Exodus writer in contexts where a user cannot
 * do so directly.
 */
class ExodusWriterFactory : public FieldWriterFactoryBase
{
public:
  /** */
  ExodusWriterFactory() {}

  /** Create a writer with the specified filename */
  RCP<FieldWriterBase> createWriter(const string& name) const 
    {return rcp(new ExodusWriter(name));}

  /** */
  virtual RCP<FieldWriterFactoryBase> getRcp() {return rcp(this);}
  
};

}


#endif
