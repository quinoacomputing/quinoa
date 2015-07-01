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

#ifndef SUNDANCE_TRIANGLEWRITER_H
#define SUNDANCE_TRIANGLEWRITER_H


#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"

namespace Sundance
{
/**
 * TriangleWriter writes a mesh or fields to a file
 * in Shewchuk's Triangle format.
 */
class TriangleWriter : public FieldWriterBase
{
public:
  /** */
  TriangleWriter(const std::string& filename="", int indexOffset=0) 
    : FieldWriterBase(filename), indexOffset_(indexOffset) {;}
    
  /** virtual dtor */
  virtual ~TriangleWriter(){;}

  /** */
  virtual void write() const ;

  /** Return a ref count pointer to self */
  virtual RCP<FieldWriterBase> getRcp() {return rcp(this);}

protected:
  /** */
  void writePoints(const std::string& filename) const ;

  /** */
  void writeCells(const std::string& filename) const ;
    
  /** */
  void writeEdges(const std::string& filename) const ;
    
  /** */
  void writeFaces(const std::string& filename) const ;
    
  /** */
  void writeHeader(const std::string& filename) const ;
    
  /** */
  void writeParallelInfo(const std::string& filename) const ;
    
  int indexOffset_;
};

/** 
 * TriangleWriterFactory produces an Triangle writer in contexts where a user cannot
 * do so directly.
 */
class TriangleWriterFactory : public FieldWriterFactoryBase
{
public:
  /** */
  TriangleWriterFactory() {}

  /** Create a writer with the specified filename */
  RCP<FieldWriterBase> createWriter(const string& name) const 
    {return rcp(new TriangleWriter(name));}

  /** */
  virtual RCP<FieldWriterFactoryBase> getRcp() {return rcp(this);}
  
};



}




#endif
