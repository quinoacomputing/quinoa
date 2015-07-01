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

#ifndef SUNDANCE_FIELDWRITERBASE_H
#define SUNDANCE_FIELDWRITERBASE_H


#include "SundanceDefs.hpp"
#include "SundanceMesh.hpp"
#include "SundanceFieldBase.hpp"

namespace Sundance
{
using namespace Teuchos;
/**
 * FieldWriterBase is a base class for objects that write fields
 * and/or meshes to a stream. 
 */
class FieldWriterBase : public Playa::Handleable<FieldWriterBase>,
                        public ObjectWithClassVerbosity<FieldWriterBase>
{
public:
  /** */
  FieldWriterBase(const std::string& filename);

  /** virtual dtor */
  virtual ~FieldWriterBase(){;}

  /** */
  void addMesh(const Mesh& mesh);

  /** add a comment */
  virtual void addCommentLine(const std::string& line) ;

  /** add a field, tagging it with the given std::string as a name */
  virtual void addField(const std::string& name, 
    const RCP<FieldBase>& field) ;

  /** */
  virtual void write() const = 0 ;

  /**  */
  virtual void impersonateParallelProc(int nProc, int rank) ;

  /** set the numerical value to be written at cells on which
   * a field is undefined. */
  void setUndefinedValue(const double& x) {undefinedValue_ = x;}

protected:

  /** */
  double undefinedValue() const {return undefinedValue_;}
  /** */
  int nProc() const ;

  /** */
  int myRank() const ;

  /** */
  const std::string& filename() const {return filename_;}

  /** */
  const Mesh& mesh() const {return mesh_;}

  /** Indicate whether the given writer subtype does anything special
   * for vector fields. Default is false, in which case
   * vectors are simply written as a list of scalars.
   */
  virtual bool supportsSpecializedVectors() const {return false;}

  const Array<string>& comments() const {return comments_;}
  Array<string>& comments() {return comments_;}

  const Array<RCP<FieldBase> >& pointScalarFields() const {return pointScalarFields_;}
  Array<RCP<FieldBase> >& pointScalarFields() {return pointScalarFields_;}

  const Array<RCP<FieldBase> >& cellScalarFields() const {return cellScalarFields_;}
  Array<RCP<FieldBase> >& cellScalarFields() {return cellScalarFields_;}

  const Array<string>& pointScalarNames() const {return pointScalarNames_;}
  Array<string>& pointScalarNames() {return pointScalarNames_;}

  const Array<string>& cellScalarNames() const {return cellScalarNames_;}
  Array<string>& cellScalarNames() {return cellScalarNames_;}

  const Array<RCP<FieldBase> >& pointVectorFields() const {return pointVectorFields_;}
  Array<RCP<FieldBase> >& pointVectorFields() {return pointVectorFields_;}

  const Array<RCP<FieldBase> >& cellVectorFields() const {return cellVectorFields_;}
  Array<RCP<FieldBase> >& cellVectorFields() {return cellVectorFields_;}

  const Array<string>& pointVectorNames() const {return pointVectorNames_;}
  Array<string>& pointVectorNames() {return pointVectorNames_;}

  const Array<string>& cellVectorNames() const {return cellVectorNames_;}
  Array<string>& cellVectorNames() {return cellVectorNames_;}

  virtual void writeCommentLine(const std::string& line) const {;}

private:
  std::string filename_;

  Mesh mesh_;

  int nProc_;

  int myRank_;

  int meshID_;

  Array<string> comments_;

  Array<RCP<FieldBase> > pointScalarFields_;
  Array<RCP<FieldBase> > cellScalarFields_;
  Array<RCP<FieldBase> > pointVectorFields_;
  Array<RCP<FieldBase> > cellVectorFields_;
  Array<string> pointScalarNames_;
  Array<string> cellScalarNames_;
  Array<string> pointVectorNames_;
  Array<string> cellVectorNames_;

  double undefinedValue_;
};

/**
 * FieldWriterFactoryBase
 */
class FieldWriterFactoryBase : public Playa::Handleable<FieldWriterFactoryBase>
{
public:
  /** Create a writer with the specified filename */
  virtual RCP<FieldWriterBase> createWriter(const string& name) const = 0 ;
}; 

}



#endif
