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

#ifndef SUNDANCE_MESHREADERBASE_H
#define SUNDANCE_MESHREADERBASE_H


#include "SundanceDefs.hpp"
#include "SundanceMeshSourceBase.hpp"
#include "Teuchos_StrUtils.hpp"

namespace Sundance
{

/**
 * MeshReaderBase is a base class for mesh sources that get a mesh
 * from a file. It provides several utilities for parsing lines
 * from mesh files. 
 */
class MeshReaderBase : public MeshSourceBase
{
public:
  /** Construct with a filename */
  MeshReaderBase(const std::string& filename,
    const MeshType& meshType,
    int verbosity,
    const MPIComm& comm)
    : MeshSourceBase(meshType, verbosity, comm), filename_(filename)
    {}

  /** Construct from a parameter list */
  MeshReaderBase(const ParameterList& params);

  /** */
  virtual ~MeshReaderBase(){;}

protected:
  /** access to the filename */
  const std::string& filename() const {return filename_;}

  /** convert a std::string to its integer value */
  int atoi(const std::string& x) const ;

  /** convert a std::string to its double value */
  double atof(const std::string& x) const ;

  /** Determine whether a line is empty */
  bool isEmptyLine(const std::string& x) const ;

  /** Open a file "fname" and check for success.
   * @param fname name of the file to be opened
   * @param description a description of the file, e.g., "node file",
   * to be included in any error messages generated.  
   **/
  RCP<std::ifstream> openFile(const std::string& fname, 
    const std::string& description) const ;

  /** 
   * Read the next non-empty, non-comment line from a stream
   * @param is the stream from which to get the line
   * @param line upon return, filled in with the line that was read
   * @param tokens array of space-separated tokens in the line
   * @param comment a character indicating that everything after it
   * is a comment
   */
  bool getNextLine(std::istream& is, std::string& line,
    Array<string>& tokens,
    char comment) const ;
private:
  std::string filename_;
  mutable int nVertexVars_;
  mutable Array<double> vertexVars_;
};
}



#endif
