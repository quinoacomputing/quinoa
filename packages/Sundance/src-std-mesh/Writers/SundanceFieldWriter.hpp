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

#ifndef SUNDANCE_FIELDWRITER_H
#define SUNDANCE_FIELDWRITER_H

#include "SundanceDefs.hpp"
#include "SundanceFieldWriterBase.hpp"
#include "PlayaHandle.hpp"

namespace Sundance
{
  /**
   * FieldWriter is the user level object for writing fields and meshes
   * to output stream. 
   *
   * <h4> Example: </h4> Write fields u0 and w0 to a VTK file "results.vtu"
   * \code
   * FieldWriter vtkWriter = new VTKWriter("results");
   * vtkWriter.addField(u0);
   * vtkWriter.addField(w0);
   * vtkWriter.write();
   * \endcode
   *
   * <h4> Example: </h4> Write verbose mesh information to cout
   * \code
   * FieldWriter writer = new VerboseFieldWriter();
   * writer.addMesh(mesh);
   * writer.write();
   * \endcode
   */
  class FieldWriter : public Playa::Handle<FieldWriterBase>
  {
  public:
    /* Boilerplate handle ctors */
    HANDLE_CTORS(FieldWriter, FieldWriterBase);

    /** add a mesh to the list of things to be written */
    void addMesh(const Mesh& mesh) const ;

    /** add a field, tagging it with the given std::string as a name */
    void addField(const std::string& name, 
      const Playa::Handle<FieldBase>& field) ;

    /** set the numerical value to be written at cells on which
     * a field is undefined. Default value is 0.0. */
    void setUndefinedValue(const double& x);

    /** 
     * Tell the writer to pretend that it is running as one of nProc processes
     * with the given rank. This is used when partitioning meshes in serial.
     */
    void impersonateParallelProc(int nProc, int rank)
      {
        ptr()->impersonateParallelProc(nProc, rank);
      }
    
    

    /** write to stream */
    void write() const ;
  };
}

#endif
