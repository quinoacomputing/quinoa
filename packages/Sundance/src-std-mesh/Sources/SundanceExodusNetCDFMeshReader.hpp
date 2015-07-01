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

#ifndef SUNDANCE_EXODUSNETCDFMESHREADER_H
#define SUNDANCE_EXODUSNETCDFMESHREADER_H

#include "SundanceDefs.hpp"
#include "SundanceMeshReaderBase.hpp"
#include "SundanceMap.hpp"
#include "Teuchos_Array.hpp"

namespace Sundance
{
  using namespace Teuchos;
  
  /**
   * ExodusNetCDFMeshReader reads a mesh from a NetCDF dump of an Exodus file.
   * This will often be less efficient than reading from an exodus file
   * directly, but does not require any proprietary libraries.
   * 
   * Utilities to dump exodus to NetCDF are available from 
   * <A HREF="http://my.unidata.ucar.edu/content/software/netcdf/index.html"> 
   * </A>
   */
  class ExodusNetCDFMeshReader : public MeshReaderBase
  {
  public:
    /** */
    ExodusNetCDFMeshReader(const std::string& filename, 
      const MeshType& meshType,
      int verbosity=0,
      const MPIComm& comm = MPIComm::world());
    /** */
    ExodusNetCDFMeshReader(const ParameterList& params);

    /** virtual dtor */
    virtual ~ExodusNetCDFMeshReader(){;}


    /** Create a mesh */
    virtual Mesh fillMesh() const ;

    /** Print a short descriptive std::string */
    virtual std::string description() const 
    {return "ExodusNetCDFMeshReader[file=" + filename() + "]";}


    /** Return a ref count pointer to self */
    virtual RCP<MeshSourceBase> getRcp() {return rcp(this);}

  private:

    
                      
  };
}

#endif
