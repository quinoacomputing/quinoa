// *****************************************************************************
/*!
  \file      src/IO/MeshWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ group for outputing mesh data to file
  \details   Charm++ group declaration used to output data associated to
     unstructured meshes to file(s). Charm++ chares (work units) send mesh and
     field data associated to mesh entities to the MeshWriter class defined here
     to write the data to file(s).
*/
// *****************************************************************************
#ifndef MeshWriter_h
#define MeshWriter_h

#include <vector>
#include <string>
#include <tuple>
#include <map>

#include "Types.hpp"
#include "Options/FieldFile.hpp"
#include "Centering.hpp"
#include "UnsMesh.hpp"

#include "NoWarning/meshwriter.decl.h"

namespace tk {

//! Charm++ group used to output particle data to file in parallel
class MeshWriter : public CBase_MeshWriter {

  public:
    //! Constructor: set some defaults that stay constant at all times
    MeshWriter( ctr::FieldFileType filetype,
                Centering bnd_centering,
                bool benchmark );

    //! Set the total number of chares
    void nchare( int n );

    //! Output unstructured mesh into file
    void write( bool meshoutput,
                bool fieldoutput,
                uint64_t itr,
                uint64_t itf,
                tk::real time,
                int chareid,
                const std::string& basefilename,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                const std::map< int, std::vector< std::size_t > >& bface,
                const std::map< int, std::vector< std::size_t > >& bnode,
                const std::vector< std::size_t >& triinpoel,
                const std::vector< std::string >& elemfieldnames,
                const std::vector< std::string >& nodefieldnames,
                const std::vector< std::vector< tk::real > >& elemfields,
                const std::vector< std::vector< tk::real > >& nodefields,
                CkCallback c );

  private:
    //! Output file format type
    const ctr::FieldFileType m_filetype;
    //! Centering to identify what boundary data to write.
    const Centering m_bndCentering;
    //! True if benchmark mode
    const bool m_benchmark;

    //! Total number chares across the whole problem
    int m_nchare;

    //! Compute filename
    std::string filename( const std::string& basefilename,
                          uint64_t itr,
                          int chareid ) const;
};

} // tk::

#endif // MeshWriter_h
