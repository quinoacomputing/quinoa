// *****************************************************************************
/*!
  \file      src/IO/MeshWriter.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ nodegroup for outputing mesh data to file
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

#include "Types.h"
#include "Options/FieldFile.h"
#include "Centering.h"
#include "UnsMesh.h"

#include "NoWarning/meshwriter.decl.h"

namespace tk {

//! Charm++ nodegroup used to output particle data to file in parallel
class MeshWriter : public CBase_MeshWriter {

  public:
    //! Constructor: set some defaults that stay constant at all times
    MeshWriter( const std::string& output_basefilename,
                ctr::FieldFileType filetype,
                Centering bnd_centering,
                bool benchmark );

    //! Set the total number of chares
    void expect( int nchare );

    //! Output unstructured mesh into file
    void writeMesh( uint64_t itr,
                    int chareid,
                    const std::vector< std::size_t >& inpoel,
                    const UnsMesh::Coords& coord,
                    const std::map< int, std::vector< std::size_t > >& bface,
                    const std::vector< std::size_t >& triinpoel,
                    const std::map< int, std::vector< std::size_t > >& bnode,
                    const std::unordered_map< std::size_t, std::size_t >& lid );

    //! Output metadata (field names) to file
    void writeMeta( uint64_t itr,
                    int chareid,
                    Centering centering,
                    const std::vector< std::string >& names );

    //! Output field data to file
    void writeFields( uint64_t itr,
                      uint64_t itf,
                      tk::real time,
                      int chareid,
                      Centering centering,
                      const std::vector< std::vector< tk::real > >& fields );

  private:
    //! String to use as the base of the filename
    const std::string m_outputBasefilename;
    //! Output file format type
    const ctr::FieldFileType m_filetype;
    //! Centering to identify what boundary data to write.
    const Centering m_bndCentering;
    //! True if benchmark mode
    const bool m_benchmark;

    //! Total number chares across the whole problem
    int m_nchare;

    //! Compute filename
    std::string filename( uint64_t itr, int chareid ) const;
};

} // tk::

#endif // MeshWriter_h
