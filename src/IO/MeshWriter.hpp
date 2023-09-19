// *****************************************************************************
/*!
  \file      src/IO/MeshWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
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
                bool benchmark,
                std::size_t nmesh );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit MeshWriter( CkMigrateMessage* m ) : CBase_MeshWriter( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Set the total number of chares
    void nchare( std::size_t meshid, int n );

    //! Output unstructured mesh into file
    void write( std::size_t meshid,
                bool meshoutput,
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
                const std::vector< std::string >& elemsurfnames,
                const std::vector< std::string >& nodesurfnames,
                const std::vector< std::vector< tk::real > >& elemfields,
                const std::vector< std::vector< tk::real > >& nodefields,
                const std::vector< std::vector< tk::real > >& elemsurfs,
                const std::vector< std::vector< tk::real > >& nodesurfs,
                const std::set< int >& outsets,
                CkCallback c );

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \note This is a Charm++ group, pup() is thus only for
    //!    checkpoint/restart.
    void pup( PUP::er &p ) override {
      p | m_filetype;
      p | m_bndCentering;
      p | m_benchmark;
      p | m_nmesh;
      p | m_nchare;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] m MeshWriter object reference
    friend void operator|( PUP::er& p, MeshWriter& m ) { m.pup(p); }
    //@}

  private:
    //! Output file format type
    ctr::FieldFileType m_filetype;
    //! Centering to identify what boundary data to write.
    Centering m_bndCentering;
    //! True if benchmark mode
    bool m_benchmark;
    //! Total number of meshes
    std::size_t m_nmesh;
    //! Total number chares across the whole problem (one per mesh)
    std::vector< int > m_nchare;

    //! Compute filename
    std::string filename( const std::string& basefilename,
                          std::size_t meshid,
                          uint64_t itr,
                          int chareid,
                          int surfid = 0 ) const;
};

} // tk::

#endif // MeshWriter_h
