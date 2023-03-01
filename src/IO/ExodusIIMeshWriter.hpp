// *****************************************************************************
/*!
  \file      src/IO/ExodusIIMeshWriter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     ExodusII mesh-based data writer
  \details   ExodusII mesh-based data writer class declaration.
*/
// *****************************************************************************
#ifndef ExodusIIMeshWriter_h
#define ExodusIIMeshWriter_h

#include <cstddef>
#include <iosfwd>
#include <vector>
#include <map>

#include "Types.hpp"
#include "UnsMesh.hpp"

namespace tk {

class UnsMesh;

//! ExodusII writer constructor modes
enum class ExoWriter { CREATE, OPEN };

//! ExodusII mesh-based data writer
//! \details Mesh writer class facilitating writing a mesh and associated
//!   mesh-based field data to a file in ExodusII format.
//! \see http://sourceforge.net/projects/exodusii
class ExodusIIMeshWriter {

  public:
    //! Constructor: create/open ExodusII file
    explicit ExodusIIMeshWriter( const std::string& filename,
                                 ExoWriter mode,
                                 int cpuwordsize = sizeof(double),
                                 int iowordsize = sizeof(double) );

    //! Destructor
    ~ExodusIIMeshWriter() noexcept;

    //! Write ExodusII mesh file taking a tk::UnsMesh object
    void writeMesh( const UnsMesh& mesh ) const;

    //! Write ExodusII mesh file taking inputs to a tk::UnsMesh object
    void writeMesh( const std::vector< std::size_t >& tetinp,
                    const UnsMesh::Coords& coord,
                    const std::map< int, std::vector< std::size_t > >& bface,
                    const std::vector< std::size_t >& triinp ) const;

    //! Write ExodusII mesh file taking inputs to a tk::UnsMesh object
    void writeMesh( const std::vector< std::size_t >& tetinp,
                    const UnsMesh::Coords& coord,
                    const std::map< int, std::vector< std::size_t > >& bnode )
      const;

    //! Write ExodusII mesh file taking inputs to a tk::UnsMesh object
    //! \tparam nnode 3 or 4, indicating a triangle or tetrahedron mesh
    //! \param[in] inpoel Element connectivity
    //! \param[in] coord Node coordinates
    template< std::size_t nnode >
    void writeMesh( const std::vector< std::size_t >& inpoel,
                    const UnsMesh::Coords& coord ) const
    {
      if (nnode == 4)
        writeMesh( tk::UnsMesh( inpoel, coord ) );
      else if (nnode == 3)
        writeMesh( tk::UnsMesh( coord, inpoel ) );
    }

    //!  Write time stamp to ExodusII file
    void writeTimeStamp( uint64_t it, tk::real time ) const;

    //! Write time values to ExodusII file
    void writeTimeValues( const std::vector< tk::real >& tv ) const;

    //! Write the names of nodal output variables to ExodusII file
    void writeNodeVarNames( const std::vector< std::string >& nv ) const;

    //! Write the names of element output variables to ExodusII file
    void writeElemVarNames( const std::vector< std::string >& ev ) const;

    //! \brief Write multiple node scalar fields to ExodusII file at multiple
    //!   time steps
    void writeNodeScalars(
      const std::vector< std::vector< std::vector< tk::real > > >& var ) const;

    //!  Write node scalar field to ExodusII file
    void writeNodeScalar( uint64_t it,
                          int varid,
                          const std::vector< tk::real >& var ) const;

    //! \brief Write multiple element scalar fields to ExodusII file at multiple
    //!   time steps
    void writeElemScalars(
      const std::vector< std::vector< std::vector< tk::real > > >& var ) const;

    //!  Write elem scalar field to ExodusII file
    void writeElemScalar( uint64_t it,
                          int varid,
                          const std::vector< tk::real >& var ) const;

    //! Write header without mesh, function overloading
    void writeHeader( const char* title, int64_t ndim, int64_t nnodes,
                      int64_t nelem, int64_t nblk, int64_t node_set,
                      int64_t side_set ) const;

    //! Write nodes without mesh, function overloading.
    void writeNodes( const std::vector< tk::real >& x, 
                     const std::vector< tk::real >& y,
                     const std::vector< tk::real >& z ) const;

    //! Write element block to ExodusII file
    void writeElemBlock( int& elclass,
                         int64_t nnpe,
                         const std::string& eltype,
                         const std::vector< std::size_t >& inpoel ) const;

  private:
    //! Write ExodusII header
    void writeHeader( const UnsMesh& mesh ) const;

    //! Write nodes coordinates to ExodusII file
    void writeNodes( const UnsMesh& mesh ) const;

    //! Write element conectivity to ExodusII file
    void writeElements( const UnsMesh& mesh ) const;

    //! Write side sets and their face connectivity to ExodusII file
    void writeSidesets( const UnsMesh& mesh ) const;

    //! Write side sets and their node list to ExodusII file
    void writeNodesets( const UnsMesh& mesh ) const;

    const std::string m_filename;          //!< File name
    int m_outFile;                         //!< ExodusII file handle
};

} // tk::

#endif // ExodusIIMeshWriter_h
