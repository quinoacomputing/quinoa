// *****************************************************************************
/*!
  \file      src/IO/Omega_h_MeshReader.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Omega_h mesh reader
  \details   Omega_h mesh reader class declaration.
*/
// *****************************************************************************
#ifndef Omega_h_MeshReader_h
#define Omega_h_MeshReader_h

#include <map>
#include <vector>
#include <array>
#include <unordered_map>

#include "Types.h"
#include "Exception.h"
#include "UnsMesh.h"

namespace tk {

//! Omega_h mesh-based data reader
//! \details Mesh reader class facilitating reading from mesh-based field data
//!   a file in Omega_h format.
//! \see https://github.com/ibaned/omega_h
class Omega_h_MeshReader {

  public:
    //! Constructor
    explicit Omega_h_MeshReader( const std::string& filename )
      : m_filename( filename ) {}

    //! Read part of the mesh (graph and coords) from Omega_h file
    //! \details Total number of PEs defaults to 1 for a single-CPU read, this
    //!    PE defaults to 0 for a single-CPU read.
    void readMeshPart( std::vector< std::size_t >& ginpoel,
                       std::vector< std::size_t >& inpoel,
                       std::vector< std::size_t >& triinp,
                       std::vector< std::size_t >& gid,
                       std::unordered_map< std::size_t, std::size_t >& lid,
                       tk::UnsMesh::Coords& coord,
                       int numpes=1, int mype=0 );

    //! Read face list of all side sets from Omega h file
    void
    readSidesetFaces( std::map< int, std::vector< std::size_t > >& bface,
                      std::map< int, std::vector< std::size_t > >& faces );

    //! Read face connectivity of a number boundary faces from Omega h file
    void readFaces( std::vector< std::size_t >& conn ) const;

    //! Read node list of all side sets from Omega h file
    std::map< int, std::vector< std::size_t > > readSidesetNodes();

    //! ...
    std::vector< std::size_t > triinpoel(
      std::map< int, std::vector< std::size_t > >& belem,
      const std::map< int, std::vector< std::size_t > >& faces,
      const std::vector< std::size_t >& ginpoel,
      const std::vector< std::size_t >& triinp ) const;

  private:
    const std::string m_filename;       //!< Input file name
};

} // tk::

#endif // Omega_h_MeshReader_h
