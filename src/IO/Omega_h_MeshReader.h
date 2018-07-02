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

//#include "Omega_h_library.hpp"

#include "Types.h"
#include "Exception.h"

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

    //! Public interface to read our chunk of the mesh graph from Omega h file
    void readGraph( std::vector< std::size_t >& ginpoel, int n, int m );

    //! Read header from Omega h file
    std::size_t readHeader();

    //! Read coordinates of a number of mesh nodes from Omega h file
    std::array< std::vector< tk::real >, 3 >
    readCoords( const std::vector< std::size_t >& gid ) const;

    //! Read face list of all side sets from Omega h file
    std::size_t
    readSidesetFaces( std::map< int, std::vector< std::size_t > >& belem,
                      std::map< int, std::vector< int > >& faceid );

    //! Read face connectivity of a number boundary faces from Omega h file
    void readFaces( std::size_t nbfac, std::vector< std::size_t >& conn ) const;

    //! Read node list of all side sets from Omega h file
    std::map< int, std::vector< std::size_t > > readSidesets();

  private:
    const std::string m_filename;       //!< Input file name
};

} // tk::

#endif // Omega_h_MeshReader_h
