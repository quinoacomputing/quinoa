// *****************************************************************************
/*!
  \file      src/IO/MeshReader.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Mesh reader class for inciter connecting to various readers
  \details   Mesh reader ckass for inciter connecting to various lower level
             mesh readers.
*/
// *****************************************************************************
#ifndef MeshReader_h
#define MeshReader_h

#include <vector>
#include <array>

#include "ExodusIIMeshReader.h"

namespace inciter {

//! Mesh reader class for inciter connecting to various readers
class MeshReader {

  public:
    //! Constructor
    //! \param[in] filename File to read mesh from
    //! \param[in] npes Total number of PEs
    //! \param[in] mype This PE
    explicit MeshReader( const std::string& filename, int npes, int mype ) :
      m_er( filename ),
      m_npes( static_cast< std::size_t >( npes ) ),
      m_mype( static_cast< std::size_t >( mype ) ) {}

    //! Read our chunk of the mesh graph (connectivity) from file
    void readGraph( std::vector< long >& gelemid,
                    std::vector< std::size_t >& ginpoel );

    //! Read coordinates of a number of mesh nodes from ExodusII file
    std::array< std::vector< tk::real >, 3 >
    readCoords( const std::vector< std::size_t > gid );

  private:
    tk::ExodusIIMeshReader m_er;        //!< ExodusII mesh reader object
    std::size_t m_npes;                 //!< Total number of PEs
    std::size_t m_mype;                 //!< This PE
};

} // inciter::

#endif // MeshReader_h
