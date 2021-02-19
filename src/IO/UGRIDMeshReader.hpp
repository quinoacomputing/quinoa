// *****************************************************************************
/*!
  \file      src/IO/UGRIDMeshReader.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     UGRID mesh reader class declaration
  \details   UGRID mesh reader class declaration. Mesh reader facilitating
             reading a mesh from a simple text file used by NASA.
  \see       http://www.simcenter.msstate.edu/software/downloads/doc/ug_io/3d_grid_file_type_ugrid.html, http://www.simcenter.msstate.edu/software/downloads/doc/aflr3/aflr3_io_summary.html
*/
// *****************************************************************************
#ifndef UGRIDMeshReader_h
#define UGRIDMeshReader_h

#include <iosfwd>

#include "Reader.hpp"

namespace tk {

class UnsMesh;

//! \brief UGRIDMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a simple text
//!    file used by NASA.
class UGRIDMeshReader : public Reader {

  public:
    //! Constructor
    explicit UGRIDMeshReader( const std::string& filename ) :
      Reader( filename ), m_nnode(0), m_ntet(0), m_ntri(0) {}

    //! Read UGRID mesh
    void readMesh( UnsMesh& mesh );

  private:
    std::size_t m_nnode;        //!< Number of nodes
    std::size_t m_ntet;         //!< Number of tetrahedra
    std::size_t m_ntri;         //!< Number of triangles

    //! Read header
    void readHeader();

    //! Read nodes
    void readNodes( UnsMesh& mesh );

    //! Read element connectivity
    void readElements( UnsMesh& mesh );
};

} // tk::

#endif // UGRIDMeshReader_h
