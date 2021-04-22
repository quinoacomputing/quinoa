// *****************************************************************************
/*!
  \file      src/IO/RDGFLOMeshReader.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     RDGFLO mesh reader class declaration
  \details   RDGFLO mesh reader class declaration. Mesh reader facilitating
             reading a mesh from a simple text file used by Prof. Hong Luo at
             North Carolina State University.
*/
// *****************************************************************************
#ifndef RDGFLOMeshReader_h
#define RDGFLOMeshReader_h

#include <iosfwd>

#include "Reader.hpp"

namespace tk {

class UnsMesh;

//! \brief RDGFLOMeshReader : tk::Reader
//! \details Mesh reader class facilitating reading a mesh from a simple text
//!    file used by Prof. Hong Luo at North Carolina State University.
class RDGFLOMeshReader : public Reader {

  public:
    //! Constructor
    explicit RDGFLOMeshReader( const std::string& filename ) :
      Reader( filename ), m_nnode(0), m_ntet(0), m_ntri(0) {}

    //! Read RDGFLO mesh
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

#endif // RDGFLOMeshReader_h
