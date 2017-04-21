// *****************************************************************************
/*!
  \file      src/IO/GmshMeshIO.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Gmsh mesh reader and writer related types
  \details   Gmsh mesh reader and writer related types.
*/
// *****************************************************************************
#ifndef GmshMeshIO_h
#define GmshMeshIO_h

namespace tk {

//! Identifiers of supported Gmsh elements
enum GmshElemType { LIN = 1,
                    TRI = 2,
                    TET = 4,
                    PNT = 15 };

//! Supported Gmsh mesh file types
enum class GmshFileType { UNDEFINED = -1,
                          ASCII = 0,
                          BINARY = 1 };

} // tk::

#endif // GmshMeshIO_h
