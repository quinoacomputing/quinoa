//******************************************************************************
/*!
  \file      src/IO/GmshMeshIO.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 08:06:18 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Gmsh mesh reader and writer related types
  \details   Gmsh mesh reader and writer related types.
*/
//******************************************************************************
#ifndef GmshMeshIO_h
#define GmshMeshIO_h

namespace quinoa {

//! Identifiers of supported Gmsh elements
enum GmshElemType { LIN = 1,
                    TRI = 2,
                    TET = 4,
                    PNT = 15 };

//! Supported Gmsh mesh file types
enum class GmshFileType { UNDEFINED = -1,
                          ASCII = 0,
                          BINARY = 1 };

} // quinoa::

#endif // GmshMeshIO_h
