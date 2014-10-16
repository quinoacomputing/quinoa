//******************************************************************************
/*!
  \file      src/IO/GmshMeshIO.h
  \author    J. Bakosi
  \date      Thu 10 Apr 2014 09:30:12 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Gmsh mesh reader and writer common
  \details   Gmsh mesh reader and writer common
*/
//******************************************************************************
#ifndef GmshMeshIO_h
#define GmshMeshIO_h

namespace quinoa {

// Identifiers of supported Gmsh elements
enum GmshElemType { LIN = 1,
                    TRI = 2,
                    TET = 4,
                    PNT = 15 };

// Gmsh mesh file types
enum class GmshFileType { UNDEFINED = -1,
                          ASCII = 0,
                          BINARY = 1 };

} // quinoa::

#endif // GmshMeshIO_h
