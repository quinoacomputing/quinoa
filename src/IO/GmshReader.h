//******************************************************************************
/*!
  \file      src/IO/GmshReader.h
  \author    J. Bakosi
  \date      Mon 10 Sep 2012 04:22:38 AM KST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh reader class declaration
  \details   Gmsh reader class declaration
*/
//******************************************************************************
#ifndef GmshReader_h
#define GmshReader_h

#include <UnsMesh.h>
#include <MeshReader.h>
#include <GmshException.h>

namespace Quinoa {

//! GmshReader : MeshReader
class GmshReader : MeshReader {

  public:
    //! Constructor
    GmshReader(string filename);

    //! Destructor
    ~GmshReader();

    //! Interface for read
    virtual void read(UnsMesh* mesh);

  private:
    //! Don't permit copy operator
    GmshReader(const GmshReader&);
    //! Don't permit assigment operator
    GmshReader& operator=(const GmshReader&);

    //! Read mandatory "$MeshFormat--$EndMeshFormat" section
    void readMeshFormat();

    //! Read "$Nodes--$EndNodes" section
    void readNodes();

    //! Read "$Elements--$EndElements" section
    void readElements();

    //! Read "$PhysicalNames--$EndPhysicalNames" section
    void readPhysicalNames();
};

} // namespace Quinoa

#endif // GmshReader_h
