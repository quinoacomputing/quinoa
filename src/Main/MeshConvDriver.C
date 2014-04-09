//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Tue 08 Apr 2014 09:29:32 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshConvDriver that drives MeshConv
  \details   MeshConvDriver that drives MeshConv
*/
//******************************************************************************

#include <make_unique.h>

#include <Factory.h>
#include <MeshConvDriver.h>
#include <MeshConv/CmdLine/Parser.h>
#include <GmshMeshReader.h>
#include <GmshMeshWriter.h>

using meshconv::MeshConvDriver;

MeshConvDriver::MeshConvDriver(int argc, char** argv, const tk::Print& print)
//******************************************************************************
//  Constructor
//! \param[in] argc      Argument count from command line
//! \param[in] argv      Argument vector from command line
//! \param[in] print     Simple pretty printer
//! \author J. Bakosi
//******************************************************************************
{
  // Parse command line into cmdline
  CmdLineParser cmdParser(argc, argv, print, m_cmdline);
}

void
MeshConvDriver::execute() const
//******************************************************************************
//  Execute: Convert Gmsh mesh to Exodus II mesh
//! \author J. Bakosi
//******************************************************************************
{
  // Create Gmsh mesh
  quinoa::GmshMesh mesh;

  // Create Gmsh mesh reader
  quinoa::GmshMeshReader
    inMesh( m_cmdline->get<tag::io, tag::input>(), mesh );

  // Read in Gmsh mesh
  inMesh.read();

  // Debug: echo element sets
  //mesh.echoElemSets();

  // Debug: output Gmsh mesh
  quinoa::GmshMeshWriter outMesh( "out.msh", mesh );
  outMesh.write();
}
