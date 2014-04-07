//******************************************************************************
/*!
  \file      src/Main/Gmsh2ExoDriver.C
  \author    J. Bakosi
  \date      Mon 07 Apr 2014 07:15:02 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh2ExoDriver that drives Gmsh2Exo
  \details   Gmsh2ExoDriver that drives Gmsh2Exo
*/
//******************************************************************************

#include <make_unique.h>

#include <Factory.h>
#include <Gmsh2ExoDriver.h>
#include <Gmsh2Exo/CmdLine/Parser.h>
#include <GmshMeshReader.h>
#include <GmshMeshWriter.h>

using gmsh2exo::Gmsh2ExoDriver;

Gmsh2ExoDriver::Gmsh2ExoDriver(int argc, char** argv, const tk::Print& print)
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
Gmsh2ExoDriver::execute() const
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
