//******************************************************************************
/*!
  \file      src/Main/MeshConvDriver.C
  \author    J. Bakosi
  \date      Tue Apr 15 07:43:00 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MeshConvDriver that drives MeshConv
  \details   MeshConvDriver that drives MeshConv
*/
//******************************************************************************

#include <Factory.h>
#include <MeshConvDriver.h>
#include <MeshConv/CmdLine/Parser.h>
#include <GmshMeshReader.h>
#include <NetgenMeshReader.h>
#include <NetgenMeshWriter.h>
#include <GmshMeshWriter.h>
#include <ExodusIIMeshWriter.h>

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
  // Create 3D unstructured mesh
  quinoa::UnsMesh mesh;

  // Create mesh reader
  quinoa::GmshMeshReader
    inMesh( m_cmdline->get<tag::io, tag::input>(), mesh );
  //quinoa::NetgenMeshReader
  //  inMesh( m_cmdline->get<tag::io, tag::input>(), mesh );

  // Read in mesh
  inMesh.read();

  // Output mesh
  //quinoa::GmshMeshWriter outMesh( "out.msh", mesh );
  //quinoa::NetgenMeshWriter outMesh( "out.mesh", mesh );
  quinoa::ExodusIIMeshWriter outMesh( "out.exo", mesh );
  outMesh.write();
}
