//******************************************************************************
/*!
  \file      src/Main/Gmsh2ExoDriver.C
  \author    J. Bakosi
  \date      Wed Mar 19 15:47:23 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh2ExoDriver that drives Gmsh2Exo
  \details   Gmsh2ExoDriver that drives Gmsh2Exo
*/
//******************************************************************************

#include <make_unique.h>

#include <Factory.h>
#include <Gmsh2ExoDriver.h>
//#include <Gmsh2Exo/CmdLine/Parser.h>

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
  IGNORE(argc);
  IGNORE(argv);
  IGNORE(print);
}

void
Gmsh2ExoDriver::execute() const
//******************************************************************************
//  Execute
//! \author J. Bakosi
//******************************************************************************
{
}
