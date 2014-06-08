//******************************************************************************
/*!
  \file      src/Main/MeshConv.C
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 04:00:11 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh to Exodus II mesh file converter
  \details   Gmsh to Exodus II mesh file converter
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <MeshConvDriver.h>
#include <TPLInfo/ExodusII.h>

namespace meshconv {

void echoTPL(const tk::Print& print)
//******************************************************************************
//  Echo TPL version informaion for libs specific to MeshConv
//! \author  J. Bakosi
//******************************************************************************
{
  print.raw("\n");
  echoExodusII(print, "ExodusII library");
}

} // meshconv::

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Gmsh to Exodus II mesh file converter
//! \author  J. Bakosi
//******************************************************************************
{
  meshconv::MeshConvDriver driver =
    tk::Main< meshconv::MeshConvDriver >
            ( argc, argv,
              "Quinoa: Mesh file converter",
              MESHCONV_EXECUTABLE,
              meshconv::echoTPL );
  driver.execute();
}
