//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Wed Mar 19 13:19:41 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <QuinoaDriver.h>
#include <TPLInfo/Silo.h>
#include <TPLInfo/HDF5.h>
#include <TPLInfo/Zlib.h>

namespace quinoa {

void echoTPL(const tk::Print& print)
//******************************************************************************
//  Echo TPL version informaion for libs specific to Quinoa
//! \author  J. Bakosi
//******************************************************************************
{
  //print.raw("\n");
  //echoZoltan(print, "Zoltan library");
  print.raw("\n");
  tk::echoSilo(print, "Silo library");
  print.raw("\n");
  tk::echoHDF5(print, "HDF5 library");
  print.raw("\n");
  tk::echoZlib(print, "Zlib compression library");
}

} // quinoa::

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Lagrangian particle hydrodynamics
//! \author  J. Bakosi
//******************************************************************************
{
  return tk::Main< quinoa::QuinoaDriver >
                 ( argc,
                   argv,
                   "Quinoa: Lagrangian particle hydrodynamics",
                   QUINOA_EXECUTABLE,
                   quinoa::echoTPL );
}
