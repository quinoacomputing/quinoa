//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 10:10:25 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
#include <TPLInfo/MKL.h>
#include <TPLInfo/Boost.h>
#include <TPLInfo/OpenMP.h>

namespace quinoa {

void echoTPL(const tk::Print& print)
//******************************************************************************
//  Echo TPL version informaion for libs specific to Quinoa
//! \author  J. Bakosi
//******************************************************************************
{
  //echoZoltan(print, "Zoltan library");
  echoOpenMP( print, "OpenMP runtime" );
#ifdef HAS_MKL
  echoMKL( print, "Intel Math Kernel Library" );
#else
  print.item( "Intel Math Kernel Library", "n/a" );
#endif
  echoBoost( print, "Boost C++ Libraries" );
  tk::echoSilo(print, "Silo library");
  tk::echoHDF5(print, "HDF5 library");
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
