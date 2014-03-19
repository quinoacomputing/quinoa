//******************************************************************************
/*!
  \file      src/Main/InitQuinoa.C
  \author    J. Bakosi
  \date      Wed Mar 19 11:48:53 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa-specific initialization for main
  \details   Quinoa-specific initialization for main
*/
//******************************************************************************

#include <InitQuinoa.h>
#include <TPLInfo/Silo.h>
#include <TPLInfo/HDF5.h>
#include <TPLInfo/Zlib.h>

void quinoa::echoTPL(const tk::Print& print)
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
