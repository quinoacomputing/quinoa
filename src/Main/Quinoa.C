//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Sat 02 Nov 2013 11:00:22 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <QuinoaDriver.h>

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Lagrangian particle hydrodynamics
//! \author  J. Bakosi
//******************************************************************************
{
  return tk::Main< quinoa::QuinoaDriver >( argc, argv,
           "Quinoa: Lagrangian particle hydrodynamics", QUINOA_EXECUTABLE);
}
