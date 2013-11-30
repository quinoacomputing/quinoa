//******************************************************************************
/*!
  \file      src/Main/Quinoa.C
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 05:01:58 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa main
  \details   Quinoa main
*/
//******************************************************************************

#include <Init.h>
#include <InitQuinoa.h>
#include <Config.h>
#include <QuinoaDriver.h>

int main(int argc, char* argv[])
//******************************************************************************
//  Quinoa: Lagrangian particle hydrodynamics
//! \author  J. Bakosi
//******************************************************************************
{
  return tk::Main< quinoa::QuinoaDriver, quinoa::echoTPL >
                 ( argc,
                   argv,
                   "Quinoa: Lagrangian particle hydrodynamics",
                   QUINOA_EXECUTABLE );
}
