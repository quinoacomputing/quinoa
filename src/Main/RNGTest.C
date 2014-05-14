//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Wed 14 May 2014 12:27:11 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa random number generator test suite
  \details   Quinoa random number generator test suite
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <RNGTestDriver.h>
#include <rngtest.decl.h>

namespace rngtest {

void echoTPL(const tk::Print& /*print*/)
//******************************************************************************
//  Echo TPL version informaion for libs specific to RNGTest
//! \author  J. Bakosi
//******************************************************************************
{
}

} // rngtest::

/*readonly*/ CProxy_Main mainProxy;

class Main : public CBase_Main {

  public:

    // Cosntructor
    Main( CkArgMsg* msg ) {
      tk::Main< rngtest::RNGTestDriver >
              ( msg->argc,
                msg->argv,
                "Quinoa: Random number generator (RNG) test suite",
                RNGTEST_EXECUTABLE,
                rngtest::echoTPL );
      delete msg;
      //mainProxy = thisProxy;
      //CkExit();
    }

    // Destructor
    //~Main() { CkExit(); }
};

#include <rngtest.def.h>
