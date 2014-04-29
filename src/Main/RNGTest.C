//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Tue Apr 29 10:54:38 2014
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

struct Main : CBase_Main {
  int numRecv;

  Main( CkArgMsg* msg ) : numRecv(0)
  {
    delete msg;

    mainProxy = thisProxy;

    const int X = 3;
    const int Y = 3;
    const int Z = 3;
    CProxy_ArrayA a1 = CProxy_ArrayA::ckNew(X);
    CProxy_ArrayB b1 = CProxy_ArrayB::ckNew(X, Y);
    CProxy_ArrayC c1 = CProxy_ArrayC::ckNew(X, Y, Z);

    a1.e();
    b1.e();
    c1.e();
  }

  void finished() {
    if (++numRecv == 3) CkExit();
  }
};

struct ArrayA : CBase_ArrayA {
  ArrayA() {
    CkPrintf("ArrayA: created element %d\n", thisIndex);
  }
  ArrayA(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy)); }
};

struct ArrayB : CBase_ArrayB {
  ArrayB() {
    CkPrintf("ArrayB: created element (%d,%d)\n", thisIndex.x, thisIndex.y);
  }
  ArrayB(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy)); }
};

struct ArrayC : CBase_ArrayC {
  ArrayC() {
    CkPrintf("ArrayB: created element (%d,%d,%d)\n", thisIndex.x, thisIndex.y, thisIndex.z);
  }
  ArrayC(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy)); }
};

struct ArrayD : CBase_ArrayD {
  ArrayD() {}
  ArrayD(CkMigrateMessage*) { }
  void e() { }
};

struct ArrayE : CBase_ArrayE {
  ArrayE() {}
  ArrayE(CkMigrateMessage*) { }
  void e() { }
};

struct ArrayF : CBase_ArrayF {
  ArrayF() {}
  ArrayF(CkMigrateMessage*) { }
  void e() { }
};

#include <rngtest.def.h>

// int main(int argc, char* argv[])
// //******************************************************************************
// //  Quinoa: Random number generator test suite
// //! \author  J. Bakosi
// //******************************************************************************
// {
//   return tk::Main< rngtest::RNGTestDriver >
//                  ( argc,
//                    argv,
//                    "Quinoa: Random number generator (RNG) test suite",
//                    RNGTEST_EXECUTABLE,
//                    rngtest::echoTPL );
// }
