//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Thu 01 May 2014 10:02:29 PM MDT
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
  int numRecv;

  public:
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

class ArrayA : public CBase_ArrayA {
  public:
  ArrayA() {
    CkPrintf("ArrayA: created element %d\n", thisIndex);
  }
  ArrayA(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy)); }
};

class ArrayB : public CBase_ArrayB {
  public:
  ArrayB() {
    CkPrintf("ArrayB: created element (%d,%d)\n", thisIndex.x, thisIndex.y);
  }
  ArrayB(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy)); }
};

class ArrayC : public CBase_ArrayC {
  public:
  ArrayC() {
    CkPrintf("ArrayB: created element (%d,%d,%d)\n", thisIndex.x, thisIndex.y, thisIndex.z);
  }
  ArrayC(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finished), mainProxy)); }
};

class ArrayD : public CBase_ArrayD {
  public:
  ArrayD() {}
  ArrayD(CkMigrateMessage*) { }
  void e() { }
};

class ArrayE : public CBase_ArrayE {
  public:
  ArrayE() {}
  ArrayE(CkMigrateMessage*) { }
  void e() { }
};

class ArrayF : public CBase_ArrayF {
  public:
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
