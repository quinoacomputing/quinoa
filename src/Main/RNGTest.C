//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Thu 05 Jun 2014 07:35:38 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa random number generator test suite
  \details   Quinoa random number generator test suite
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <RNGTestDriver.h>
#include <rngtest.decl.h>
#include <PUPUtil.h>

namespace rngtest {

void echoTPL(const tk::Print& /*print*/)
//******************************************************************************
//  Echo TPL version informaion for libs specific to RNGTest
//! \author  J. Bakosi
//******************************************************************************
{
}

// Global-scope data: initialized by the main chare and distibuted to all PEs by
// the Charm++ runtime system. Though semantically not const, all these global
// data are considered read-only. See also
// http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.

CProxy_Main mainProxy;

#ifdef HAS_MKL
tk::ctr::RNGMKLParameters g_rngmklparam;
#endif

tk::ctr::RNGSSEParameters g_rngsseparam;

std::vector< tk::ctr::RNGType > g_selectedrng;

//! Pack/Unpack selected RNGs
inline void operator|( PUP::er& p, std::vector< tk::ctr::RNGType >& s )
{ tk::pup_vector( p, s ); }

std::map< tk::ctr::RNGType, tk::RNG > g_rng;

//! PUP std::vector< tk::RNG >
inline
void operator|( PUP::er& p, std::map< tk::ctr::RNGType, tk::RNG >& rng ) {
  tk::RNGFactory rngfactory;
  tk::RNGDriver rngdriver;
  rngdriver.initFactory( rngfactory, tk::Paradigm().ompNthreads(),
                         #ifdef HAS_MKL
                         g_rngmklparam,
                         #endif
                         g_rngsseparam );
  rng = rngdriver.createSelected( rngfactory, g_selectedrng );
}

class Main : public CBase_Main {

  public:

    // Cosntructor
    Main( CkArgMsg* msg ) : numRecv(0) {
      tk::Main< RNGTestDriver >
              ( msg->argc,
                msg->argv,
                "Quinoa: Random number generator (RNG) test suite",
                RNGTEST_EXECUTABLE,
                echoTPL );
      delete msg;
      mainProxy = thisProxy;

//      CProxy_test t = CProxy_test::ckNew(0);
//      CProxy_test q = CProxy_test::ckNew(1);

      CProxy_ArrayA a1 = CProxy_ArrayA::ckNew(3);
      CProxy_ArrayB b1 = CProxy_ArrayB::ckNew(3, 3);
      CProxy_ArrayC c1 = CProxy_ArrayC::ckNew(3, 3, 3);
      a1.e();
      b1.e();
      c1.e();
    }

    void finalize() {
      if (++numRecv == 3) CkExit();
    }

  private:
    int numRecv;
};

//class test : public CBase_test {
//  public:
//    test() {
//      CkPrintf("Hello, my PE is %d\n", CkMyPe());
//      mainProxy.finalize();
//    }
//};

struct ArrayA : CBase_ArrayA {
  ArrayA() {
    //CkPrintf("ArrayA: created element %d\n", thisIndex);
  }
  ArrayA(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finalize), mainProxy)); }
};
struct ArrayB : CBase_ArrayB {
  ArrayB() {
    //CkPrintf("ArrayB: created element (%d,%d)\n", thisIndex.x, thisIndex.y);
  }
  ArrayB(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finalize), mainProxy)); }
};

struct ArrayC : CBase_ArrayC {
  ArrayC() {
    //CkPrintf("ArrayB: created element (%d,%d,%d)\n", thisIndex.x, thisIndex.y, thisIndex.z);
  }
  ArrayC(CkMigrateMessage*) { }
  void e() { contribute(CkCallback(CkReductionTarget(Main, finalize), mainProxy)); }
};

} // rngtest::

#include <rngtest.def.h>
