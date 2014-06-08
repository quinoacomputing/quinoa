//******************************************************************************
/*!
  \file      src/Main/RNGTest.C
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 01:49:09 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa random number generator test suite
  \details   Quinoa random number generator test suite
*/
//******************************************************************************

#include <Init.h>
#include <Config.h>
#include <Paradigm.h>
#include <RNG.h>
#include <RNGStack.h>
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

// Global-scope data. Initialized by the main chare and distibuted to all PEs by
// the Charm++ runtime system. Though semantically not const, all these global
// data should be considered read-only. See also http://charm.cs.illinois.edu/
// manuals/html/charm++/manual.html. The data below is global-scope because they
// must be available to all PEs which could be on different machines. Of course,
// in a previous non-Charm++ design, most of this data was held at class-level,
// but since the generators in g_rng must be possible to be called from
// global-scope, as external generators to TestU01, it is easier to make
// g_rng global-scope, as well the additional data required to initialize it.
// All other data below is required for serializing during migration of g_rng
// across the network. Note that the container (std::map) holding tk::RNG
// objects uses value semantics which is safer and generally less error-prone
// than reference semantics. At the same time tk::RNG is used in a polymorphic
// fashion with various classes that adhere to the concepts required by Concept
// defined inside tk::RNG. tk::RNG does not define a default, i.e.,
// non-templated constructor, since then the "derived" class object could not be
// initialized rendering the class tk::RNG empty-constructed, which invites
// abuse and ill-defined behavior. As such, the "derived" class type comes
// through the constructors and thus would not be available for a pack/unpack
// migrator required by Charm++ from within. Templating the class tk::RNG is not
// an option since then we could not hold tk::RNG objects in a simple
// std::vector. As a result of the above requirements, tk::RNG is migrated (here
// in global-scope) by reinstantiating RNGStack, which reinstatiates the RNG
// factory, from which the RNGs selected by the user are instantiated. Also,
// RNGFactory associates tk::ctr::RNG ids (enum class values) to function
// pointers (std::function objects pointing to tk::RNG constructors bound with
// their arguments). Since function pointers cannot simply be serialized and
// migrated via the network, they must also be recreated on remote machines.
// This initial migration of global-scope data is done by the Charm++ runtime
// once the main chare constructor is finished -- see the RNGTestDriver
// constructor, which initializes the data required for the migration).

CProxy_Main mainProxy;

#ifdef HAS_MKL
tk::ctr::RNGMKLParameters g_rngmklparam;
#endif

tk::ctr::RNGSSEParameters g_rngsseparam;

std::vector< tk::ctr::RNGType > g_selectedrng;

//! Pack/Unpack selected RNGs
inline void operator|( PUP::er& p, std::vector< tk::ctr::RNGType >& s )
{ tk::pup_vector( p, s ); }

std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

//! PUP std::vector< tk::RNG >
inline void operator|( PUP::er& p,
  std::map< tk::ctr::RawRNGType, tk::RNG >& rng )
{
  tk::RNGFactory factory;
  tk::RNGStack stack;
  stack.initFactory( factory, tk::Paradigm().ompNthreads(),
                     #ifdef HAS_MKL
                     g_rngmklparam,
                     #endif
                     g_rngsseparam );
  rng = stack.createSelected( factory, g_selectedrng );
}

//! Charm++ main chare
class Main : public CBase_Main {

  public:
    // Constructor
    Main( CkArgMsg* msg ) :
      m_driver( tk::Main< RNGTestDriver >( msg->argc, msg->argv,
                "Quinoa: Random number generator (RNG) test suite",
                 RNGTEST_EXECUTABLE,
                 echoTPL ) )
    {
      delete msg;
      mainProxy = thisProxy;
      // Fire up an asynchronous execute object, which when created at some
      // future point in time will call back to this->execute(). This is
      // necessary so that this->execute() can access already migrated
      // global-scope data.
      CProxy_execute::ckNew();
    }

    void execute() { m_driver.execute(); }
    void finalize() { CkExit(); }

  private:
    RNGTestDriver m_driver;
};

//! Charm++ chare execute: by the time this object is constructed, the Charm++
//! runtime system has finished migrating all global-scoped readonly objects
//! which happens after the main chare constructor has finished.
class execute : public CBase_execute {
 public:
   execute() {
     mainProxy.execute();
     mainProxy.finalize();
   }
};

} // rngtest::

#include <rngtest.def.h>
