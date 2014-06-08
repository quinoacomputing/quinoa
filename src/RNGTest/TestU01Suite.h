//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:21:28 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 random number generator test suite
  \details   TestU01 random number generator test suite
*/
//******************************************************************************
#ifndef TestU01Suite_h
#define TestU01Suite_h

#include <StatTest.h>

namespace rngtest {

extern std::vector< tk::ctr::RNGType > g_selectedrng;

//! TestU01 random number generator test suite used polymorphically with Battery
//! Suite is a policy specific to a TestU01 test suite
template< class Suite >
class TestU01Suite {

  public:
    //! Constructor
    explicit TestU01Suite()
    { for (const auto& s : g_selectedrng) Suite().addTests( m_tests, s ); }

    //! Run battery of RNG tests
    void run() const {}

    //! Print list of registered statistical tests
    void print( const RNGTestPrint& pr ) const
    { pr.names< StatTest >( m_tests, ntest() ); }

    //! Evaluate a statistical test
    void evaluate( std::size_t id ) const {}

    //! Return number of statistical tests
    std::size_t ntest() const {
      return !g_selectedrng.empty() ? m_tests.size() / g_selectedrng.size() : 0;
    }

    //! Return number of statistics produced
    std::size_t nstat() const {
      std::size_t nstats = 0;
      for (const auto& t : m_tests) nstats += t.nstat();
      return !g_selectedrng.empty() ? nstats / g_selectedrng.size() : 0;
    }

  private:
    std::vector< StatTest > m_tests;    //!< Statistical tests
};

} // rngtest::

#endif // TestU01Suite_h
