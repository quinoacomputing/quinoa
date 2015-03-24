//******************************************************************************
/*!
  \file      src/RegTest/Regression.h
  \author    J. Bakosi
  \date      Tue 24 Mar 2015 08:33:21 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Regression test class declaration
  \details   Regression test class declaration.
*/
//******************************************************************************
#ifndef Regression_h
#define Regression_h

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <pstreams/pstream.h>

#include <regression.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace regtest {

//! Charm++ chare class for regression tests
template< class Proxy >
class Regression : public CBase_Regression< Proxy > {

  public:
    //! Constructor: run test then call back to host proxy to evaluate it
    //! \param[in] proxy Host proxy to call back to after test has been run
    //! \param[in] testname Name of the test
    explicit Regression( Proxy& proxy, const std::string& testcmd ) {

      redi::ipstream is( testcmd );

      std::vector< std::string > output { testcmd };
      std::string line;
      while ( std::getline( is, line ) ) output.push_back( line );

      proxy.evaluate( output );
    }
};

} // regtest::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include <regression.def.h>
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // Regression_h
