//******************************************************************************
/*!
  \file      src/RegTest/RegSuite.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 10:13:46 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Regression test suite class declaration
  \details   Regression test suite class declaration.
*/
//******************************************************************************
#ifndef RegSuite_h
#define RegSuite_h

#include <vector>
#include <iosfwd>
#include <cstddef>

#include "RegTestPrint.h"
#include "RegTest/CmdLine/CmdLine.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "regsuite.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace tk { class Print; }

namespace regtest {

//! Regression test suite
class RegSuite : public CBase_RegSuite {

  public:
    //! Constructor
    explicit RegSuite( const ctr::CmdLine& cmdline );

    //! Evaluate a regression test
    void evaluate( std::vector< std::string > status );

  private:
    RegTestPrint m_print;       //!< Pretty printer
    std::size_t m_ntest;        //!< Number of tests to run
    std::size_t m_nrun;         //!< Number of tests ran
    std::size_t m_ncomplete;    //!< Number of completed tests
    std::size_t m_nfail;        //!< Number of failed tests
    std::size_t m_nskip;        //!< Number of skipped tests
    std::size_t m_nwarn;        //!< Number of tests with a warning
    std::size_t m_nexcp;        //!< Number of tests with an exception

    //! Fire up a test
    void spawntest( const std::string& t );

    //! Evaluate a single regression test
    void evaluate( std::vector< std::string > status,
                   std::size_t& ncomplete,
                   std::size_t& nwarn,
                   std::size_t& nskip,
                   std::size_t& nexcp,
                   std::size_t& nfail );

    //! Echo final assessment after the full regression test suite has finished
    void assess( const tk::Print& print,
                 std::size_t nfail,
                 std::size_t nwarn,
                 std::size_t nskip,
                 std::size_t nexcp,
                 std::size_t ncomplete );

};

} // regtest::

#endif // RegSuite_h
