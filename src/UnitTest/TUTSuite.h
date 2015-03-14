//******************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.h
  \author    J. Bakosi
  \date      Thu 12 Mar 2015 10:11:31 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Template Unit Test suite class declaration
  \details   Template Unit Test suite class declaration. In principle there can
    be unit test suites other than this one which uses the Template Unit Test
    library.
*/
//******************************************************************************
#ifndef TUTSuite_h
#define TUTSuite_h

#include <tut/tut.hpp>

#include <UnitTestPrint.h>
#include <UnitTest/CmdLine/CmdLine.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <tutsuite.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace unittest {

//! Template Unit Test unit test suite
class TUTSuite : public CBase_TUTSuite {

  public:
    //! Constructor
    explicit TUTSuite( const ctr::CmdLine& cmdline );

    //! Evaluate a unit test
    void evaluate( std::vector< std::string > status );

  private:
    UnitTestPrint m_print;      //!< Pretty printer
    std::size_t m_nmpi;         //!< Number of MPI test groups
    std::size_t m_nrun;         //!< Number of tests ran (including dummies)
    std::size_t m_ngroup;       //!< Number of test groups
    std::size_t m_ncomplete;    //!< Number of completed tests
    std::size_t m_nfail;        //!< Number of failed tests
    std::size_t m_nskip;        //!< Number of skipped tests
    std::size_t m_nwarn;        //!< Number of tests with a warning
    std::size_t m_nexcp;        //!< Number of tests with an exception
};

} // unittest::

#endif // TUTSuite_h
