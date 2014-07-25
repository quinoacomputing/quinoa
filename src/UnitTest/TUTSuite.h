//******************************************************************************
/*!
  \file      src/UnitTest/TUTSuite.h
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 01:07:36 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Template Unit Test suite
  \details   Template Unit Test suite
*/
//******************************************************************************
#ifndef TUTSuite_h
#define TUTSuite_h

#include <tut/tut.hpp>

#include <UnitTestPrint.h>
#include <UnitTest/CmdLine/CmdLine.h>
#include <tutsuite.decl.h>

namespace unittest {

//! Template Unit Test unit test suite
class TUTSuite : public CBase_TUTSuite {

  public:
    //! Constructor
    explicit TUTSuite( const ctr::CmdLine& cmdline );

    //! Evaluate a unit test
    void evaluate();

  private:
    //! Custom Template Unit Test reporter
    struct Reporter : tut::callback {
      void test_completed( const tut::test_result& tr ) override {
        std::cout << tr.test << "=" << tr.result << std::endl; }
    };

    UnitTestPrint m_print;            //!< Pretty printer
    int m_maxTestsInGroup;            //!< Max. number of tests per group to run
    std::size_t m_ncomplete;          //!< Number of completed tests
};

} // unittest::

#endif // TUTSuite_h
