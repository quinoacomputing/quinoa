//******************************************************************************
/*!
  \file      src/Main/UnitTestPrint.h
  \author    J. Bakosi
  \date      Thu 24 Jul 2014 11:35:13 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     UnitTest's printer
  \details   UnitTest's printer
*/
//******************************************************************************
#ifndef UnitTestPrint_h
#define UnitTestPrint_h

#include <sstream>

#include <Types.h>
#include <Print.h>
#include <flip_map.h>

namespace unittest {

//! UnitTestPrint : Print
class UnitTestPrint : public tk::Print {

  public:
    //! Constructor
    explicit UnitTestPrint( std::ostream& str = tk::null,
                            std::ostream& qstr = std::cout ) :
      Print( str, qstr ) {}

    //! Print one-liner info for test. Columns:
    //! [done/total/failed]
    //!   - done: number of tests completed so far
    //!   - total: total number of tests
    //!   - failed: number of failed tests
    //! name of the test group
    //! name of the test
    //! result of test: "pass" or "fail"
    void test( std::size_t ncomplete, std::size_t ntest ) const
    {
      std::stringstream ss;
      ss << "[" << ncomplete << "/" << ntest << "/" << "nfailed" << "] ";
      stream() <<
        m_item_widename_value_fmt % m_item_indent % ss.str() % "status";
    }

};

} // unittest::

#endif // UnitTestPrint_h
