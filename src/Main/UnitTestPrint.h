//******************************************************************************
/*!
  \file      src/Main/UnitTestPrint.h
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 04:54:27 PM MDT
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
#include <Exception.h>

namespace unittest {

//! UnitTestPrint : Print
class UnitTestPrint : public tk::Print {

  public:
    //! Constructor
    explicit UnitTestPrint( std::ostream& str = tk::null,
                            std::ostream& qstr = std::cout ) :
      Print( str, qstr ), m_ncomplete(0), m_nfail(0) {}

    //! Print unit tests header (with legend)
    void unithead( const std::string& title, std::size_t ngroup ) const {
      std::stringstream ss;
      ss << title << " (from " << ngroup << " test groups)";
      m_qstream << m_section_title_fmt % m_section_indent
                                       % m_section_bullet
                                       % ss.str();
      m_qstream << m_section_underline_fmt
                   % m_section_indent
                   % std::string( m_section_indent.size() + 2 + ss.str().size(),
                                 '-');
      raw< tk::QUIET >
         ( m_item_indent + "Legend: [done/failed] group:test : result\n\n" );
    }

    //! Print one-liner info for test. Columns:
    //! [done/failed]
    //!   - done: number of tests completed so far
    //!   - failed: number of failed tests so far
    //! name of the test group
    //! name of the test
    //! result (with additional info if failed)
    //! Assumed fields for status:
    //!   status[0]: test group name
    //!   status[1]: test name
    //!   status[2]: result (tut::test_result::result_type as string)
    //!   status[3]: exception message for failed test
    //!   status[4]: exception type id for failed test
    void test( const std::vector< std::string >& status ) {
      if (status[2] != "8") {             // if not dummy
        if (status[2] != "0") ++m_nfail;  // count number of failed tests
        std::stringstream ss;
        ss << "[" << ++m_ncomplete << "/" << m_nfail << "] " << status[0] << ":"
           << status[1];
        (status[2] == "0" ? m_stream : m_qstream) <<
          m_item_widename_value_fmt % m_item_indent % ss.str()
                                    % result( status[2], status[3], status[4] );
      }
    }

  private:
    std::size_t m_ncomplete;            //!< Count number of completed tests
    std::size_t m_nfail;                //!< Count number of failed tests

    //! Return human-readable test result based on result code
    std::string result( const std::string& code,
                        const std::string& msg,
                        const std::string& ex ) const
    {
      if (code == "0") return "ok";
      else if (code == "1") return "fail: " + msg;
      else if (code == "2") return "except: " + msg + ex;
      else if (code == "3") return "warning: " + msg;
      else if (code == "4") return "terminate: " + msg;
      else if (code == "5") return "ex_ctor: " + msg + ex;
      else if (code == "6") return "rethrown: " + msg + ex;
      else if (code == "7") return "skipped: " + msg;
      else if (code == "8") return "dummy";
      else Throw( "No such unit test result code found" );
    }
};

} // unittest::

#endif // UnitTestPrint_h
