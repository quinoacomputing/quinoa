// *****************************************************************************
/*!
  \file      src/Main/UnitTestPrint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     UnitTest's printer
  \details   UnitTest's printer
*/
// *****************************************************************************
#ifndef UnitTestPrint_h
#define UnitTestPrint_h

#include <sstream>

#include "Types.hpp"
#include "Print.hpp"
#include "Exception.hpp"

namespace unittest {

//! UnitTestPrint : tk::Print
class UnitTestPrint : public tk::Print {

  public:
    //! Constructor
    //! \param[in] screen Screen output filename
    //! \param[in,out] str Verbose stream
    //! \param[in] mode Open mode for screen output file, see
    //!   http://en.cppreference.com/w/cpp/io/ios_base/openmode
    //! \param[in,out] qstr Quiet stream
    //! \see tk::Print::Print
    explicit UnitTestPrint( const std::string& screen = {},
                            std::ostream& str = std::clog,
                            std::ios_base::openmode mode = std::ios_base::out,
                            std::ostream& qstr = std::cout ) :
      Print( screen, str, mode, qstr ) {}

    //! Print unit tests header (with legend)
    //! \param[in] t Section title
    //! \param[in] group String attempting to match unit test groups
    void unithead( const std::string& t, const std::string& group ) const {
      std::string g = group.empty() ? "all" : group;
      m_stream << m_section_title_fmt % m_section_indent
                                      % m_section_bullet
                                      % t;
      m_stream << m_section_underline_fmt
                  % m_section_indent
                  % std::string( m_section_indent.size() + 2 + t.size(),
                                '-' );
      raw( m_item_indent + "Groups: " + g + " (use -g str to match groups)\n" +
           m_item_indent + "Legend: [done/failed] group:test : result\n\n" );
    }

    //! Print one-liner info for test
    //! \details Columns:
    //!   [done/failed]
    //!   - done: number of tests completed so far
    //!   - failed: number of failed tests so far
    //!   name of the test group
    //!   name of the test
    //!   result (with additional info if failed)
    //!   Assumed fields for status:
    //!   - status[0]: test group name
    //!   - status[1]: test name
    //!   - status[2]: result (tut::test_result::result_type as string)
    //!   - status[3]: exception message for failed test
    //!   - status[4]: exception type id for failed test
    void test( std::size_t ncomplete,
               std::size_t nfail,
               const std::vector< std::string >& status )
    {
      if (status[2] != "8") {             // if not dummy
        std::stringstream ss;
        ss << "[" << ncomplete << "/" << nfail << "] " << status[0] << ":"
           << status[1];
        m_stream <<
          m_item_widename_value_fmt % m_item_indent % ss.str()
                                    % result( status[2], status[3], status[4] )
          << std::flush;
      }
    }

  private:
    //! Return human-readable test result based on result code
    //! \param[in] code Result code
    //! \param[in] msg Message to append
    //! \param[in] ex Expection message to attach to exceptions cases
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
