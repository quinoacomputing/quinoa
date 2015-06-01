//******************************************************************************
/*!
  \file      src/Main/RegTestPrint.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:40:10 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     RegTest's pretty printer
  \details   RegTest's pretty printer
*/
//******************************************************************************
#ifndef RegTestPrint_h
#define RegTestPrint_h

#include <sstream>

#include "Types.h"
#include "Print.h"
#include "Exception.h"

namespace regtest {

//! RegTestPrint : tk::Print
class RegTestPrint : public tk::Print {

  public:
    //! Constructor
    //! \param[inout] str Verbose stream
    //! \param[inout] qstr Quiet stream
    //! \see tk::Print::Print
    //! \author J. Bakosi
    explicit RegTestPrint( std::ostream& str = std::clog,
                           std::ostream& qstr = std::cout ) :
      Print( str, qstr ) {}

    //! Print regression tests header (with legend)
    //! \param[in] title Section title
    //! \param[in] group String attempting to match regression test groups
    //! \author J. Bakosi
    void reghead( const std::string& title, const std::string& group ) const {
      std::string g = group.empty() ? "all" : group;
      m_stream << m_section_title_fmt % m_section_indent
                                      % m_section_bullet
                                      % title;
      m_stream << m_section_underline_fmt
                  % m_section_indent
                  % std::string( m_section_indent.size() + 2 + title.size(),
                                '-' );
      raw( m_item_indent + "Groups: " + g + " (use -g str to match groups)\n" +
           m_item_indent + "Legend: [done/failed] group:test : result\n\n" );
    }

    //! \brief Print one-liner info for test.
    //! \details Columns:
    //!   [done/failed]
    //!   - done: number of tests completed so far
    //!   - failed: number of failed tests so far
    //!   name of the test
    //!   result (with additional info if failed)
    //!   Assumed fields for status:
    //!   - status[0]: test name
    //!   - status[1...]: test output
    //! \author J. Bakosi
    void test( std::size_t ncomplete,
               std::size_t nfail,
               const std::vector< std::string >& status )
    {
      std::cout << "B--------------\n";
      for (auto& line : status) std::cout << line << std::endl;
      std::cout << "E--------------\n";
//       std::stringstream ss;
//       ss << "[" << ncomplete << "/" << nfail << "] " << status[0] << ":"
//          << status[1];
//       (status[2] == "0" ? m_stream : m_qstream) <<
//         m_item_widename_value_fmt % m_item_indent % ss.str()
//                                   % result( status[2], status[3], status[4] );
    }

  private:
    //! Return human-readable test result based on result code
    //! \param[in] code Result code
    //! \param[in] msg Message to append
    //! \param[in] ex Expection message to attach to exceptions cases
    //! \author J. Bakosi
    std::string result( const std::string& code,
                        const std::string& msg,
                        const std::string& ex ) const
    {
      return "ok";
//       if (code == "0") return "ok";
//       else if (code == "1") return "fail: " + msg;
//       else if (code == "2") return "except: " + msg + ex;
//       else if (code == "3") return "warning: " + msg;
//       else if (code == "4") return "terminate: " + msg;
//       else if (code == "5") return "ex_ctor: " + msg + ex;
//       else if (code == "6") return "rethrown: " + msg + ex;
//       else if (code == "7") return "skipped: " + msg;
//       else if (code == "8") return "dummy";
//       else Throw( "No such unit test result code found" );
    }
};

} // regtest::

#endif // RegTestPrint_h
