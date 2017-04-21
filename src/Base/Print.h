// *****************************************************************************
/*!
  \file      src/Base/Print.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     General purpose pretty printer functionality
  \details   This file contains general purpose printer functions. Using the
    functions defined here provides formatting, and a consistent look with
    simple client-side code. For formatting, the Boost Format library is used,
    see http://www.boost.org/doc/libs/release/libs/format.
*/
// *****************************************************************************
#ifndef Print_h
#define Print_h

#include <iostream>
#include <sstream>
#include <iomanip>
#include <list>
#include <cmath>
#include <array>

#include "NoWarning/format.h"

#include "Timer.h"
#include "Exception.h"
#include "Has.h"

namespace tk {

//! Output verbosity. C-style enum as this is used for template argument.
enum Style { QUIET=0, VERBOSE=1 };

//! Pretty printer base. Contains general purpose printer functions. Using the
//! functions defined here provides formatting, and a consistent look with
//! simple client-side code. For formatting, the Boost Format library is used,
//! see http://www.boost.org/doc/libs/release/libs.
class Print {

  public:
    //! Constructor: Quiet output by default, only stuff written to qstr shown.
    //! Instantiate with str = std::cout for verbose output. Any member function
    //! can be called by overriding the default stream via the template
    //! argument, Style, a C-style enum. Note: By default, str == std::clog.
    //! This is used to initialize str to a local stringstream into which all
    //! verbose output goes by default, i.e., it will not be shown. This
    //! solution is chosen instead of trickery with null-streams, as
    //! boost:formatted output into null-streams caused invalid reads in
    //! valgrind. This way quiet output (formatted or not) simply goes into a
    //! local stringstream. In other words, the default argument to str,
    //! std::clog, is only used to detect whether client code passed a default
    //! argument or not: if it did not, the string stream is used for verbose
    //! output, if it did, the specified stream is used for the verbose output.
    //! \param[in,out] str Verbose stream
    //! \param[in,out] qstr Quiet stream
    //! \author J. Bakosi
    explicit Print( std::ostream& str = std::clog,
                    std::ostream& qstr = std::cout ) :
      m_null(),
      m_stream( str.rdbuf() == std::clog.rdbuf() ? m_null : str ),
      m_qstream( qstr ) {}

    //! Save pointer to stream. This function, used in conjunction with reset(),
    //! can be used to pass streams around. This is not possible in general,
    //! since streams are not copyable. See this in action in, e.g.,
    //! Control/Walker/CmdLine/Parser.C.
    //! \return The internal stream buffer of the stream
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    std::streambuf* save() const { return stream<s>().rdbuf(); }

    //! Reset stream to streambuf given. This function, used in conjunction with
    //! save(), can be used to pass streams around. This is not possible in
    //! general, since streams are not copyable. See this in action in, e.g.,
    //! Control/Walker/CmdLine/Parser.C.
    //! \param[in] buf Stream buffer of a stream
    //! \return The internal stream buffer of the stream
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    std::streambuf* reset( std::streambuf* buf ) {
      if (stream<s>().rdbuf() == std::cout.rdbuf())
        m_qstream << "Warning: overwriting std::cout! Doing as requested...\n";
      return stream<s>().rdbuf( buf );
    }

    //! Operator << for printing any type to the verbose stream.
    //! \param[in] os Reference to pretty printer object
    //! \param[in] t Reference to an arbitrary object of type T. T must define
    //! operator<< for std::ostream-compatible streams.
    //! \return The internal stream buffer of the stream
    //! \author J. Bakosi
    template< typename T >
    friend const Print& operator<<( const Print& os, const T& t )
    { os.m_stream << t; return os; }

    //! Operator % for printing any type to the quiet stream.
    //! \param[in] os Reference to pretty printer object
    //! \param[in] t Reference to an arbitrary object of type T. T must define
    //! operator<< for std::ostream-compatible streams.
    //! \return The internal stream buffer of the stream
    //! \author J. Bakosi
    template< typename T >
    friend const Print& operator%( const Print& os, const T& t )
    { os.m_qstream << t; return os; }

    //! Operator % for a function pointer taking ostream returning ostream.
    //! This is so that several of operators of % can be chained together.
    //! \param[in] os Reference to pretty printer object
    //! \param[in] pf Function pointer taking a reference to std::ostream and
    //!   returning a reference to std::ostream
    //! \return Reference to pretty printer object
    //! \author J. Bakosi
    friend const Print& operator%( const Print& os,
      std::ostream& (*pf)(std::ostream&) ) { os.m_qstream << pf; return os; }

    //! Operator << for a function pointer taking ostream returning ostream.
    //! This is so that several of operators of << can be chained together.
    //! \param[in] os Reference to pretty printer object
    //! \param[in] pf Function pointer taking a reference to std::ostream and
    //!   returning a reference to std::ostream
    //! \return Reference to pretty printer object
    //! \author J. Bakosi
    friend const Print& operator<<( const Print& os,
      std::ostream& (*pf)(std::ostream&) ) { os.m_stream << pf; return os; }

    //! Formatted print of part header: title.
    //! \param[in] t Part title to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void part( const std::string& t ) const {
      using std::operator+;
      std::size_t half_length = t.size()/2;
      std::string left( half_length+1, '-' );
      std::string right( t.size()%2 ? half_length+1 : half_length, '-' );
      std::string underline( left + " o " + right );
      std::string upper( t );
      std::transform( begin(t), end(t), begin(upper), ::toupper );
      upper = "< " + upper + " >";
      stream<s>() << m_part_fmt % upper;
      stream<s>() << m_part_underline_fmt % underline;
    }

    //! Formatted print of section header: t.
    //! \param[in] t Section title to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void section( const std::string& t ) const {
      stream<s>() << m_section_title_fmt % m_section_indent % m_section_bullet
                     % t;
      stream<s>() << m_section_underline_fmt % m_section_indent
               % std::string( m_section_indent.size() + 2 + t.size(), '-' );
    }

    //! Formatted print of section header: title : value.
    //! \param[in] name Section title to be printed
    //! \param[in] value Section value to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void section( const std::string& name, const std::string& value ) const {
      stream<s>() << m_section_title_value_fmt % m_section_indent
                     % m_section_bullet % name % value;
      stream<s>() << m_section_underline_fmt % m_section_indent
                     % std::string( m_section_indent.size() + 3 + name.size() +
                                    value.size(), '-' );
    }

    //! Formatted print of subsection header: title.
    //! \param[in] t Subsection title to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void subsection( const std::string& t ) const {
      stream<s>() << m_subsection_title_fmt % m_subsection_indent
                     % m_subsection_bullet % t;
    }

    //! Formatted print of title.
    //! \param[in] value Title string to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void title( const std::string& value ) const {
      // clean up white spaces and format title with no indent or line-break
      auto t = splitLines( value, "", "", 10000 );
      stream<s>() << m_section_title_value_fmt % m_section_indent
                     % m_section_bullet % "Title" % t;
      stream<s>() << m_section_underline_fmt % m_section_indent
                     % std::string( m_section_indent.size()+8+t.size(), '-' );
    }

    //! Formatted print of item: name.
    //! \param[in] name Item name to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void item( const std::string& name ) const
    { stream<s>() << m_item_name_fmt % m_item_indent % name; }

    //! Formatted print of item: name : value
    //! \param[in] name Item name to be printed
    //! \param[in] value Item value to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE, typename T >
    void item( const std::string& name, const T& value ) const
    { stream<s>() << m_item_name_value_fmt % m_item_indent % name % value; }

    //! Formatted print of item: h:m:s.
    //! \param[in] name Item name to be printed
    //! \param[in] watch Watch (in hours, minutes, seconds) to be printed as
    //!   item value
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void item( const std::string& name, const tk::Timer::Watch& watch ) const {
      stream<s>() << m_item_name_watch_fmt % m_item_indent % name
                   % watch.hrs.count() % watch.min.count() % watch.sec.count();
    }

    //! Formatted print of a performance statistic (an item of a list)
    //! \param[in] name Performance statistic name to be printed
    //! \param[in] value Performance statistic value
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void perfitem( const std::string& name, tk::real value ) const
    { stream<s>() << m_item_name_perf_fmt % m_item_indent % name % value; }

    //! Formatted print of a list: name: entries...
    //! \param[in] name Name of a section (consisting of a list) to be printed
    //! \param[in] entries Container of type Container whose elements to be
    //!   printed. Container must be iterable, e.g., possible to be used in a
    //!   range-based for loop. \see
    //!   http://en.cppreference.com/w/cpp/language/range-for
    //! \author J. Bakosi
    template< Style s = VERBOSE, class Container >
    void list( const std::string& name, const Container& entries ) const {
      if (!entries.empty()) {
        section<s>( name );
        for (auto& e : entries)
          stream<s>() << m_list_item_fmt % m_item_indent % e;
      }
    }

    //! Formatted print of a list: name: option names...
    //! \param[in] t Title of the section containing a list
    //! \param[in] factory Factory (an std::map) whose values are printed
    //!   interpreted as options (classes deriving from Toggle), defining the
    //!   name querying member function name().
    //! \author J. Bakosi
    template< class Option, Style s = VERBOSE, class Factory >
    void list( const std::string& t, const Factory& factory ) const {
      if ( !factory.empty() ) {
        section<s>( t );
        Option option;
        for (const auto& f : factory)
          stream<s>() << m_list_item_fmt % m_item_indent % option.name(f.first);
      }
    }

    //! Formatted print of elapsed times
    //! \param[in] t Title of section containing a list of elapsed times
    //! \param[in] clock std::vector of strings (clock names) and associated
    //!   timers which could be in various formats as long as there is a
    //!   corresponding item() overload that can apply operator << for outputing
    //!   their value to an output stream. Examples of allowed ClockFormats are:
    //!   tk::Timer::Watch, which is a struct containing a timestamp in h:m:s
    //!   format, and the return value of Timer::dsec(), which is a tk::real.
    //! \author J. Bakosi
    template< Style s = VERBOSE, class ClockFormat >
    void time( const std::string& t,
               const std::vector<
                 std::pair< std::string, ClockFormat > >& clock ) const
    {
      section<s>( t );
      for (const auto& c : clock) item<s>( c.first, c.second );
    }

    //! Formatted print of performance statistics
    //! \param[in] t Title of section containing a list of performance stats
    //! \param[in] stat std::vector of strings (names of a performance
    //!   statistics) and associated values.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void perf( const std::string& t,
               const std::vector< std::pair< std::string, tk::real > >& stat )
    const
    {
      if (!stat.empty()) {
        section<s>( t );
        for (const auto& c : stat) perfitem<s>( c.first, c.second );
      }
    }

    //! Formatted print of a note
    //! \param[in] msg Message to print as a note
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void note( const std::string& msg ) const
    { stream<s>() << m_note_fmt % m_item_indent % msg; }

    //! Echo formatted print of a diagnostics message
    //! \param[in] msg Message to print as a diagnostics message
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void diag( const std::string& msg ) const
    { stream<s>() << m_diag_fmt % msg << std::flush; }

    //! Start formatted print of a diagnostics message
    //! \param[in] msg First part of message to print as a diagnostics message
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void diagstart( const std::string& msg ) const
    { stream<s>() << m_diag_start_fmt % msg << std::flush; }

    //! Finish formatted print of a diagnostics message
    //! \param[in] msg Last part of message to print as a diagnostics message
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void diagend( const std::string& msg ) const
    { stream<s>() << m_diag_end_fmt % msg << std::flush; }

    //! Echo formatted print of a progress message
    //! \param[in] prefix Strings to output prefixing the progress report
    //! \param[in] done Array of integers indicating how many have been done
    //! \param[in] max Array of integers indicating how many to be done
    //! \param[in] progress_size Size of previous progress report (to overwrite)
    //! \details All input arrays are the same size. The prefix strings
    //!   are optional, i.e., they can be empty strings. The function generates
    //!   an output to the stream configured in the following fashion:
    //!   <pre1>[<done1>/<max1>], <pre2>[<done2>/<max2>], ..., e.g.,
    //!   r:[1/3], b[2/8]. Whenever this function is called, a number of
    //!   backspaces are put into the stream so that the new progress report
    //!   string overwrites the old one. In order to backtrack the correct
    //!   amount, the length of the old progress report is stored (by whatever
    //!   object holds us) and passed in by reference in progress_size, which is
    //!   overwritten here once it has been used for backtracking. Therefore,
    //!   for restarting a new series of progress reports, this variable must be
    //!   zeroed. Also, it is best to not to interleave multiple tasks, because
    //!   even if a different progress_size is kept for each, there is no regard
    //!   as to which line we output to in the stream. In other words, multiple
    //!   task outputs will be intermingled, leading to confusing screen output.
    //! \author J. Bakosi
    template< std::size_t N, Style s = VERBOSE >
    void progress( const std::array< std::string, N >& prefix,
                   const std::array< int, N >& done,
                   const std::array< int, N >& max,
                   std::size_t& progress_size ) const
    {
      // lambda to determine the number of digits in an integer
      auto numdig = []( int i ) -> std::size_t {
        return i > 0 ?
          static_cast< std::size_t >( std::log10(static_cast<double>(i)) ) + 1
          : 1; };
      // Backspace so that new progress can overwrite old one
      stream<s>() << std::string( progress_size, '\b' );
      std::stringstream ss;
      auto ip = prefix.cbegin();
      auto id = done.cbegin();
      auto im = max.cbegin();
      progress_size = 0;
      while (ip != prefix.cend()) {
        // Compute new length of progress string
        progress_size += 4 + ip->size() + numdig(*id) + numdig(*im);
        // Construct and output new progress string to stream
        ss << *ip << ":[" << *id << '/' << *im << ']';
        ++ip; ++id; ++im;
        // if next subprogress is not the last one, put in a comma
        if (ip != prefix.cend()) {
          ss << ", ";
          progress_size += 2;
        } else {
          ss << ' ';
          ++progress_size;
        }
      }
      stream<s>() << m_progress_fmt % ss.str() << std::flush;
    }

    //! \brief Formatted print of help of one-liners on all command-line
    //!   parameters or control file keywords
    //! \param[in] executable Name of executable to output help for
    //! \param[in] pool std::map of keywords and their associated information
    //! \param[in] msg Message to print after exectuable in the title
    //! \param[in] pfx Prefix in front of alias, double prefix in front of
    //!   keyword
    //! \author J. Bakosi
    template< Style s = VERBOSE, class Help >
    void help( std::string executable,
               const Help& pool,
               const std::string& msg,
               const std::string& pfx = "" ) const
    {
      stream<s>() << m_help_title_fmt % executable % msg;
      for (const auto& keyword : pool) {
        const auto& info = keyword.second;
        const auto& alias = info.alias;
        const auto& expt = info.expt;
        stream<s>() << m_help_item_fmt
                       % std::string( ( alias ? pfx + *alias + ", " : "") +
                                        pfx + pfx + keyword.first )
                       % (expt ? *expt : "")
                       % info.shortDescription;
      }
    }

    //! \brief Formatted print of verbose help on a single command-line
    //!   parameter or control file keywords
    //! \param[in] executable Name of executable to output help for
    //! \param[in] kw Keyword help struct on which help is to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE, class HelpKw >
    void helpkw( std::string executable, const HelpKw& kw ) const {
      Assert( !kw.keyword.empty(), "Empty keyword in Print::helpkw()" );
      const auto& info = kw.info;
      const auto& alias = info.alias;
      const auto& expt = info.expt;
      const auto& choices = info.choices;
      // print keyword title
      if (kw.cmd)
        stream<s>() << m_helpkw_cmd_title_fmt
                       % executable
                       % (alias ? "-" + *alias + ", " : "")
                       % kw.keyword;
      else
        stream<s>() << m_helpkw_ctr_title_fmt
                       % executable
                       % kw.keyword;
      // print short description
      stream<s>() << m_description_fmt
                     % splitLines( info.shortDescription, m_subsection_indent );
      // print long description
      stream<s>() << m_description_fmt
                     % splitLines( info.longDescription, m_subsection_indent );
      // print expected type description
      if (expt)
        stream<s>() << m_description_fmt
                       % splitLines( *expt, m_subsection_indent,
                                     "Expected type: " );
      // print expected valied choices
      if (choices)
        stream<s>() << m_description_fmt
                    % splitLines( *choices, m_subsection_indent,
                                  "Expected valid choices: ");
    }

    //! Print end of a part
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void endpart() const { stream<s>() << std::endl; }

    //! Print end of subsection
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void endsubsection() const { stream<s>() << '\n'; }

    //! Print raw data to stream.
    //! \param[in] r Arbitrary data of arbitrary type as long as it defines
    //!   operator << for std::ostream.
    //! \author J. Bakosi
    template< Style s = VERBOSE, typename T >
    void raw( const T& r ) const { stream<s>() << r; }

    //! Return verbose or quiet stream depending on style template argument.
    //! Non-const version.
    //! \return Reference to underlying std::ostream.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    std::ostream& stream() noexcept { return s ? m_stream : m_qstream; }

    //! Return verbose or quiet stream depending on style template argument.
    //! Const version.
    //! \return Reference to underlying std::ostream.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    std::ostream& stream() const noexcept { return s ? m_stream : m_qstream; }

    //! Function object for echoing policies to screen
    //! \author J. Bakosi
    struct echoPolicies {
      //! Need to store reference to host class whose data we operate on
      const Print* const m_host;
      //! Constructor: store host object pointer
      echoPolicies( const Print* const host ) : m_host( host ) {}
      //! Function call operator templated on the type that echos a policy
      template< typename U > void operator()( U ) {
        static_assert( tk::HasTypedefCode< typename U::info >::value,
                       "Policy code undefined for keyword" );
        // Print policy code - policy name
        m_host->raw( m_host->m_item_indent + "   " +
                     *U::code() + " - " + U::info::name() + '\n' );

      }
    };

    //! Print Inciter header. Text ASCII Art Generator used for executable
    //! names: http://patorjk.com/software/taag, Picture ASCII Art Generator
    //! used for converting the logo text "Quinoa": http://picascii.com.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void headerInciter() const {
       stream<s>() << R"(
      ,::,`                                                            `.
   .;;;'';;;:                                                          ;;#
  ;;;@+   +;;;  ;;;;;,   ;;;;. ;;;;;, ;;;;      ;;;;   `;;;;;;:        ;;;
 :;;@`     :;;' .;;;@,    ,;@, ,;;;@: .;;;'     .;+;. ;;;@#:';;;      ;;;;'
 ;;;#       ;;;: ;;;'      ;:   ;;;'   ;;;;;     ;#  ;;;@     ;;;     ;+;;'
.;;+        ;;;# ;;;'      ;:   ;;;'   ;#;;;`    ;#  ;;@      `;;+   .;#;;;.
;;;#        :;;' ;;;'      ;:   ;;;'   ;# ;;;    ;# ;;;@       ;;;   ;# ;;;+
;;;#        .;;; ;;;'      ;:   ;;;'   ;# ,;;;   ;# ;;;#       ;;;:  ;@  ;;;
;;;#        .;;' ;;;'      ;:   ;;;'   ;#  ;;;;  ;# ;;;'       ;;;+ ;',  ;;;@
;;;+        ,;;+ ;;;'      ;:   ;;;'   ;#   ;;;' ;# ;;;'       ;;;' ;':::;;;;
`;;;        ;;;@ ;;;'      ;:   ;;;'   ;#    ;;;';# ;;;@       ;;;:,;+++++;;;'
 ;;;;       ;;;@ ;;;#     .;.   ;;;'   ;#     ;;;;# `;;+       ;;# ;#     ;;;'
 .;;;      :;;@  ,;;+     ;+    ;;;'   ;#      ;;;#  ;;;      ;;;@ ;@      ;;;.
  ';;;    ;;;@,   ;;;;``.;;@    ;;;'   ;+      .;;#   ;;;    :;;@ ;;;      ;;;+
   :;;;;;;;+@`     ';;;;;'@    ;;;;;, ;;;;      ;;+    +;;;;;;#@ ;;;;.   .;;;;;;
     .;;#@'         `#@@@:     ;::::; ;::::      ;@      '@@@+   ;:::;    ;::::::
    :;;;;;;.     .___              .__  __
   .;@+@';;;;;;' |   | ____   ____ |__|/  |_  ___________
    `     '#''@` |   |/    \_/ ___\|  \   __\/ __ \_  __ \
                 |   |   |  \  \___|  ||  | \  ___/|  | \/
                 |___|___|  /\___  >__||__|  \___  >__|
                          \/     \/              \/)"
      << std::endl;
    }

    //! Print RNGTest header. Text ASCII Art Generator used for executable
    //! names: http://patorjk.com/software/taag, Picture ASCII Art Generator
    //! used for converting the logo text "Quinoa": http://picascii.com.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void headerRNGTest() const {
       stream<s>() << R"(
      ,::,`                                                            `.
   .;;;'';;;:                                                          ;;#
  ;;;@+   +;;;  ;;;;;,   ;;;;. ;;;;;, ;;;;      ;;;;   `;;;;;;:        ;;;
 :;;@`     :;;' .;;;@,    ,;@, ,;;;@: .;;;'     .;+;. ;;;@#:';;;      ;;;;'
 ;;;#       ;;;: ;;;'      ;:   ;;;'   ;;;;;     ;#  ;;;@     ;;;     ;+;;'
.;;+        ;;;# ;;;'      ;:   ;;;'   ;#;;;`    ;#  ;;@      `;;+   .;#;;;.
;;;#        :;;' ;;;'      ;:   ;;;'   ;# ;;;    ;# ;;;@       ;;;   ;# ;;;+
;;;#        .;;; ;;;'      ;:   ;;;'   ;# ,;;;   ;# ;;;#       ;;;:  ;@  ;;;
;;;#        .;;' ;;;'      ;:   ;;;'   ;#  ;;;;  ;# ;;;'       ;;;+ ;',  ;;;@
;;;+        ,;;+ ;;;'      ;:   ;;;'   ;#   ;;;' ;# ;;;'       ;;;' ;':::;;;;
`;;;        ;;;@ ;;;'      ;:   ;;;'   ;#    ;;;';# ;;;@       ;;;:,;+++++;;;'
 ;;;;       ;;;@ ;;;#     .;.   ;;;'   ;#     ;;;;# `;;+       ;;# ;#     ;;;'
 .;;;      :;;@  ,;;+     ;+    ;;;'   ;#      ;;;#  ;;;      ;;;@ ;@      ;;;.
  ';;;    ;;;@,   ;;;;``.;;@    ;;;'   ;+      .;;#   ;;;    :;;@ ;;;      ;;;+
   :;;;;;;;+@`     ';;;;;'@    ;;;;;, ;;;;      ;;+    +;;;;;;#@ ;;;;.   .;;;;;;
     .;;#@'         `#@@@:     ;::::; ;::::      ;@      '@@@+   ;:::;    ;::::::
    :;;;;;;.     __________ _______    __________________             __
   .;@+@';;;;;;' \______   \\      \  /  _____\__    _______   ______/  |_
    `     '#''@`  |       _//   |   \/   \  ___ |    |_/ __ \ /  ___\   __\
                  |    |   /    |    \    \_\  \|    |\  ___/ \___ \ |  |
                  |____|_  \____|__  /\______  /|____| \___  /____  >|__|
                         \/        \/        \/            \/     \/)"
      << std::endl;
    }

    //! Print UnitTest header. Text ASCII Art Generator used for executable
    //! names: http://patorjk.com/software/taag, Picture ASCII Art Generator
    //! used for converting the logo text "Quinoa": http://picascii.com.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void headerUnitTest() const {
       stream<s>() << R"(
      ,::,`                                                            `.
   .;;;'';;;:                                                          ;;#
  ;;;@+   +;;;  ;;;;;,   ;;;;. ;;;;;, ;;;;      ;;;;   `;;;;;;:        ;;;
 :;;@`     :;;' .;;;@,    ,;@, ,;;;@: .;;;'     .;+;. ;;;@#:';;;      ;;;;'
 ;;;#       ;;;: ;;;'      ;:   ;;;'   ;;;;;     ;#  ;;;@     ;;;     ;+;;'
.;;+        ;;;# ;;;'      ;:   ;;;'   ;#;;;`    ;#  ;;@      `;;+   .;#;;;.
;;;#        :;;' ;;;'      ;:   ;;;'   ;# ;;;    ;# ;;;@       ;;;   ;# ;;;+
;;;#        .;;; ;;;'      ;:   ;;;'   ;# ,;;;   ;# ;;;#       ;;;:  ;@  ;;;
;;;#        .;;' ;;;'      ;:   ;;;'   ;#  ;;;;  ;# ;;;'       ;;;+ ;',  ;;;@
;;;+        ,;;+ ;;;'      ;:   ;;;'   ;#   ;;;' ;# ;;;'       ;;;' ;':::;;;;
`;;;        ;;;@ ;;;'      ;:   ;;;'   ;#    ;;;';# ;;;@       ;;;:,;+++++;;;'
 ;;;;       ;;;@ ;;;#     .;.   ;;;'   ;#     ;;;;# `;;+       ;;# ;#     ;;;'
 .;;;      :;;@  ,;;+     ;+    ;;;'   ;#      ;;;#  ;;;      ;;;@ ;@      ;;;.
  ';;;    ;;;@,   ;;;;``.;;@    ;;;'   ;+      .;;#   ;;;    :;;@ ;;;      ;;;+
   :;;;;;;;+@`     ';;;;;'@    ;;;;;, ;;;;      ;;+    +;;;;;;#@ ;;;;.   .;;;;;;
     .;;#@'         `#@@@:     ;::::; ;::::      ;@      '@@@+   ;:::;    ;::::::
    :;;;;;;.      ____ ___      .__  __ ___________              __
   .;@+@';;;;;;' |    |   \____ |__|/  |\__    ___/___   _______/  |_
    `     '#''@` |    |   /    \|  \   __\|    |_/ __ \ /  ___/\   __\
                 |    |  /   |  \  ||  |  |    |\  ___/ \___ \  |  |
                 |______/|___|  /__||__|  |____| \___  >____  > |__|
                              \/                     \/     \/)"
      << std::endl;
    }

    //! Print MeshConv header. Text ASCII Art Generator used for executable
    //! names: http://patorjk.com/software/taag, Picture ASCII Art Generator
    //! used for converting the logo text "Quinoa": http://picascii.com.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void headerMeshConv() const {
      stream<s>() << R"(
      ,::,`                                                            `.
   .;;;'';;;:                                                          ;;#
  ;;;@+   +;;;  ;;;;;,   ;;;;. ;;;;;, ;;;;      ;;;;   `;;;;;;:        ;;;
 :;;@`     :;;' .;;;@,    ,;@, ,;;;@: .;;;'     .;+;. ;;;@#:';;;      ;;;;'
 ;;;#       ;;;: ;;;'      ;:   ;;;'   ;;;;;     ;#  ;;;@     ;;;     ;+;;'
.;;+        ;;;# ;;;'      ;:   ;;;'   ;#;;;`    ;#  ;;@      `;;+   .;#;;;.
;;;#        :;;' ;;;'      ;:   ;;;'   ;# ;;;    ;# ;;;@       ;;;   ;# ;;;+
;;;#        .;;; ;;;'      ;:   ;;;'   ;# ,;;;   ;# ;;;#       ;;;:  ;@  ;;;
;;;#        .;;' ;;;'      ;:   ;;;'   ;#  ;;;;  ;# ;;;'       ;;;+ ;',  ;;;@
;;;+        ,;;+ ;;;'      ;:   ;;;'   ;#   ;;;' ;# ;;;'       ;;;' ;':::;;;;
`;;;        ;;;@ ;;;'      ;:   ;;;'   ;#    ;;;';# ;;;@       ;;;:,;+++++;;;'
 ;;;;       ;;;@ ;;;#     .;.   ;;;'   ;#     ;;;;# `;;+       ;;# ;#     ;;;'
 .;;;      :;;@  ,;;+     ;+    ;;;'   ;#      ;;;#  ;;;      ;;;@ ;@      ;;;.
  ';;;    ;;;@,   ;;;;``.;;@    ;;;'   ;+      .;;#   ;;;    :;;@ ;;;      ;;;+
   :;;;;;;;+@`     ';;;;;'@    ;;;;;, ;;;;      ;;+    +;;;;;;#@ ;;;;.   .;;;;;;
     .;;#@'         `#@@@:     ;::::; ;::::      ;@      '@@@+   ;:::;    ;::::::
    :;;;;;;.        _____                .__    _________
   .;@+@';;;;;;'   /     \   ____   _____|  |__ \_   ___ \  ____   _______  __
    `     '#''@`  /  \ /  \_/ __ \ /  ___|  |  \/    \  \/ /  _ \ /    \  \/ /
                 /    Y    \  ___/ \___ \|   Y  \     \___(  <_> |   |  \   /
                 \____|__  /\___  /____  |___|  /\______  /\____/|___|  /\_/
                         \/     \/     \/     \/        \/            \/)"
      << std::endl;
    }

    //! Print Walker header. Text ASCII Art Generator used for executable names:
    //! http://patorjk.com/software/taag, Picture ASCII Art Generator used for
    //! converting the logo text "Quinoa": http://picascii.com.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void headerWalker() const {
      stream<s>() << R"(
      ,::,`                                                            `.
   .;;;'';;;:                                                          ;;#
  ;;;@+   +;;;  ;;;;;,   ;;;;. ;;;;;, ;;;;      ;;;;   `;;;;;;:        ;;;
 :;;@`     :;;' .;;;@,    ,;@, ,;;;@: .;;;'     .;+;. ;;;@#:';;;      ;;;;'
 ;;;#       ;;;: ;;;'      ;:   ;;;'   ;;;;;     ;#  ;;;@     ;;;     ;+;;'
.;;+        ;;;# ;;;'      ;:   ;;;'   ;#;;;`    ;#  ;;@      `;;+   .;#;;;.
;;;#        :;;' ;;;'      ;:   ;;;'   ;# ;;;    ;# ;;;@       ;;;   ;# ;;;+
;;;#        .;;; ;;;'      ;:   ;;;'   ;# ,;;;   ;# ;;;#       ;;;:  ;@  ;;;
;;;#        .;;' ;;;'      ;:   ;;;'   ;#  ;;;;  ;# ;;;'       ;;;+ ;',  ;;;@
;;;+        ,;;+ ;;;'      ;:   ;;;'   ;#   ;;;' ;# ;;;'       ;;;' ;':::;;;;
`;;;        ;;;@ ;;;'      ;:   ;;;'   ;#    ;;;';# ;;;@       ;;;:,;+++++;;;'
 ;;;;       ;;;@ ;;;#     .;.   ;;;'   ;#     ;;;;# `;;+       ;;# ;#     ;;;'
 .;;;      :;;@  ,;;+     ;+    ;;;'   ;#      ;;;#  ;;;      ;;;@ ;@      ;;;.
  ';;;    ;;;@,   ;;;;``.;;@    ;;;'   ;+      .;;#   ;;;    :;;@ ;;;      ;;;+
   :;;;;;;;+@`     ';;;;;'@    ;;;;;, ;;;;      ;;+    +;;;;;;#@ ;;;;.   .;;;;;;
     .;;#@'         `#@@@:     ;::::; ;::::      ;@      '@@@+   ;:::;    ;::::::
    :;;;;;;.      __      __        .__   __
  .;@+@';;;;;;'  /  \    /  \_____  |  | |  | __ ___________
    `     '#''@` \   \/\/   /\__  \ |  | |  |/ // __ \_  __ \
                  \        /  / __ \|  |_|    <\  ___/|  | \/
                   \__/\  /  (____  /____/__|_ \\___  >__|
                        \/        \/          \/    \/)"
      << std::endl;
    }

  protected:
    //! Bullets
    const char m_section_bullet = '*';
    const char m_subsection_bullet = '<';
    //! Indents
    const std::string m_section_indent = " ";
    const std::string m_subsection_indent =
      std::operator+(m_section_indent,"  ");
    const std::string m_item_indent = std::operator+(m_subsection_indent,"  ");

    //! Format strings. See http://www.boost.org/doc/libs/release/libs/format.
    using format = boost::format;
    mutable format m_header_fmt = format("%|=80|\n");
    mutable format m_part_fmt = format("\n%|=80|\n");
    mutable format m_section_title_fmt = format("\n%s%c %s:\n");
    mutable format m_section_title_value_fmt = format("\n%s%c %s: %s\n");
    mutable format m_subsection_title_fmt = format("%s%c %s >\n");
    mutable format m_list_item_fmt = format("%s%-30s\n");
    mutable format m_note_fmt = format("%s%-30s\n");
    mutable format m_diag_fmt = format("Quinoa> %s\n");
    mutable format m_diag_start_fmt = format("Quinoa> %s ");
    mutable format m_diag_end_fmt = format("%s\n");
    mutable format m_progress_fmt = format("%s");
    mutable format m_help_title_fmt = format("\n%s %s\n");
    mutable format m_help_item_fmt = format("%20s%11s %s\n");
    mutable format m_helpkw_cmd_title_fmt =
              format("\n%s command-line keyword %s--%s\n\n");
    mutable format m_helpkw_ctr_title_fmt =
              format("\n%s control file keyword '%s'\n\n");
    mutable format m_helpkw_fmt = format("%s%s\n\n%s%s\n\n");
    mutable format m_description_fmt = format("%s\n\n");
    mutable format m_item_name_fmt = format("%s%-30s : ");
    mutable format m_item_name_value_fmt = format("%s%-30s : %s\n");
    mutable format m_item_name_watch_fmt = format("%s%-75s : %d:%d:%d\n");
    mutable format m_item_name_perf_fmt = format("%s%-75s : %s\n");
    mutable format m_item_widename_value_fmt = format("%s%-75s : %s\n");
    mutable format m_part_underline_fmt = format("      %|=68|\n");
    mutable format m_section_underline_fmt = format("%s%s\n");

    // Stream objects
    std::stringstream m_null;   //!< Default verbose stream
    std::ostream& m_stream;     //!< Verbose stream
    std::ostream& m_qstream;    //!< Quiet stream

  private:
    //! \brief Clean up whitespaces and format a long string into multiple lines
    //! \param[in] str String to format
    //! \param[in] name String to insert before string to output
    //! \param[in] indent String to use as identation
    //! \param[in] width Width in characters to insert newlines for output
    //! \author J. Bakosi
    //! \see http://stackoverflow.com/a/6892562
    //! \see http://stackoverflow.com/a/8362145
    // TODO A line longer than 'width' will cause a hang!
    std::string splitLines( std::string str,
                            std::string indent,
                            std::string name = "",
                            std::size_t width = 80 ) const {
      // remove form feeds, line feeds, carriage returns, horizontal tabs,
      // vertical tabs, see http://en.cppreference.com/w/cpp/string/byte/isspace
      str.erase(
        std::remove_if( str.begin(), str.end(),
                        []( char x ){ return std::isspace( x ) && x != ' '; } ),
        str.end() );
      // remove duplicate spaces
      str.erase(
        std::unique( str.begin(), str.end(),
                     []( char a, char b ){ return a == b && a == ' '; } ),
        str.end() );
      // format str to 'width'-character-long lines with indent
      str.insert( 0, indent + name );
      std::size_t currIndex = width - 1;
      std::size_t sizeToElim;
      while ( currIndex < str.length() ) {
        const std::string whitespace = " ";
        currIndex = str.find_last_of( whitespace, currIndex + 1 );
        if ( currIndex == std::string::npos ) break;
        currIndex = str.find_last_not_of( whitespace, currIndex );
        if ( currIndex == std::string::npos ) break;
        sizeToElim =
          str.find_first_not_of( whitespace, currIndex + 1 ) - currIndex - 1;
        str.replace( currIndex + 1, sizeToElim , "\n" + indent );
        currIndex += width + indent.length() + 1;
      }
      return str;
    }
};

} // tk::

#endif // Print_h
