//******************************************************************************
/*!
  \file      src/Base/Print.h
  \author    J. Bakosi
  \date      Tue 17 Feb 2015 04:02:08 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     General purpose pretty printer functionality
  \details   This file contains general purpose printer functions. Using the
    functions defined here provides formatting, and a consistent look with
    simple client-side code. For formatting, the Boost Format library is used,
    see http://www.boost.org/doc/libs/release/libs/format.
*/
//******************************************************************************
#ifndef Print_h
#define Print_h

#include <iostream>
#include <sstream>
#include <iomanip>
#include <list>
#include <map>

#include <boost/format.hpp>

#include <Timer.h>
#include <Exception.h>

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
    //! \param[inout] str Verbose stream
    //! \param[inout] qstr Quiet stream
    //! \author J. Bakosi
    explicit Print( std::ostream& str = std::clog,
                    std::ostream& qstr = std::cout ) :
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

    //! Operator << for a function pointer taking ostream returning ostream.
    //! This is so that several of operators of << can be chained together.
    //! \param[in] os Reference to pretty printer object
    //! \param[in] pf Function pointer taking a reference to std::ostream and
    //!   returning a reference to std::ostream
    //! \return Reference to pretty printer object
    //! \author J. Bakosi
    friend const Print& operator<<( const Print& os,
      std::ostream& (*pf)(std::ostream&) ) { os.m_stream << pf; return os; }

    //! Operator % for a function pointer taking ostream returning ostream.
    //! This is so that several of operators of % can be chained together.
    //! \param[in] os Reference to pretty printer object
    //! \param[in] pf Function pointer taking a reference to std::ostream and
    //!   returning a reference to std::ostream
    //! \return Reference to pretty printer object
    //! \author J. Bakosi
    friend const Print& operator%( const Print& os,
      std::ostream& (*pf)(std::ostream&) ) { os.m_qstream << pf; return os; }

    //! Formatted print of part header: title.
    //! \param[in] title Part title to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void part( const std::string& title ) const {
      using std::operator+;
      std::size_t half_length = title.size()/2 + 1;
      std::string str( half_length, '-' );
      std::string underline( str + " o " + str );
      std::string upper( title );
      std::transform( begin(title), end(title), begin(upper), ::toupper );
      upper = "< " + upper + " >";
      stream<s>() << m_part_fmt % upper;
      stream<s>() << m_part_underline_fmt % underline;
    }

    //! Formatted print of section header: title.
    //! \param[in] title Section title to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void section( const std::string& title ) const {
      stream<s>() << m_section_title_fmt % m_section_indent % m_section_bullet
                     % title;
      stream<s>() << m_section_underline_fmt % m_section_indent
               % std::string( m_section_indent.size() + 2 + title.size(), '-' );
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
    //! \param[in] title Subsection title to be printed
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void subsection( const std::string& title ) const {
      stream<s>() << m_subsection_title_fmt % m_subsection_indent
                     % m_subsection_bullet % title;
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
    //! \param[in] title Title of the section containing a list
    //! \param[in] factory Factory (an std::map) whose values are printed
    //!   interpreted as options (classes deriving from Toggle), defining the
    //!   name querying member function name().
    //! \author J. Bakosi
    template< class Option, Style s = VERBOSE, class Factory >
    void list( const std::string& title, const Factory& factory ) const {
      if ( !factory.empty() ) {
        section<s>( title );
        Option option;
        for (const auto& f : factory)
          stream<s>() << m_list_item_fmt % m_item_indent % option.name(f.first);
      }
    }

    //! Formatted print of elapsed times
    //! \param[in] title Title of section containing a list of elapsed times
    //! \param[in] clock std::map of strings (clock names) and associated timers
    //!   which could be in various formats as long as there is a corresponding
    //!   item() overload that can apply operator << for outputing their value
    //!   to an output stream. Examples of allowed ClockFormats are:
    //!   tk::Timer::Watch, which is a struct containing a timestamp in h:m:s
    //!   format, and the return value of Timer::dsec(), which is a tk::real.
    //! \author J. Bakosi
    template< Style s = VERBOSE, class ClockFormat >
    void time( const std::string& title,
               const std::map< std::string, ClockFormat >& clock ) const
    {
      section<s>( title );
      for (const auto& c : clock) item<s>( c.first, c.second );
      endsubsection<s>();
    }

    //! Formatted print of a note
    //! \param[in] msg Message to print as a note
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void note( const std::string& msg ) const
    { stream<s>() << m_note_fmt % m_section_indent % msg; }

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
      stream<s>() << m_description_fmt % splitLines(info.shortDescription);
      // print long description
      stream<s>() << m_description_fmt % splitLines(info.longDescription);
      // print expected type description
      if (expt)
        stream<s>() << m_description_fmt % splitLines(*expt, "Expected type: ");
      // print expected valied choices
      if (choices)
        stream<s>() << m_description_fmt
                    % splitLines(*choices, "Expected valid choices: ");
    }

    //! Print end of a part
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void endpart() const { stream<s>() << '\n'; }

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

    //! Print Quinoa header. Text ASCII Art Generator used for executable
    //! names: http://patorjk.com/software/taag, Picture ASCII Art Generator
    //! used for converting the logo text "Quinoa": http://picascii.com.
    //! \author J. Bakosi
    template< Style s = VERBOSE >
    void headerQuinoa() const {
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
    :;;;;;;.
   .;@+@';;;;;;'
    `     '#''@`                                                               )"
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
    :;;;;;;.     __________ _______    __________________
   .;@+@';;;;;;' \______   \\      \  /  _____\__    _______   ______/  |_
    `     '#''@`  |       _//   |   \/   \  ___ |    |_/ __ \ /  ___\   __\
                  |    |   /    |    \    \_\  \|    |\  ___/ \___ \ |  |
                  |____|_  \____|__  /\______  /|____| \___  /____  >|__|
                         \/        \/        \/            \/     \/         )"
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
                              \/                     \/     \/            )"
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
                         \/     \/     \/     \/        \/            \/       )"
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
                        \/        \/          \/    \/                       )"
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
    mutable format m_item_name_watch_fmt = format("%s%-65s : %d:%d:%d\n");
    mutable format m_item_widename_value_fmt = format("%s%-65s : %s\n");
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
    //! \param[in] witdth Width in characters to insert newlines for output
    //! \author J. Bakosi
    //! \see http://stackoverflow.com/a/6892562
    //! \see http://stackoverflow.com/a/8362145
    // TODO A line longer than 'width' will cause a hang!
    std::string splitLines( std::string str,
                            std::string name = "",
                            std::size_t width = 80 ) const
    {
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
      // format str to 'witdh'-character-long lines with indent
      const auto& indent = m_subsection_indent;
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
