//******************************************************************************
/*!
  \file      src/Base/Print.h
  \author    J. Bakosi
  \date      Sun 24 Aug 2014 09:11:16 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Print
  \details   Print
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

namespace tk {

//! Output verbosity. C-style enum as this is used for template argument.
enum Style { QUIET=0, VERBOSE=1 };

//! Print base
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
    explicit Print( std::ostream& str = std::clog,
                    std::ostream& qstr = std::cout ) :
      m_stream( str.rdbuf() == std::clog.rdbuf() ? m_null : str ),
      m_qstream( qstr ) {}

    //! Save pointer to stream
    template< Style s = VERBOSE >
    std::streambuf* save() const { return stream<s>().rdbuf(); }

    //! Reset stream to streambuf
    template< Style s = VERBOSE >
    std::streambuf* reset( std::streambuf* buf ) {
      if (stream<s>().rdbuf() == std::cout.rdbuf())
        m_qstream << "Warning: overwriting std::cout! Doing as requested...\n";
      return stream<s>().rdbuf( buf );
    }

    //! Operator << for printing any type to the verbose stream
    template< typename T >
    friend const Print& operator<<( const Print& os, const T& t )
    { os.m_stream << t; return os; }

    //! Operator % for printing any type to the quiet stream
    template< typename T >
    friend const Print& operator%( const Print& os, const T& t )
    { os.m_qstream << t; return os; }

    //! Operator << for a function pointer taking ostream returning ostream.
    //! This is so that several of operators of << can be chained together.
    friend const Print& operator<<( const Print& os,
      std::ostream& (*pf)(std::ostream&) ) { os.m_stream << pf; return os; }

    //! Operator % for a function pointer taking ostream returning ostream.
    //! This is so that several of operators of % can be chained together.
    friend const Print& operator%( const Print& os,
      std::ostream& (*pf)(std::ostream&) ) { os.m_qstream << pf; return os; }

    //! Print part header: title
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

    //! Print section header: title
    template< Style s = VERBOSE >
    void section( const std::string& title ) const {
      stream<s>() << m_section_title_fmt % m_section_indent % m_section_bullet
                     % title;
      stream<s>() << m_section_underline_fmt % m_section_indent
               % std::string( m_section_indent.size() + 2 + title.size(), '-' );
    }

    //! Print section header: title : value
    template< Style s = VERBOSE >
    void section( const std::string& name, const std::string& value ) const {
      stream<s>() << m_section_title_value_fmt % m_section_indent
                     % m_section_bullet % name % value;
      stream<s>() << m_section_underline_fmt % m_section_indent
                     % std::string( m_section_indent.size() + 3 + name.size() +
                                    value.size(), '-' );
    }

    //! Print subsection header: title
    template< Style s = VERBOSE >
    void subsection( const std::string& title ) const {
      stream<s>() << m_subsection_title_fmt % m_subsection_indent
                     % m_subsection_bullet % title;
    }

    //! Print item: name
    template< Style s = VERBOSE >
    void item( const std::string& name ) const
    { stream<s>() << m_item_name_fmt % m_item_indent % name; }

    //! Print item: name : value
    template< Style s = VERBOSE, typename T >
    void item( const std::string& name, const T& value ) const
    { stream<s>() << m_item_name_value_fmt % m_item_indent % name % value; }

    //! Print item: h:m:s
    template< Style s = VERBOSE >
    void item( const std::string& name, const tk::Timer::Watch& watch ) const {
      stream<s>() << m_item_name_watch_fmt % m_item_indent % name
                   % watch.hrs.count() % watch.min.count() % watch.sec.count();
    }

    //! Print list: name: entries...
    template< Style s = VERBOSE, class Container >
    void list( const std::string& name, const Container& entries ) const {
      if (!entries.empty()) {
        section<s>( name );
        for (auto& e : entries)
          stream<s>() << m_list_item_fmt % m_item_indent % e;
      }
    }

    //! Print list: name: option names...
    template< class Option, Style s = VERBOSE, class Factory >
    void list( const std::string& title, const Factory& factory ) const {
      if ( !factory.empty() ) {
        section<s>( title );
        Option option;
        for (const auto& f : factory)
          stream<s>() << m_list_item_fmt % m_item_indent % option.name(f.first);
      }
    }

    //! Print elapsed time
    template< Style s = VERBOSE, class ClockFormat >
    void time( const std::string& title,
               const std::map< std::string, ClockFormat >& clock ) const
    {
      section<s>( title );
      for (const auto& c : clock) item<s>( c.first, c.second );
      endsubsection<s>();
    }

    //! Print note
    template< Style s = VERBOSE >
    void note( const std::string& msg ) const
    { stream<s>() << m_note_fmt % m_section_indent % msg; }

    //! Print end of part
    template< Style s = VERBOSE >
    void endpart() const { stream<s>() << '\n'; }

    //! Print end of subsection
    template< Style s = VERBOSE >
    void endsubsection() const { stream<s>() << '\n'; }

    //! Print raw
    template< Style s = VERBOSE, typename T >
    void raw( const T& r ) const { stream<s>() << r; }

    //! Return verbose or quiet stream depending on style template argument
    template< Style s = VERBOSE >
    std::ostream& stream() noexcept { return s ? m_stream : m_qstream; }

    //! Return verbose or quiet stream depending on style template argument
    template< Style s = VERBOSE >
    std::ostream& stream() const noexcept { return s ? m_stream : m_qstream; }

    //! Print Quinoa header
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

    //! Print RNGTest header
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
    :;;;;;;.       __________ _______    __________________
   .;@+@';;;;;;'   \______   \\      \  /  _____\__    _______   ______/  |_
    `     '#''@`    |       _//   |   \/   \  ___ |    |_/ __ \ /  ___\   __\
                    |    |   /    |    \    \_\  \|    |\  ___/ \___ \ |  |
                    |____|_  \____|__  /\______  /|____| \___  /____  >|__|
                           \/        \/        \/            \/     \/         )"
      << std::endl;
    }

    //! Print UnitTest header
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
    :;;;;;;.           ____ ___      .__  __ ___________              __
   .;@+@';;;;;;'      |    |   \____ |__|/  |\__    ___/___   _______/  |_
    `     '#''@`      |    |   /    \|  \   __\|    |_/ __ \ /  ___/\   __\
                      |    |  /   |  \  ||  |  |    |\  ___/ \___ \  |  |
                      |______/|___|  /__||__|  |____| \___  >____  > |__|
                                   \/                     \/     \/            )"
      << std::endl;
    }

    //! Print MeshConv header
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

  protected:
    //! Bullets
    const char m_section_bullet = '*';
    const char m_subsection_bullet = '<';
    //! Indents
    const std::string m_section_indent = " ";
    const std::string m_subsection_indent =
      std::operator+(m_section_indent,"  ");
    const std::string m_item_indent = std::operator+(m_subsection_indent,"  ");

    //! Format strings
    using format = boost::format;
    mutable format m_header_fmt = format("%|=80|\n");
    mutable format m_part_fmt = format("\n%|=80|\n");
    mutable format m_section_title_fmt = format("\n%s%c %s:\n");
    mutable format m_section_title_value_fmt = format("\n%s%c %s: %s\n");
    mutable format m_subsection_title_fmt = format("%s%c %s >\n");
    mutable format m_list_item_fmt = format("%s%-30s\n");
    mutable format m_note_fmt = format("%s%-30s\n");
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
};

} // tk::

#endif // Print_h
