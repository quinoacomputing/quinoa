//******************************************************************************
/*!
  \file      src/Control/Keyword.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 02:19:40 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Basic keywords recognized by all parsers
  \details   Basic keywords recognized by all parsers
*/
//******************************************************************************
#ifndef Keyword_h
#define Keyword_h

#ifndef __INTEL_COMPILER
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wshadow"
#endif
#include <pegtl.hh>
#ifndef __INTEL_COMPILER
  #pragma GCC diagnostic pop
#endif

namespace kw {

//! A keyword is a struct that collects the information that makes up a keyword.
//! The last template parameter is a list of integers, specifying the
//! case-sensitive characters of the keyword. The keyword must be at least one
//! character long, but otherwise its length is only limited by the compiler's
//! recursion handling of variadic templates.
template< typename Info, int Char, int... Chars >
struct keyword {

  //! Accessor to keyword as pegtl::string
  using pegtl_string = pegtl::string<Char, Chars...>;

  //! Accessor to keyword as std::tring
  std::string string() const {
    return std::string( (sizeof...(Chars)) ?
                        (pegtl::escaper<Char, Chars...>::result()) :
                        (pegtl::escape(Char)) );
  }

  //! Accessor to name (more human readable but still short description)
  const char* name() const { return Info::name(); }

  //! Accessor to help (a few lines' worth description of the keyword
  const char* help() const { return Info::help(); }
};

//! A cmdline_keyword is a keyword that adds a character alias to a keyword.
template< typename Info, int Alias, int Char, int... Chars >
struct cmdline_keyword : keyword< Info, Char, Chars... > {

  //! Accessor to keyword alias character as pegtl::string
  using pegtl_alias = pegtl::one< Alias >;

  //! Accessor to keyword alias character as st::string
  std::string alias() const { return pegtl::escape( Alias ); }
};

// This will go away once all the keywords are documented
struct undefined_info {
  static const char* name() { return "undefined"; }
  static const char* help() { return "Undefined."; }
};

} // kw::

#endif // Keyword_h
