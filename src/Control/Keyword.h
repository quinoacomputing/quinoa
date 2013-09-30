//******************************************************************************
/*!
  \file      src/Control/Keyword.h
  \author    J. Bakosi
  \date      Mon Sep 30 08:50:46 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Basic keywords recognized by all parsers
  \details   Basic keywords recognized by all parsers
*/
//******************************************************************************
#ifndef Keyword_h
#define Keyword_h

#include <pegtl.hh>

namespace quinoa {
namespace grm {
namespace kw {

//! A keyword is struct that collects the information that makes up a keyword.
//! The last template parameter is a list of integers, specifying the
//! case-sensitive characters of the keyword. The keyword must be at least one
//! character long, but otherwise its length is only limited by the compiler's
//! recursion handling of variadic templates. The namespace pegtl::ascii
//! contains some predefined helper definitions of upper and lower-case letters,
//! for all other characters, just quote it according to the C-language rules,
//! e.g., '_' is the underscore.
template< typename Info, int Char, int... Chars >
struct keyword {

  //! Accessor to keyword as pegtl::string
  using pegtl_string = pegtl::string<Char, Chars...>;

  //! Accessor to keyword as std::tring
  std::string string() const {
    return std::move(std::string( (sizeof...(Chars)) ?
                                  (pegtl::escaper<Char, Chars...>::result()) :
                                  (pegtl::escape(Char)) ));
  }

  //! Accessor to name (more human readable but still short description)
  const char* name() const { return Info::name(); }

  //! Accessor to help (a few lines' worth description of the keyword
  const char* help() const { return Info::help(); }
};

} // kw::
} // grm::
} // quinoa::

#endif // Keyword_h
