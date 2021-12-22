// *****************************************************************************
/*!
  \file      src/Control/Keyword.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Generic definition of a keyword
  \details   Generic definition of all keywords - both command-line arguments
    and control file keywords.
*/
// *****************************************************************************
#ifndef Keyword_h
#define Keyword_h

#include <optional>

#include "NoWarning/set.hpp"
#include "NoWarning/pegtl.hpp"

#include "Has.hpp"
#include "Escaper.hpp"

namespace tk {

//! Helper to declare a set of command line keywords
//! \details This ensures that a compile-error is generated if there is no alias
//!    defined for the keyword, and also if the aliases are non-unique.
template< typename... T >
class cmd_keywords {
  public:
    using set = brigand::set< T... >;
  private:
    template< typename K > using alias = typename K::info::alias::type;
    using aliases = brigand::set< alias<T>... >;
};

} // tk::

namespace kw {

using namespace tao;

//! \brief Keyword alias helper
//! \details This struct is used to define both a type and a value for a keyword
//!   alias, which is a single character. Used for command-line arguments, e.g.,
//!   --help, -h, where 'h' is the alias for keyword 'help'.
//! \see Control/Keywords.h
template< int Char >
struct Alias {
  using type = pegtl::one< Char >;
  static const int value = Char;
};

//! \brief Generic definition of a keyword
//! \details A keyword is a struct that collects the information that makes up a
//!    keyword. The requirement on the first template argument, Info, is that it
//!    must define the name(), shortDescription(,) and longDescription() member
//!    functions returning compile-time (static) std::strings. The
//!    shortDescription() member function is used to return a short description
//!    of what the keyword is used for, while the longDescription() member
//!    function is used for a longer, e.g., a paragraph-long, description on
//!    what the keyword can be used for and how it can and should be used. The
//!    last template parameter is a pegtl string, a list of character constants,
//!    specifying the case-sensitive characters that make up the keyword, which
//!    is then matched by the parser. The keyword must be at least one character
//!    long, but otherwise its length is only limited by the compiler's
//!    recursion handling capability of variadic templates. While the name(),
//!    shortDescription() and longDescription() member functions of Info are
//!    required, there are also optional ones, such as
//!    Info::exptect::description(), which, if defined, must also be static and
//!    must return a std::string, describing the type the particular keyword
//!    expects during parsing. This is optional since not every keyword expects
//!    a value (or values) of a particular type. For example, the keyword 'end'
//!    is simply used to close a block in the input file, and what follows does
//!    not have a relationship to the keyword. A counterexample is is 'title',
//!    which expects a double-quoted string immediately after the keyword
//!    title'.
//! \see For example client-code and more detailed documentation on the possible
//!    fields, see Control/Keywords.h.
template< typename Info, typename > struct keyword;
template< typename Info, char... Chars >
struct keyword< Info, pegtl::string< Chars... > > {

  //! Accessor to keyword as pegtl::string
  using pegtl_string = pegtl::string< Chars... >;

  //! Accessor to keyword as std::string
  //! \return Keyword as std::string
  static std::string string() { return kw::escaper< Chars... >::result(); }

  //! Accessor to required short name of a keyword
  //! \return Name of keyword as std::string
  static std::string name() { return Info::name(); }

  //! Accessor to required short description of a keyword
  //! \return Short help as std::string
  static std::string shortDescription() { return Info::shortDescription(); }

  //! Accessor to required long description of a keyword
  //! \return Long help as std::string
  static std::string longDescription() { return Info::longDescription(); }

  //! Bring template argument 'Info' to scope as 'info'
  //! \details This is used to access, e.g., Info::alias, etc., if exist.
  //! \see tk::grm::alias
  //! \see tk::grm::readcmd
  using info = Info;

  //! Alias accessor for keyword
  //! \return An initialized (or uninitialized) std::optional< std::string >
  //! \details Though an alias is only a single character, it returns it as
  //!   std::string since pegtl::string returns std::string.
  template< typename T = Info >
  static std::optional< std::string > alias() {
    if constexpr( tk::HasTypedef_alias_v< T > )
      return std::string( 1, static_cast<char>( Info::alias::value ) );
    else
      return std::nullopt;
  }

  //! Expected type description accessor for keyword
  //! \return An initialized (or uninitialized) std::optional< std::string >
  template< typename T = Info >
  static std::optional< std::string > expt() {
    if constexpr( tk::HasFunction_expect_description_v< T > )
      return Info::expect::description();
    else
      return std::nullopt;
  }

  //! Expected choices description accessor for a keyword
  //! \return An initialized (or uninitialized) std::optional< std::string >
  template< typename T = Info >
  static std::optional< std::string > choices() {
    if constexpr( tk::HasFunction_expect_choices_v< T > )
      return Info::expect::choices();
    else
      return std::nullopt;
  }

  //! Expected lower bound accessor for a keyword
  //! \return An initialized (or uninitialized) std::optional< std::string >
  template< typename T = Info >
  static std::optional< std::string > lower() {
    if constexpr( tk::HasVar_expect_lower_v< T > )
      return std::to_string( Info::expect::lower );
    else
      return std::nullopt;
  }

  //! Expected upper bound accessor for a keyword
  //! \return An initialized (or uninitialized) std::optional< std::string >
  template< typename T = Info >
  static std::optional< std::string > upper() {
    if constexpr( tk::HasVar_expect_upper_v< T > )
      return std::to_string( Info::expect::upper );
    else
      return std::nullopt;
  }
};

} // kw::

#endif // Keyword_h
