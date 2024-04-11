// *****************************************************************************
/*!
  \file      src/Base/Types.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Toolkit-level type definitions
  \details   Toolkit-level type definitions
*/
// *****************************************************************************
#ifndef Types_h
#define Types_h

#include <cstddef>
#include <string>
#include <optional>

namespace tk {

//! Real number type used throughout the whole code.
// TODO Test with single precision and possibly others.
using real = double;

//! uint type used throughout the whole code.
using ncomp_t = std::size_t;

//! Struct for storing info
struct info_t {
  std::string kw;
  std::string sDescr;
  std::string lDescr;
  std::string tDescr;

  // constructor
  info_t(
    std::string k,
    std::string s,
    std::string l,
    std::string t ) :
    kw(k),
    sDescr(s),
    lDescr(l),
    tDescr(t)
    {}

  // Accessors
  std::string keyword() const {return kw;}
  std::string shortDescription() const {return sDescr;}
  std::string longDescription() const {return lDescr;}
  std::string typeDescription() const {return tDescr;}
};

//! Struct for storing keyword in the input deck
struct entry_t {
  info_t info;

  // constructor
  entry_t(
    std::string kw,
    std::string sDescr,
    std::string lDescr,
    std::string tDescr="" ) :
    info(kw, sDescr, lDescr, tDescr)
    {}

  //! Accessor to keyword as std::string
  //! \return Keyword as std::string
  const std::string string() const { return info.keyword(); }

  //! Accessor to required short name of a keyword
  //! \return Name of keyword as std::string
  const std::string name() const { return info.keyword(); }

  //! Accessor to required short description of a keyword
  //! \return Short help as std::string
  const std::string shortDescription() const { return info.shortDescription(); }

  //! Accessor to required long description of a keyword
  //! \return Long help as std::string
  const std::string longDescription() const { return info.longDescription(); }

  //! Alias accessor for keyword
  //! \return Null. This is a punt to fit in the existing HelpFactory
  //!   infrastructure; needs to be removed.
  const std::optional< std::string > alias() const { return std::nullopt; }

  //! Expected type description accessor for keyword
  //! \return Type description for keyword
  const std::optional< std::string > expt() const {
    if (!info.typeDescription().empty())
      return info.typeDescription();
    else
      return std::nullopt;
  }

  //! Expected lower bound accessor for a keyword
  //! \return Null.
  const std::optional< std::string > lower() const { return std::nullopt; }

  //! Expected upper bound accessor for a keyword
  //! \return Null.
  const std::optional< std::string > upper() const { return std::nullopt; }

  //! Expected choices description accessor for a keyword
  //! \return Null.
  const std::optional< std::string > choices() const { return std::nullopt; }

  //! \brief Less-than operator for ordering, used by, e.g., std::set::insert
  //! \param[in] en entry_t to compare
  //! \return Boolean indicating if en is less than 'this'
  bool operator< ( const entry_t& en ) const {
    if (info.kw < en.info.kw)
      return true;
    else
      return false;
  }
};

} // tk::

#endif // Types_h
