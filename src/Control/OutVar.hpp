// *****************************************************************************
/*!
  \file      src/Control/OutVar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Types and associated functions to deal with output variables
  \details   Types and associated functions to deal with output variables.
*/
// *****************************************************************************
#ifndef OutVar_h
#define OutVar_h

#include "Keywords.hpp"
#include "Centering.hpp"

namespace tk {
namespace ctr {

//! Output variable
struct OutVar {
  using ncomp_t = kw::ncomp::info::expect::type;

  char var;             //!< Variable name
  ncomp_t field;        //!< Field ID
  Centering centering;  //!< Centering
  std::string name;     //!< Human readable name

  /** @name Pack/Unpack: Serialize OutVar object for Charm++ */
  ///@{
  //! Pack/Unpack serialize member function
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  void pup( PUP::er& p ) {
    p | var;
    p | field;
    p | centering;
    p | name;
  }
  //! Pack/Unpack serialize operator|
  //! \param[in,out] p Charm++'s PUP::er serializer object reference
  //! \param[in,out] v OutVar object reference
  friend void operator|( PUP::er& p, OutVar& v ) { v.pup(p); }
  ///@}

  //! Constructor: initialize all state data
  //! \param[in] v Variable name
  //! \param[in] f Field ID
  //! \param[in] n Human readable name
  explicit OutVar( char v = 0,
                   ncomp_t f = 0,
                   Centering c = Centering::NODE,
                   const std::string& n = {} ) :
    var(v), field(f), centering(c), name(n) {}

  //! \brief Less-than operator for ordering, used by, e.g., std::set::insert
  //! \param[in] outvar OutVar to compare
  //! \return Boolean indicating if outvar is less than 'this'
  bool operator< ( const OutVar& outvar ) const {
    if (var < outvar.var)
      return true;
    else if (var == outvar.var && field < outvar.field)
      return true;
    else if (var == outvar.var && field == outvar.field &&
             centering < outvar.centering)
      return true;
    else if (var == outvar.var && field == outvar.field &&
             centering == outvar.centering && name < outvar.name)
      return true;
    else
      return false;
  }
};

//! \brief Pack/Unpack: Namespace-scope serialize OutVar object for Charm++
//! \param[in,out] p Charm++'s PUP::er serializer object reference
//! \param[in,out] v OutVar object reference
inline void pup( PUP::er& p, OutVar& v ) { v.pup(p); }

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-function"
#endif

//! Operator << for writing OutVar to output streams
//! \param[in,out] os Output stream to write to
//! \param[in] outvar OutVar to write
//! \return Updated output stream
static std::ostream& operator<< ( std::ostream& os, const OutVar& outvar ) {
  if (outvar.name.empty())
    os << outvar.var << outvar.field+1;
  else
    os << outvar.name;
  return os;
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#endif

} // ctr::
} // tk::

#endif // OutVar_h
