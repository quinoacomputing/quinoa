// *****************************************************************************
/*!
  \file      src/DiffEq/DiffEqFactory.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Differential equations factory
  \details   This file declares the type for a differential equations factory.
*/
// *****************************************************************************
#ifndef DiffEqFactory_h
#define DiffEqFactory_h

#include <map>
#include <functional>

#include "DiffEq.h"
#include "Factory.h"
#include "SystemComponents.h"
#include "Walker/Options/DiffEq.h"

namespace walker {

//! Differential equation factory: keys associated to their constructors
using DiffEqFactory =
  std::map< ctr::DiffEqKey,
            std::function< DiffEq(const tk::ctr::ncomp_type&) > >;

//! \brief Function object for registering a differential equation into the
//!   differential equation factory
//! \details This functor is repeatedly called by brigand's cartesian_product,
//!   sweeping all combinations of the differential equation policies. The
//!   purpose of template template is to simplify client code as that will
//!   not have to specify the template arguments of the template argument
//!   (the policies of Eq), since we can figure it out here. See also
//!   http://stackoverflow.com/a/214900
template< template< class, class > class Eq >
struct registerDiffEq {
  //! Need to store the reference to factory we are registering into
  DiffEqFactory& factory;
  //! Need to store which differential equation we are registering
  const ctr::DiffEqType type;
  //! Constructor, also count number of unique equation types registered
  explicit registerDiffEq( DiffEqFactory& f,
                           ctr::DiffEqType t,
                           std::set< ctr::DiffEqType >& eqTypes ) :
    factory( f ), type( t ) { eqTypes.insert( t ); }
  //! \brief Function call operator called with tk::cartesian_product for
  //!   each unique sequence of policy combinations
  template< typename U > void operator()( brigand::type_<U> ) {
    // Get Initialization policy: first type of brigand::list U
    using InitPolicy = typename brigand::front< U >;
    // Get coefficients policy: last type of brigand::list U
    using CoeffPolicy = typename brigand::back< U >;
    // Build differential equation key
    ctr::DiffEqKey key{ type, InitPolicy::type(), CoeffPolicy::type() };
    // Register equation (with policies given by brigand::list U) into
    // factory
    tk::recordModelLate< DiffEq, Eq< InitPolicy, CoeffPolicy > >
                       ( factory, key, static_cast<tk::ctr::ncomp_type>(0) );
  }
};

//! \brief Convert and return values from vector as string
//! \param[in] v Vector whose components to return as a string
//! \return Concatenated string of values read from a vector
template< typename V >
std::string parameters( const V& v ) {
  std::stringstream s;
  s << "{ ";
  for (auto p : v) s << p << ' ';
  s << "}";
  return s.str();
}

//! \brief Return names of options (tk::Toggle) from vector as a string
//! \param[in] opt Option instance (inheriting from tk::Toggle)
//! \param[in] v Option vector whose names of components to return
//! \return Concatenated string of option names read from option vector
template< class Option, class OptTypeVec >
std::string options( const Option& opt, const OptTypeVec& v ) {
  std::stringstream s;
  s << "{ ";
  for (auto o : v) s << opt.name(o) << ' ';
  s << "}";
  return s.str();
}

//! \brief Insert spike information (used to specify delta PDFs) into info
//!   vector
//! \param[in,out] nfo Info vector of string-pairs to insert to
//! \param[in] spike Vector of vectors specifying spike info
template< typename Info, typename VV >
void spikes( Info& nfo, const VV& spike ) {
  std::size_t i = 0;
  for (const auto& s : spike)
    nfo.emplace_back( "delta spikes [comp" + std::to_string(++i) + ":" +
                        std::to_string( s.size()/2 ) + "]",
                      parameters( s ) );
}

//! \brief Insert betapdf information (used to specify beta PDFs) into info
//!   vector
//! \param[in,out] nfo Info vector of string-pairs to insert to
//! \param[in] betapdf Vector of vectors specifying betapdf info
template< typename Info, typename VV >
void betapdfs( Info& nfo, const VV& betapdf ) {
  std::size_t i = 0;
  for (const auto& s : betapdf)
    nfo.emplace_back( "beta pds [comp" + std::to_string(++i) + "]",
                      parameters( s ) );
}

} // walker::

#endif // DiffEqFactory_h
