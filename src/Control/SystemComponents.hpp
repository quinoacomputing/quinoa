// *****************************************************************************
/*!
  \file      src/Control/SystemComponents.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Operations on numbers of scalar components of systems of equations
  \details   Operations on numbers of scalar components of systems of equations,
    e.g. multiple equation sets of a physical model or a set of (stochastic
    ordinary or partial differential) equations of different kinds.

    _Problem:_ We are given a type that is a tk::tuple::tagged_tuple that
    contains an arbitrary number of std::vectors of integers. The number of
    vectors are fixed at compile-time (accessed via tags) but their (differing)
    length is only known at run-time (after parsing user input). What we need
    are functions that operate on this data structure and return, e.g., the
    total number of components in the whole system or the offset in the whole
    data structure for a given tag. The functions should thus be able to operate
    on a list of types, i.e., a double for-loop over all tags and associated
    vectors - one at compile-time and the other one at run-time.

    _Example:_ An example is to define storage for systems of stochastic
    differential equations as in walker::ctr::ncomps. In
    [walker](walker.html) various types of stochastic differential
    equations are implemented and numerically integrated in time in order to
    compute and learn about their statistical behavior. Since different types of
    systems can be parameterized differently, there can be several systems of
    the same type (the number is user-defined so only known at run-time). For
    example, it is possible to numerically integrate two systems of Dirichlet
    equations, see DiffEq/Dirichlet.h, one with 4, the other one with 5 scalar
    variables, parameterized differently, and estimate their arbitrary coupled
    statistics and joint PDFs. This requires that each type of equation system
    has a vector of integers storing the number of scalar variables.

    _Solution:_ Looping through a list of types is done using the [Brigand]
    (https://github.com/edouarda/brigand)'s [MetaProgramming Library].
    Such operations on types happen at compile-time, i.e., the code runs inside
    the compiler and only its result gets compiled into code to be run at
    run-time. Advantages are abstraction and generic code that is independent
    of the size and order of the tags in the tuple (and the associated vectors).
    Though it is not of primary concern for this data structure, this solution
    also generates efficient code, since part of the algorithm (for computing
    the offset, the total number of components, etc.) runs inside the compiler.
    All of this is kept generic, so it does not need to be changed when a new
    pair of tag and system of equations are added to the underlying tuple.
    Thus, _the main advantage is that changing the order of the vectors in the
    tuple or adding new ones at arbitrary locations will not require change to
    client-code._

    _An alternative approach:_ Of course this could have been simply done with a
    purely run-time vector of vector of integers. However, that would not have
    allowed a tag-based access to the different systems of equations. (Sure, one
    can define an enum for equation indexing, but that just adds code to what is
    already need to be changed if a change happens to the underlying data
    structure.) As a consequence, a component in a vector of vector would have
    had to be accessed as something like, ncomps[2][3], which is more
    error-prone to read than ncomps< tag::dirichlet >( 3 ). Furthermore, that
    design would have resulted in double-loops at run-time, instead of letting
    the compiler loop over the types at compile-time and do the run-time loops
    only for what is only known at run-time. Additionally, and most importantly,
    _any time a new system is introduced (or the order is changed), all
    client-code would have to be changed._

    For example client-code, see walker::Dirichlet::Dirichlet() in
    DiffEq/Dirichlet.h, or walker::Integrator::Integrator in Walker/Integrator.h.
*/
// *****************************************************************************
#ifndef SystemComponents_h
#define SystemComponents_h

#include <vector>
#include <algorithm>
#include <numeric>
#include <string>

#include <brigand/sequences/list.hpp>
#include <brigand/algorithms/wrap.hpp>

#include "NoWarning/flatten.hpp"
#include "NoWarning/transform.hpp"

#include "TaggedTuple.hpp"
#include "StatCtr.hpp"
#include "Keywords.hpp"
#include "Tags.hpp"

namespace tk {
//! Toolkit control, general purpose user input to internal data transfer
namespace ctr {

//! Inherit type of number of components from keyword 'ncomp'
using ncomp_t = kw::ncomp::info::expect::type;

//! \brief Map associating number of scalar components to dependent variables
//!   for systems
//! \details This map associates the number of properties (scalar components)
//!   of systems of differential equations for all scalar components of a
//!   system of systems.
//! \note We use a case-insensitive character comparison functor to be
//!   consistent with OffsetMap.
using NcompMap = std::map< char, ncomp_t, CaseInsensitiveCharLess >;

//! Helper for converting a brigand::list to a tagged_tuple
template< typename... T >
using tagged_tuple_wrapper = typename tk::TaggedTuple< brigand::list<T...> >;

//! Helper for converting a brigand::list to a tagged_tuple
template< typename L >
using as_tagged_tuple = brigand::wrap< L, tagged_tuple_wrapper >;

//! Number of components storage as a vector for a system of equations
//! \details This is only a helper class, defining a type 'type' for
//!    brigand::apply, so it can be used for defining a base for ncomponents
struct ComponentVector : public std::vector< ncomp_t > {
  using type = std::vector< ncomp_t >;
};

//! \brief Number of components storage
//! \details All this trickery with template meta-programming allows the code
//!   below to be generic. As a result, adding a new component requires adding a
//!   single line (a tag and its type) to the already existing list, see
//!   typedefs 'ncomps'. The member functions, doing initialization, computing
//!   the number of total components, the offset for a given tag, and computing
//!   the offset map, need no change -- even if the order of the number of
//!   components change.
template< typename... Tags >
class ncomponents : public
  // tk::tuple::tagged_tuple< tag1, vec1, tag2, vec2, ... >
  as_tagged_tuple< brigand::flatten< brigand::transform< brigand::list<Tags...>,
    brigand::bind< brigand::list, brigand::_1, ComponentVector > > > > {

  private:
    //! Access vector of number of components of an eq system as const-ref
    template< typename Eq >
    constexpr const auto& vec() const { return this->template get< Eq >(); }

    //! Access vector of number of components of an eq system as reference
    template< typename Eq >
    constexpr auto& vec() { return this->template get< Eq >(); }

  public:
    //! Default constructor: set defaults to zero for all number of components
    ncomponents() {
      ( ... , std::fill( begin(vec<Tags>()), end(vec<Tags>()), 0 ) );
    }

    //! Compute total number of components in the systems of systems configured
    //! \return Total number of components
    ncomp_t nprop() const noexcept {
      return (... + std::accumulate(begin(vec<Tags>()), end(vec<Tags>()), 0u));
    }

    //! \brief Compute map of number of properties (scalar components)
    //!   associated to dependent variables
    //! \param[in] d Input deck to operate on
    //! \return Map of number of properties associated to dependent variables
    template< class InputDeck >
    NcompMap ncompmap( const InputDeck& d ) const {
      NcompMap map;
      ( ... , [&](){
        const auto& depvar = d.template get< tag::param, Tags, tag::depvar >();
        const auto& ncomps = d.template get< tag::component >();
        const auto& ncvec = ncomps.template get<Tags>();
        Assert( ncvec.size() == depvar.size(), "ncompsize != depvarsize" );
        ncomp_t c = 0;
        for (auto v : depvar) map[ v ] = ncvec[c++]; }() );
      return map;
    }

    //! \brief Return vector of dependent variables + component id for all
    //!   equations configured
    //! \param[in] deck Input deck to operate on
    //! \return Vector of dependent variables + comopnent id for all equations
    //!   configured. The length of this vector equals the total number of
    //!   components configured, see nprop(), containing the depvar + the
    //!   component index relative to the given equation. E.g., c1, c2, u1, u2,
    //!   u3, u4, u5.
    template< class InputDeck >
    std::vector< std::string > depvar( const InputDeck& deck ) const {
      std::vector< std::string > d;
      ( ..., [&](){
        const auto& dveq = deck.template get< tag::param, Tags, tag::depvar >();
        const auto& nceq = deck.template get< tag::component, Tags >();
        Assert( dveq.size() == nceq.size(), "Size mismatch" );
        std::size_t e = 0;
        for (auto v : dveq) {
          for (std::size_t c=0; c<nceq[e]; ++c )
            d.push_back( v + std::to_string(c+1) );
          ++e;
        } }() );
      return d;
    }
};

} // ctr::
} // tk::

#endif // SystemComponents_h
