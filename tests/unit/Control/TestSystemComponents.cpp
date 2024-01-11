// *****************************************************************************
/*!
  \file      tests/unit/Control/TestSystemComponents.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit tests for Control/SystemComponents
  \details   Unit tests for Control/SystemComponents
*/
// *****************************************************************************

#include <brigand/algorithms/for_each.hpp>
#include <brigand/sequences/list.hpp>

#include "NoWarning/tut.hpp"

#include "TUTConfig.hpp"
#include "SystemComponents.hpp"
#include "TaggedTuple.hpp"
#include "Tags.hpp"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! Tags for unit tests
struct eq1 {};
struct eq2 {};

//! All tests in group inherited from this base
struct SystemComponents_common {
  // Typedef simple set of equations of two types
  using ncomps = tk::ctr::ncomponents< eq1, eq2 >;
  // Typedef vector of all equation tags
  using eqs = brigand::list< eq1, eq2 >;

  // Functor verifying the number of components
  struct testncomp {
    const ncomps& m_host;
    const std::vector< tk::ctr::ncomp_t > m_comps{ 2, 3, 3, 2 };
    std::size_t m_c;
    explicit testncomp( const ncomps& host ) : m_host( host ), m_c( 0 ) {}
    template< typename U > void operator()( brigand::type_<U> ) {
      for (const auto& c : m_host.get< U >())
        ensure_equals( "number of components for c=" + std::to_string(c),
                       c, m_comps[m_c++] );
    }
  };
};

//! Test group shortcuts
using SystemComponents_group =
  test_group< SystemComponents_common, MAX_TESTS_IN_GROUP >;
using SystemComponents_object = SystemComponents_group::object;

//! Define test group
static SystemComponents_group SystemComponents( "Control/SystemComponents" );

//! Test definitions for group

//! Test that number of components are correct
template<> template<>
void SystemComponents_object::test< 1 >() {
  set_test_name( "number of components" );

  // Instantiate number of components object
  ncomps nc;

  // Add a couple of systems of eq1 equations with components 2 and 3, resp.
  nc.get< eq1 >().push_back( 2 );
  nc.get< eq1 >().push_back( 3 );

  // Add a couple of systems of eq2 equations with components 3 and 2, resp.
  nc.get< eq2 >().push_back( 3 );
  nc.get< eq2 >().push_back( 2 );

  // Test number of components of all equations
  brigand::for_each< eqs >( testncomp( nc ) );
}

//! Test the total number of components are correct
template<> template<>
void SystemComponents_object::test< 3 >() {
  set_test_name( "total number of components" );

  // Instantiate number of components object
  ncomps nc;

  // Add a couple of systems of eq1 equations with components 2 and 3, resp.
  nc.get< eq1 >().push_back( 2 );
  nc.get< eq1 >().push_back( 3 );

  // Add a couple of systems of eq2 equations with components 3 and 2, resp.
  nc.get< eq2 >().push_back( 3 );
  nc.get< eq2 >().push_back( 2 );

  // Test offsets of all equations
  ensure_equals( "total number of components", nc.nprop(0), 5 );
  ensure_equals( "total number of components", nc.nprop(1), 5 );
}

} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
