// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Control/TestSystemComponents.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unit tests for Control/SystemComponents
  \details   Unit tests for Control/SystemComponents
*/
// *****************************************************************************

#include <brigand/algorithms/for_each.hpp>
#include <brigand/sequences/list.hpp>

#include "NoWarning/tut.h"

#include "TUTConfig.h"
#include "SystemComponents.h"
#include "TaggedTuple.h"
#include "Tags.h"
#include "Control.h"

#ifndef DOXYGEN_GENERATING_OUTPUT

namespace tut {

//! Tags for unit tests
struct eq1 {};
struct eq2 {};

//! All tests in group inherited from this base
struct SystemComponents_common {
  // Typedef simple set of equations of two types
  using ncomps = tk::ctr::ncomponents<
    eq1, std::vector< tk::ctr::ncomp_type >,
    eq2, std::vector< tk::ctr::ncomp_type >
  >;
  // Typedef vector of all equation tags
  using eqs = brigand::list< eq1, eq2 >;

  // Functor verifying the number of components
  struct testncomp {
    const ncomps& m_host;
    const std::vector< tk::ctr::ncomp_type > m_comps{ 2, 3, 3, 2 };
    std::size_t m_c;
    testncomp( const ncomps& host ) : m_host( host ), m_c( 0 ) {}
    template< typename U > void operator()( brigand::type_<U> ) {
      for (const auto& c : m_host.get< U >())
        ensure_equals( "number of components for c=" + std::to_string(c),
                       c, m_comps[m_c++] );
    }
  };

  // Functor verifying the offset
  struct testoffset {
    const ncomps& m_host;
    const std::vector< tk::ctr::ncomp_type > m_offs{ 0, 2, 5, 8 };
    std::size_t m_c;
    testoffset( const ncomps& host ) : m_host( host ), m_c( 0 ) {}
    template< typename U > void operator()( brigand::type_<U> ) {
      for (std::size_t c=0; c<m_host.get< U >().size(); ++c)
        ensure_equals( "offset for c=" + std::to_string(c),
                       m_host.offset< U >(c), m_offs[m_c++] );
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

//! Test that number of components are correct
template<> template<>
void SystemComponents_object::test< 2 >() {
  set_test_name( "offsets" );

  // Instantiate number of components object
  ncomps nc;

  // Add a couple of systems of eq1 equations with components 2 and 3, resp.
  nc.get< eq1 >().push_back( 2 );
  nc.get< eq1 >().push_back( 3 );

  // Add a couple of systems of eq2 equations with components 3 and 2, resp.
  nc.get< eq2 >().push_back( 3 );
  nc.get< eq2 >().push_back( 2 );

  // Test offsets of all equations
  brigand::for_each< eqs >( testoffset( nc ) );
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
  ensure_equals( "total number of components", nc.nprop(), 10 );
}

//! Test that offsetmap builds a the correct map of offsets
template<> template<>
void SystemComponents_object::test< 4 >() {
  set_test_name( "offsetmap" );

  // Instantiate number of components object
  ncomps nc;

  // Add a couple of systems of eq1 equations with components 2 and 3, resp.
  nc.get< eq1 >().push_back( 2 );
  nc.get< eq1 >().push_back( 3 );

  // Add a couple of systems of eq2 equations with components 3 and 2, resp.
  nc.get< eq2 >().push_back( 3 );
  nc.get< eq2 >().push_back( 2 );

  // Dependent variables (as the only parameters) for two equation systems
  using eq1_parameters = tk::tuple::tagged_tuple<
    tag::depvar, std::vector< char >
  >;
  using eq2_parameters = tk::tuple::tagged_tuple<
    tag::depvar, std::vector< char >
  >;

  // Parameters for two equation systems
  using parameters = tk::tuple::tagged_tuple<
    eq1, eq1_parameters,
    eq2, eq2_parameters
  >;

  // Crate mock input deck with two systems of eqations with dependent variables
  tk::Control< tag::component, ncomps,
               tag::param, parameters > deck;

  // Assign ncomps to deck
  deck.get< tag::component >() = std::move( nc );

  // Assign character codes for dependent variables
  deck.get< tag::param, eq1, tag::depvar >() = { 'a', 'b' };
  deck.get< tag::param, eq2, tag::depvar >() = { 'c', 'd' };

  // Test offsetmap in deck
  ensure( "offsetmap",
          deck.get< tag::component >().offsetmap( deck ) ==
            tk::ctr::OffsetMap{ {'a',0}, {'b',2}, {'c',5}, {'d',8} } );
}

//! Test that ncompmap builds a the correct map of number of components
template<> template<>
void SystemComponents_object::test< 5 >() {
  set_test_name( "ncompmap" );

  // Instantiate number of components object
  ncomps nc;

  // Add a couple of systems of eq1 equations with components 2 and 3, resp.
  nc.get< eq1 >().push_back( 2 );
  nc.get< eq1 >().push_back( 5 );

  // Add a couple of systems of eq2 equations with components 3 and 2, resp.
  nc.get< eq2 >().push_back( 3 );
  nc.get< eq2 >().push_back( 8 );

  // Dependent variables (as the only parameters) for two equation systems
  using eq1_parameters = tk::tuple::tagged_tuple<
    tag::depvar, std::vector< char >
  >;
  using eq2_parameters = tk::tuple::tagged_tuple<
    tag::depvar, std::vector< char >
  >;

  // Parameters for two equation systems
  using parameters = tk::tuple::tagged_tuple<
    eq1, eq1_parameters,
    eq2, eq2_parameters
  >;

  // Crate mock input deck with two systems of eqations with dependent variables
  tk::Control< tag::component, ncomps,
               tag::param, parameters > deck;

  // Assign ncomps to deck
  deck.get< tag::component >() = std::move( nc );

  // Assign character codes for dependent variables
  deck.get< tag::param, eq1, tag::depvar >() = { 'a', 'b' };
  deck.get< tag::param, eq2, tag::depvar >() = { 'c', 'd' };

  // Test ncompmap in deck
  ensure( "ncompmap",
          deck.get< tag::component >().ncompmap( deck ) ==
            tk::ctr::OffsetMap{ {'a',2}, {'b',5}, {'c',3}, {'d',8} } );
}


} // tut::

#endif  // DOXYGEN_GENERATING_OUTPUT
