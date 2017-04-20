// *****************************************************************************
/*!
  \file      src/UnitTest/tests/Control/TestSystemComponents.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Unit tests for Control/SystemComponents
  \details   Unit tests for Control/SystemComponents
*/
// *****************************************************************************
#ifndef test_SystemComponents_h
#define test_SystemComponents_h

#include "NoWarning/tut.h"

#include "SystemComponents.h"

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
  using eqs = boost::mpl::vector< eq1, eq2 >;

  // Functor verifying the number of components
  struct testncomp {
    const ncomps& m_host;
    const std::vector< tk::ctr::ncomp_type > m_comps{ 2, 3, 3, 2 };
    std::size_t m_c;
    testncomp( const ncomps& host ) : m_host( host ), m_c( 0 ) {}
    template< typename U > void operator()( U ) {
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
    template< typename U > void operator()( U ) {
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
//! \author J. Bakosi
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
  boost::mpl::for_each< eqs >( testncomp( nc ) );
}

//! Test that number of components are correct
//! \author J. Bakosi
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
  boost::mpl::for_each< eqs >( testoffset( nc ) );
}

//! Test the total number of components are correct
//! \author J. Bakosi
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

//! Test that offsetmap builds a linear map of offsets
//! \author J. Bakosi
template<> template<>
void SystemComponents_object::test< 4 >() {
  set_test_name( "linear offsetmap" );

  // Instantiate number of components object
  ncomps nc;

  // Add a couple of systems of eq1 equations with components 2 and 3, resp.
  nc.get< eq1 >().push_back( 2 );
  nc.get< eq1 >().push_back( 3 );

  // Add a couple of systems of eq2 equations with components 3 and 2, resp.
  nc.get< eq2 >().push_back( 3 );
  nc.get< eq2 >().push_back( 2 );

  // Create vector of vectors of dependent variables denoted by characters
  std::vector< std::vector< char > > depvars;
  depvars.push_back( { 'a', 'b' } );    // for the two eq1 systems
  depvars.push_back( { 'c', 'd' } );    // for the two eq2 systems

  // Test if offsetmap is linear
  ensure( "linear offsetmap",
          nc.offsetmap( depvars ) ==
            tk::ctr::OffsetMap{ {'a',0}, {'b',1}, {'c',2}, {'d',3} } );
}

} // tut::

#endif // test_SystemComponents_h
