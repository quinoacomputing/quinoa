//******************************************************************************
/*!
  \file      src/RNGTest/TestU01.h
  \author    J. Bakosi
  \date      Fri 13 Mar 2015 12:27:44 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     TestU01 statistical test
  \details   TestU01 statistical test
*/
//******************************************************************************
#ifndef TestU01_h
#define TestU01_h

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <charm++.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <testu01.decl.h>
#include <TestU01Props.h>

namespace rngtest {

//! TestU01 statistical test used polymorphically with rngtest::StatTest
template< class TestU01Props >
class TestU01 : public CBase_TestU01< TestU01Props > {

  public:
    using Props = TestU01Props;
    using Proxy = CProxy_TestU01< TestU01Props >;

    //! Constructor
    //! \param[in] props TestU01 properties
    explicit TestU01( const TestU01Props& props ) : m_props( props ) {}

    //! Query and contribute number of results/test, i.e., p-values
    void npval() { m_props.proxy().npval( m_props.npval() ); }

    //! Query and contribute test name(s)
    void names() { m_props.proxy().names( m_props.names() ); }

    //! Run test then evaluate it
    void run() { m_props.proxy().evaluate( m_props.run() ); }

    //! Query and contribute test run time measured in seconds
    void time() { m_props.proxy().time( m_props.time() ); }

  private:
    TestU01Props m_props;               //!< TestU01 test properties
};

} // rngtest::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include <testu01.def.h>
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // TestU01_h
