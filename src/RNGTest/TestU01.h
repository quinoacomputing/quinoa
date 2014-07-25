//******************************************************************************
/*!
  \file      src/RNGTest/TestU01.h
  \author    J. Bakosi
  \date      Thu 24 Jul 2014 03:13:52 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     TestU01 statistical test
  \details   TestU01 statistical test
*/
//******************************************************************************
#ifndef TestU01_h
#define TestU01_h

#include <charm++.h>

#include <testu01.decl.h>
#include <TestU01Props.h>

namespace rngtest {

//! TestU01 statistical test used polymorphically with StatTest
template< class TestU01Props >
class TestU01 : public CBase_TestU01< TestU01Props > {

  public:
    using Props = TestU01Props;
    using Proxy = CProxy_TestU01< TestU01Props >;

    //! Constructor
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

#define CK_TEMPLATES_ONLY
#include <testu01.def.h>
#undef CK_TEMPLATES_ONLY

#endif // TestU01_h
