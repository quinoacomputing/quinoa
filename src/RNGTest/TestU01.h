//******************************************************************************
/*!
  \file      src/RNGTest/TestU01.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 05:25:50 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 statistical tests
  \details   TestU01 statistical tests
*/
//******************************************************************************
#ifndef TestU01_h
#define TestU01_h

#include <charm++.h>

#include <TestU01Props.h>
#include <testu01.decl.h>
#include <Msg.h>

namespace rngtest {

//! TestU01 statistical test used polymorphically with StatTest
template< class TestU01Props >
class TestU01 : public CBase_TestU01< TestU01Props > {

  public:
    using Props = TestU01Props;
    using Proxy = CProxy_TestU01< TestU01Props >;

    //! Constructor, CHARM: should move
    explicit TestU01( const TestU01Props& props ) : m_props( props ) {}

    //! Run test
    void run( std::size_t id ) {
      std::cout << "TestU01 run\n";
      //m_props.run();
      m_props.proxy().evaluate( id );
    }

    //! Contribute number of results/test, i.e., p-values
    void npval( CkFuture f ) {
      using Msg = tk::Msg< std::size_t >;
      Msg* n = new Msg( m_props.npval() );
      CkSendToFuture( f, n );
    }

    //! Contribute test name
    void name( CkFuture f ) const {
      using Msg = tk::StringsMsg;
      Msg* n = new Msg( m_props.name() );
      CkSendToFuture( f, n );
    }

//     //! RNG enum accessor
//     const tk::ctr::RNGType& rng() const { return m_id; }
// 
//     //! RNG id accessor
//     std::size_t id() const { return 0; }
// 
//     //! Query whether test is failed
//     bool fail( std::size_t p ) const {
//       if ( (m_pvals[p] <= gofw_Suspectp) || (m_pvals[p] >= 1.0-gofw_Suspectp) )
//         return true;
//       else
//         return false;
//     }
// 
//     //! Return number of failed tests
//     std::size_t nfail() const {
//       std::size_t f = 0;
//       for (std::size_t p = 0; p<m_npval; ++p) if (fail(p)) ++f;
//       return f;
//     }
// 
//     //! p-value accessor
//     double pval( std::size_t p ) const { return m_pvals[p]; }
// 
//     //! Return human-readable p-value (based on TestU01::bbattery.c::WritePval)
//     std::string pvalstr( std::size_t p ) const {
//       std::stringstream ss;
//       double val = m_pvals[p];
//       if (val < gofw_Suspectp) {
//         if ((val >= 0.01) && (val <= 0.99))
//           ss << val;
//         else if (val < gofw_Epsilonp)
//           ss << "eps";
//         else if (val < 0.01)
//           ss << val;
//         else if (val >= 1.0 - gofw_Epsilonp1)
//           ss << "1 - eps1";
//         else if (val < 1.0 - 1.0e-4)
//           ss << val;
//         else
//           ss << 1.0 - val;
//       } else if (val > 1.0 - gofw_Suspectp) {
//         if (val >= 1.0 - gofw_Epsilonp1)
//           ss << "1 - eps1";
//         else if (val >= 1.0 - 1.0e-4)
//           ss << "1 - " << 1.0 - val;
//         else
//           ss << val;
//       }
//       return ss.str();
//     }

//     //! Copy assignment
//     TestU01& operator=( const TestU01& x) {
//       m_props = x.m_props;
//       return *this;
//     }
//     //! Copy constructor: in terms of copy assignment
//     TestU01( const TestU01& x ) { operator=(x); }
//     //! Move assignment
//     TestU01& operator=( TestU01&& ) noexcept = default;
//     //! Move constructor
//     TestU01( TestU01&& ) noexcept = default;

  private:
    TestU01Props m_props;               //!< TestU01 test properties
};

} // rngtest::

#define CK_TEMPLATES_ONLY
#include <testu01.def.h>
#undef CK_TEMPLATES_ONLY

#endif // TestU01_h
