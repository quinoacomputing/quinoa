// *****************************************************************************
/*!
  \file      src/RNGTest/TestU01Props.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     TestU01 statistical test properties class
  \details   This file defines a generic TestU01 statistical test properties
    class, used to initialize, interface, and evaluate TestU01 RNG statistical
    tests by the TestU01 library.
*/
// *****************************************************************************
#ifndef TestU01Props_h
#define TestU01Props_h

#include <vector>
#include <functional>

#include "Macro.h"
#include "RNG.h"
#include "Timer.h"
#include "Options/RNG.h"
#include "TestU01Util.h"
#include "TestStack.h"
#include "TestU01Wrappers.h"
#include "QuinoaConfig.h"

namespace rngtest {

extern TestStack g_testStack;

//! \brief TestU01 properties used to initialize TestU01 tests.
//! \details This class is used to abstract away the initialization for TestU01
//!   statistical tests. Needed because of the move semantics and variadic
//!   templates, since Charm++ chares do not support these yet. TestU01Props is
//!   therefore not a chare, but TestU01, initialized with a TestU01Props
//!   object, is. Note that TestU01Props still needs to be migratable, i.e.,
//!   defines the pack/unpack operator, as it is an argument to chare TestU01's
//!   constructor.
//! \author J. Bakosi
template< class Test,                 //!< Statistical test tag, struct{}
          class Proxy,                //!< Host proxy type
          class Result,               //!< Results type
          Result* (*Creator)(void),   //!< Results creator function pointer type
          void (*Deleter)(Result *),  //!< Results deleter function pointer type
          typename... Ts >            //!< Extra test runner argument types
class TestU01Props {

  public:
    //! Test extra arguments type
    using Xargs = std::tuple< Ts... >;
    //! Test runner function pointer type
    using RunFn = std::vector<double> (*)( unif01_Gen*, Result*, const Xargs& );
    //! TestU01 results type with a custom deleter by TestU01
    using ResultPtr = TestU01Ptr< Result, Deleter >;

    //! \brief Default constructor taking no arguments
    //! \details Required as migratable. Called by Charm++. Initialize what we
    //!    can.
    explicit TestU01Props() :
      m_proxy(),
      m_rng( tk::ctr::RNGType::NO_RNG ),
      m_names(),
      m_xargs(),
      m_gen( nullptr ),
      m_runner( g_testStack.TestU01.runner.get<Test>() ),
      m_res( ResultPtr(Creator()) ),
      m_time( 0.0 ) {}

    //! \brief Initializer constructor
    //! \details None of the state data is const since the this class is
    //!   designed to be migratable over the network by the Charm++ runtime
    //!   system.
    //! \param[in] host Host proxy facilitating call-back to host object chare.
    //! \param[in] rng Random number generator ID enum to be tested
    //! \param[in] n Vector of statisical test names (can be more than one
    //!   associated with a given test, since a test can contain more than one
    //!   statistical test evaluation, yielding multiple p-values)
    //! \param[in] gen Raw function pointer to TestU01 statistical test
    //! \param[in] xargs Extra arguments to test-run
    explicit TestU01Props( Proxy& host,
                           tk::ctr::RNGType rng,
                           std::vector< std::string >&& n,
                           unif01_Gen* gen,
                           Ts&&... xargs ) :
      m_proxy( host ),
      m_rng( rng ),
      m_names( std::move(n) ),
      m_xargs( std::forward<Ts>(xargs)... ),
      m_gen( gen ),
      m_runner( g_testStack.TestU01.runner.get<Test>() ),
      m_res( ResultPtr(Creator()) ),
      m_time( 0.0 ) {}

    //! Copy assignment
    TestU01Props& operator=( const TestU01Props& x) {
      m_proxy = x.m_proxy;
      m_rng = x.m_rng;
      m_names = x.m_names;
      m_xargs = x.m_xargs;
      m_gen = x.m_gen;
      m_runner = x.m_runner;
      m_res = ResultPtr(Creator());
      m_time = x.m_time;
      return *this;
    }

    //! Copy constructor: in terms of copy assignment
    TestU01Props( const TestU01Props& x ) { operator=(x); }

    //! Move assignment
    TestU01Props& operator=( TestU01Props&& ) = default;
    //! Move constructor
    TestU01Props( TestU01Props&& ) = default;

    /** @name Pack/Unpack: Serialize Term object for Charm++ */
    ///@{
    //! Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er& p ) {
      p | m_proxy;
      p | m_rng;
      p | m_names;
      p | m_xargs;
      if (p.isUnpacking()) {
        m_gen = g_testStack.TestU01.generator( m_rng );
        m_runner = g_testStack.TestU01.runner.get< Test >();
        m_res = ResultPtr( Creator() );
      }
      p | m_time;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] c TestU01Props object reference
    friend void operator|( PUP::er& p, TestU01Props& c ) { c.pup(p); }
    ///@}

    //! Host proxy accessor
    //! \return Host proxy
    Proxy& proxy() noexcept { return m_proxy; }

    //! Number of results/test (i.e., p-values) accessor
    //! \return Number of p-values this test yields
    std::size_t npval() const { return m_names.size(); }

    //! Test name(s) accessor
    //! \return Vector of test names (there can be more than one)
    std::vector< std::string > names() const { return m_names; }

    //! \brief Run test and return its status as a vector of vector of strings.
    //! \details The return type could potentially be a more specific type,
    //!   e.g., a struct with fields RNGType, and a vector of strings for the
    //!   p-values, and the names. Instead, we return a more generic vector of
    //!   vector of strings type. This helps keeping the corresponding argument
    //!   to Battery::evaluate() generic, which may come in handy when other
    //!   test libraries are added in the future. The return value is thus the
    //!   following in a vector of vectors:
    //! - 0: Name(s) of statistics (note that a test may produce more than one
    //!   statistics and thus p-values)
    //! - 1: p-value strings: "pass" or "fail, p-value = ..." for all tests run
    //!   (note that a test may produce more than one statistics and thus
    //!   p-values)
    //! - 2: RNG name used to run the test: sub-vector length = 1
    std::vector< std::vector< std::string > > run() {
      // Run and time statistical test
      tk::Timer timer;
      const auto pvals = m_runner( m_gen, m_res.get(), m_xargs );
      m_time = timer.dsec();
      // Construct status
      std::vector< std::string > pvalstrs;
      for (std::size_t p=0; p<m_names.size(); ++p) {
        if ( fail(pvals[p]) )
          pvalstrs.emplace_back( "fail, p-value = " + pval( pvals[p] ) );
        else
          pvalstrs.emplace_back( "pass" );
      }
      return { m_names, pvalstrs, { tk::ctr::RNG().name(m_rng) } };
    }

    //! Test time accessor
    //! \return Time it took to run the test with the associated RNG name
    std::pair< std::string, tk::real > time()
    { return { tk::ctr::RNG().name(m_rng), m_time }; }

  private:
    //! \brief Pack/Unpack TestU01 external generator pointer
    //! \details Admittedly, the code below is ugly and looks stupid at first
    //!   sight. However, this is a translation of runtime information (a
    //!   user-selected RNG that this statistical test exercises) to
    //!   compile-time information: associating an RNG id from an enum class,
    //!   tk::ctr::RNGType::value, to a compile-time constant, underlying_type
    //!   value, facilitating a different global-scope TestU01 external
    //!   generator with code reuse. Note that createTestU01Gen() must be
    //!   global-scope as it is used to create a different external generator to
    //!   TestU01, depending on the RNG. Templating createTestU01Gen on the id
    //!   enables the compiler to generate a different wrapper for a different
    //!   RNG facilitating simultaneous calls to any or all wrappers as they are
    //!   unique functions.
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] g Reference to raw function pointer to TestU01
    //!   statistical test
    void pup( PUP::er& p, unif01_Gen*& g ) {
      IGNORE(p);
      using tk::ctr::RNGType;
      using tk::ctr::raw;
      const auto& rngname = tk::ctr::RNG().name(m_rng);
      Assert( !rngname.empty(), "Empty RNG name not allowed");
      Assert( m_rng != RNGType::NO_RNG, "No RNG type not allowed");
      #ifdef HAS_MKL
      if (m_rng == RNGType::MKL_MCG31)
        g = createTestU01Gen< raw(RNGType::MKL_MCG31) >( rngname );
      else if (m_rng == RNGType::MKL_R250)
        g = createTestU01Gen< raw(RNGType::MKL_R250) >( rngname );
      else if (m_rng == RNGType::MKL_MRG32K3A)
        g = createTestU01Gen< raw(RNGType::MKL_MRG32K3A) >( rngname );
      else if (m_rng == RNGType::MKL_MCG59)
        g = createTestU01Gen< raw(RNGType::MKL_MCG59) >( rngname );
      else if (m_rng == RNGType::MKL_WH) 
        g = createTestU01Gen< raw(RNGType::MKL_WH) >( rngname );
      else if (m_rng == RNGType::MKL_MT19937)
        g = createTestU01Gen< raw(RNGType::MKL_MT19937) >( rngname );
      else if (m_rng == RNGType::MKL_MT2203)
        g = createTestU01Gen< raw(RNGType::MKL_MT2203) >( rngname );
      else if (m_rng == RNGType::MKL_SFMT19937)
        g = createTestU01Gen< raw(RNGType::MKL_SFMT19937) >( rngname );
      else if (m_rng == RNGType::MKL_SOBOL)
        g = createTestU01Gen< raw(RNGType::MKL_SOBOL) >( rngname );
      else if (m_rng == RNGType::MKL_NIEDERR)
        g = createTestU01Gen< raw(RNGType::MKL_NIEDERR) >( rngname );
      //else if (m_rng == RNGType::MKL_IABSTRACT)
      //  g = createTestU01Gen< raw(RNGType::MKL_IABSTRACT) >( rngname );
      //else if (m_rng == RNGType::MKL_DABSTRACT)
      //  g = createTestU01Gen< raw(RNGType::MKL_DABSTRACT) >( rngname );
      //else if (m_rng == RNGType::MKL_SABSTRACT)
      //  g = createTestU01Gen< raw(RNGType::MKL_SABSTRACT) >( rngname );
      else if (m_rng == RNGType::MKL_NONDETERM)
        g = createTestU01Gen< raw(RNGType::MKL_NONDETERM) >( rngname );
      else
      #endif
      #ifdef HAS_RNGSSE2
      if (m_rng == RNGType::RNGSSE_GM19)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM19) >( rngname );
      else if (m_rng == RNGType::RNGSSE_GM29)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM29) >( rngname );
      else if (m_rng == RNGType::RNGSSE_GM31)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM31) >( rngname );
      else if (m_rng == RNGType::RNGSSE_GM55)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM55) >( rngname );
      else if (m_rng == RNGType::RNGSSE_GM61)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM61) >( rngname );
      else if (m_rng == RNGType::RNGSSE_GQ581)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GQ581) >( rngname );
      else if (m_rng == RNGType::RNGSSE_GQ583)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GQ583) >( rngname );
      else if (m_rng == RNGType::RNGSSE_GQ584)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GQ584) >( rngname );
      else if (m_rng == RNGType::RNGSSE_MT19937)
        g = createTestU01Gen< raw(RNGType::RNGSSE_MT19937) >( rngname );
      else if (m_rng == RNGType::RNGSSE_LFSR113)
        g = createTestU01Gen< raw(RNGType::RNGSSE_LFSR113) >( rngname );
      else if (m_rng == RNGType::RNGSSE_MRG32K3A)
        g = createTestU01Gen< raw(RNGType::RNGSSE_MRG32K3A) >( rngname );
      else
      #endif
      if (m_rng == RNGType::R123_THREEFRY)
        g = createTestU01Gen< raw(RNGType::R123_THREEFRY) >( rngname );
      else if (m_rng == RNGType::R123_PHILOX)
        g = createTestU01Gen< raw(RNGType::R123_PHILOX) >( rngname );
    }

    //! Query whether test is failed
    //! \param[in] val p-value returned from test
    //! \return true if the RNG has failed the statistical test
    bool fail( double val ) const {
      if ((val <= gofw_Suspectp) || (val >= 1.0-gofw_Suspectp))
        return true;
      else
        return false;
    }

    //! Return human-readable p-value (ala TestU01::bbattery.c::WritePval)
    //! \param[in] val p-value returned from test
    //! \return Human-readable p-value (ala TestU01::bbattery.c::WritePval)
    std::string pval( double val ) const {
      std::stringstream ss;
      if (val < gofw_Suspectp) {

        if ((val >= 0.01) && (val <= 0.99))
          ss << val;
        else if (val < gofw_Epsilonp)
          ss << "eps";
        else if (val < 0.01)
          ss << val;
        else if (val >= 1.0 - gofw_Epsilonp1)
          ss << "1 - eps1";
        else if (val < 1.0 - 1.0e-4)
          ss << val;
        else
          ss << 1.0 - val;

      } else if (val > 1.0 - gofw_Suspectp) {

        if (val >= 1.0 - gofw_Epsilonp1)
          ss << "1 - eps1";
        else if (val >= 1.0 - 1.0e-4)
          ss << "1 - " << 1.0 - val;
        else
          ss << val;

      }
      return ss.str();
    }

    Proxy m_proxy;                      //!< Host proxy
    tk::ctr::RNGType m_rng;             //!< RNG id
    std::vector< std::string > m_names; //!< Name(s) of tests
    Xargs m_xargs;                      //!< Extra args for run()
    unif01_Gen* m_gen;                  //!< Pointer to generator
    RunFn m_runner;                     //!< Test runner function
    ResultPtr m_res;                    //!< TestU01 results
    tk::real m_time;                    //!< Test run time measured in seconds
};

} // rngtest::

#endif // TestU01Props_h
