//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Props.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 10:31:15 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 statistical test properties
  \details   TestU01 statistical test properties
*/
//******************************************************************************
#ifndef TestU01Props_h
#define TestU01Props_h

#include <vector>
#include <functional>

#include <RNG.h>
#include <Options/RNG.h>
#include <TestU01Util.h>
#include <TestStack.h>

namespace rngtest {

extern TestStack g_testStack;

template< tk::ctr::RawRNGType id >
extern unif01_Gen* createTestU01Gen( const std::string& name );

//! TestU01 properties used to initialize TestU01 tests. Used to abstract away
//! the initialization for TestU01 statistical tests. Needed because of the move
//! semantics and variadic templates since Charm++ chares do not support these
//! yet. TestU01Props is therefore not a chare, but TestU01, initialized with a
//! TestU01Props object, is. Note that TestU01Props still needs to be
//! migratable, i.e., defines the pack/unpack operator, as it is an argument to
//! chare TestU01's constructor.
template< class Test,                   //!< Statistical test tag, struct{}
          class Proxy,                  //!< Host proxy
          class Result,                 //!< Results type
          Result* (*Creator)(void),     //!< Results creator
          void (*Deleter)(Result *),    //!< Results deleter
          typename... Ts >              //!< Extra runner args types
class TestU01Props {

  public:
    //! Test extra arguments type
    using Xargs = std::tuple< Ts... >;
    //! Test runner function pointer type
    using RunFn = std::vector<double> (*)( unif01_Gen*, Result*, const Xargs& );
    //! TestU01 results type with a custom deleter by TestU01
    using ResultPtr = TestU01Ptr< Result, Deleter >;

    //! Constructor taking no args, required as migratable, init what we can
    explicit TestU01Props() :
      m_rng( tk::ctr::RNGType::NO_RNG ),
      m_rngname( tk::ctr::RNG().name(m_rng) ),
      m_runner( g_testStack.TestU01.runner.get<Test>() ),
      m_res( ResultPtr(Creator()) ) {}

    //! Constructor
    explicit TestU01Props( Proxy& proxy,
                           tk::ctr::RNGType rng,
                           std::vector< std::string >&& names,
                           unif01_Gen* gen,
                           Ts&&... xargs ) :
      m_proxy( proxy ),
      m_rng( rng ),
      m_rngname( tk::ctr::RNG().name(rng) ),
      m_names( std::move(names) ),
      m_xargs( std::forward<Ts>(xargs)... ),
      m_npval( m_names.size() ),
      m_pvals( m_names.size(), -1.0 ),
      m_gen( gen ),
      m_runner( g_testStack.TestU01.runner.get<Test>() ),
      m_res( ResultPtr(Creator()) ) {}

    //! Copy assignment
    TestU01Props& operator=( const TestU01Props& x) {
      m_proxy = x.m_proxy;
      m_rng = x.m_rng;
      m_rngname = x.m_rngname;
      m_names = x.m_names;
      m_xargs = x.m_xargs;
      m_npval = x.m_npval;
      m_pvals = x.m_pvals;
      m_gen = x.m_gen;
      m_runner = x.m_runner;
      m_res = ResultPtr(Creator());
      return *this;
    }
    //! Copy constructor: in terms of copy assignment
    TestU01Props( const TestU01Props& x ) { operator=(x); }
    //! Move assignment
    TestU01Props& operator=( TestU01Props&& ) = default;
    //! Move constructor
    TestU01Props( TestU01Props&& ) = default;

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      p | m_proxy;
      using tk::operator|;
      p | m_rng;
      p | m_rngname;
      p | m_names;
      p | m_xargs;
      p | m_npval;
      p | m_pvals;
      if (p.isUnpacking()) {
        pup( p, m_gen );
        m_runner = g_testStack.TestU01.runner.get<Test>();
        m_res = ResultPtr(Creator());
      }
    }
    friend void operator|( PUP::er& p, TestU01Props& c ) { c.pup(p); }

    //! Host proxy accessor
    Proxy& proxy() noexcept { return m_proxy; }

    //! Run test
    void run() { m_pvals = m_runner( m_gen, m_res.get(), m_xargs ); }

    //! Number of results/test (i.e., p-values) accessor
    std::size_t npval() const { return m_npval; }

    //! Test name(s) accessor
    std::string name() const {
      // serialize vector of strings into a single semi-colon-separated string
      std::string names;
      for (const auto& n : m_names) names += n + ";";
      return names;
    }

  private:
    //! Pack/Unpack TestU01 external generator pointer
    void pup( PUP::er& p, unif01_Gen*& g ) {
      // Admittedly, the code below is ugly and looks stupid at first sight.
      // However, this is a translation of runtime information (a user-selected
      // RNG that this statistical test exercises) to compile-time information:
      // associating an RNG id from an enum class, tk::ctr::RNGType::value, to a
      // compile-time constant, underlying_type value, facilitating a different
      // global-scope TestU01 external generator with code reuse. Note that
      // createTestU01Gen must be global-scope as it is used to create a
      // different external generator to TestU01, depending on the RNG.
      // Templating createTestU01Gen on the id enables the compiler generate a
      // different wrapper for a different RNG facilitating simultaneous calls
      // to any or all wrappers as they are unique functions.
      using tk::ExceptType;
      using tk::ctr::RNGType;
      using tk::ctr::raw;
      Assert(!m_rngname.empty(), ExceptType::FATAL, "Empty RNG name not allowed");
      Assert(m_rng != RNGType::NO_RNG, ExceptType::FATAL, "No RNG type allowed");
      #ifdef HAS_MKL
      if (m_rng == RNGType::MKL_MCG31)
        g = createTestU01Gen< raw(RNGType::MKL_MCG31) >( m_rngname );
      else if (m_rng == RNGType::MKL_R250)
        g = createTestU01Gen< raw(RNGType::MKL_R250) >( m_rngname );
      else if (m_rng == RNGType::MKL_MRG32K3A)
        g = createTestU01Gen< raw(RNGType::MKL_MRG32K3A) >( m_rngname );
      else if (m_rng == RNGType::MKL_MCG59)
        g = createTestU01Gen< raw(RNGType::MKL_MCG59) >( m_rngname );
      else if (m_rng == RNGType::MKL_WH) 
        g = createTestU01Gen< raw(RNGType::MKL_WH) >( m_rngname );
      else if (m_rng == RNGType::MKL_MT19937)
        g = createTestU01Gen< raw(RNGType::MKL_MT19937) >( m_rngname );
      else if (m_rng == RNGType::MKL_MT2203)
        g = createTestU01Gen< raw(RNGType::MKL_MT2203) >( m_rngname );
      else if (m_rng == RNGType::MKL_SFMT19937)
        g = createTestU01Gen< raw(RNGType::MKL_SFMT19937) >( m_rngname );
      else if (m_rng == RNGType::MKL_SOBOL)
        g = createTestU01Gen< raw(RNGType::MKL_SOBOL) >( m_rngname );
      else if (m_rng == RNGType::MKL_NIEDERR)
        g = createTestU01Gen< raw(RNGType::MKL_NIEDERR) >( m_rngname );
      //else if (m_rng == RNGType::MKL_IABSTRACT)
      //  g = createTestU01Gen< raw(RNGType::MKL_IABSTRACT) >( m_rngname );
      //else if (m_rng == RNGType::MKL_DABSTRACT)
      //  g = createTestU01Gen< raw(RNGType::MKL_DABSTRACT) >( m_rngname );
      //else if (m_rng == RNGType::MKL_SABSTRACT)
      //  g = createTestU01Gen< raw(RNGType::MKL_SABSTRACT) >( m_rngname );
      else if (m_rng == RNGType::MKL_NONDETERM)
        g = createTestU01Gen< raw(RNGType::MKL_NONDETERM) >( m_rngname );
      #endif
      else if (m_rng == RNGType::RNGSSE_GM19)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM19) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_GM29)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM29) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_GM31)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM31) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_GM55)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM55) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_GM61)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GM61) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_GQ581)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GQ581) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_GQ583)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GQ583) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_GQ584)
        g = createTestU01Gen< raw(RNGType::RNGSSE_GQ584) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_MT19937)
        g = createTestU01Gen< raw(RNGType::RNGSSE_MT19937) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_LFSR113)
        g = createTestU01Gen< raw(RNGType::RNGSSE_LFSR113) >( m_rngname );
      else if (m_rng == RNGType::RNGSSE_MRG32K3A)
        g = createTestU01Gen< raw(RNGType::RNGSSE_MRG32K3A) >( m_rngname );
    }

    Proxy m_proxy;                      //!< Host proxy
    tk::ctr::RNGType m_rng;             //!< RNG id
    std::string m_rngname;              //!< name of RNG
    std::vector< std::string > m_names; //!< name(s) of tests
    Xargs m_xargs;                      //!< extra args for run()
    std::size_t m_npval;                //!< number of p-values produced
    std::vector< double > m_pvals;      //!< p-value(s)
    unif01_Gen* m_gen;                  //!< pointer to generator
    RunFn m_runner;                     //!< test runner function
    ResultPtr m_res;                    //!< testU01 results
};

} // rngtest::

#endif // TestU01Props_h
