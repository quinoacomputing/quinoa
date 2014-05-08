//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.h
  \author    J. Bakosi
  \date      Wed 07 May 2014 10:43:01 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 random number generator test suite
  \details   TestU01 random number generator test suite
*/
//******************************************************************************
#ifndef TestU01Suite_h
#define TestU01Suite_h

extern "C" {
  #include <unif01.h>
  #include <sres.h>
  #include <sstring.h>
  #include <sknuth.h>
  #include <swalk.h>
  #include <smarsa.h>
  #include <snpair.h>
  #include <scomp.h>
  #include <sspectral.h>
}

#include <Battery.h>
#include <TestU01Util.h>
#include <testu01suite.decl.h>

namespace rngtest {

//! TestChareContainer
class TestChareContainer : public CBase_TestChareContainer {

  public:
    //! Constructor
    explicit TestChareContainer() {}

    //! Migrator
    TestChareContainer( CkMigrateMessage* ) {}

  private:
    //! Don't permit copy constructor
    TestChareContainer(const TestChareContainer&) = delete;
    //! Don't permit copy assigment
    TestChareContainer& operator=(const TestChareContainer&) = delete;
    //! Don't permit move constructor
    TestChareContainer(TestChareContainer&&) = delete;
    //! Don't permit move assigment
    TestChareContainer& operator=(TestChareContainer&&) = delete;
};


//! Global pointers to RNGs (in TestU01Suite.C)
extern std::vector< std::unique_ptr< tk::RNG > > g_rng;


//! TestU01 random number generator test suite
class TestU01Suite : public Battery {

  public:
    //! Run battery of RNG tests
    void run() override;

    //! Print list of registered statistical tests
    void print() const override;

    //! Return number of statistical tests in battery
    std::size_t ntest() const override {
      return m_numRNGs ? m_tests.size() / m_numRNGs : 0;
    }

    //! Return number of statistics produced by battery
    std::size_t nstat() const override {
      return m_numRNGs ? m_npval / m_numRNGs : 0;
    }

  protected:
    //! Constructor
    explicit TestU01Suite(const Base& base, const std::string& name);

    //! Destructor
    // ICC should be default
    ~TestU01Suite() override {}

    //! TestU01 external generator type with a custom deleter by TestU01
    using Gen01Ptr = TestU01Ptr< unif01_Gen, unif01_DeleteExternGen01 >;

    using Pvals = StatTest::Pvals;

    //! Statistical tests wrappers
    static Pvals BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res,
               const std::tuple<long, long, int, long, int, int>& xargs );
    static Pvals Collision( unif01_Gen* gen, sknuth_Res2* res,
               const std::tuple<long, long, int, long, int>& xargs );
    static Pvals Gap( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, double, double>& xargs );
    static Pvals SimpPoker( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, int, int>& xargs );
    static Pvals CouponCollector( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals MaxOft( unif01_Gen* gen, sknuth_Res1* res,
               const std::tuple<long, long, int, int, int, gofw_TestType,
                                gofw_TestType>& xargs );
    static Pvals WeightDistrib( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, long, double, double>& xargs );
    static Pvals MatrixRank( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, int, int, int>& xargs );
    static Pvals HammingIndep( unif01_Gen* gen, sstring_Res* res,
               const std::tuple<long, long, int, int, int, int>& xargs );
    static Pvals RandomWalk1( unif01_Gen* gen, swalk_Res* res,
               const std::tuple<long, long, int, int, long, long>& xargs );
    static Pvals SerialOver( unif01_Gen* gen, sres_Basic* res,
               const std::tuple<long,long, int, long, int>& xargs );
    static Pvals CollisionOver( unif01_Gen* gen, smarsa_Res* res,
               const std::tuple<long, long, int, long, int>& xargs );
    static Pvals ClosePairs( unif01_Gen* gen, snpair_Res* res,
               const std::tuple<long, long, int, int, int, int, int>& xargs );
    static Pvals ClosePairsBitMatch( unif01_Gen* gen, snpair_Res* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals Run( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals Permutation( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals CollisionPermut( unif01_Gen* gen, sknuth_Res2* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals SampleProd( unif01_Gen* gen, sres_Basic* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals SampleMean( unif01_Gen* gen, sres_Basic* res,
               const std::tuple<long, long, int>& xargs );
    static Pvals SampleCorr( unif01_Gen* gen, sres_Basic* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals AppearanceSpacings( unif01_Gen* gen, sres_Basic* res,
               const std::tuple<long, long, long, int, int, int>& xargs );
    static Pvals SumCollector( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, double>& xargs );
    static Pvals Savir2( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, long, int>& xargs );
    static Pvals GCD( unif01_Gen* gen, smarsa_Res2* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals LinearComp( unif01_Gen* gen, scomp_Res* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals LempelZiv( unif01_Gen* gen, sres_Basic* res,
               const std::tuple<long, int, int, int>& xargs );
    static Pvals Fourier3( unif01_Gen* gen, sspectral_Res* res,
               const std::tuple<long, int, int, int>& xargs );
    static Pvals LongestHeadRun( unif01_Gen* gen, sstring_Res2* res,
               const std::tuple<long, long, int, int, long>& xargs );
    static Pvals PeriodsInStrings( unif01_Gen* gen, sres_Chi2* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals HammingWeight2( unif01_Gen* gen, sres_Basic* res,
               const std::tuple<long, long, int, int, long>& xargs );
    static Pvals HammingCorr( unif01_Gen* gen, sstring_Res* res,
               const std::tuple<long, long, int, int, int>& xargs );
    static Pvals StringRun( unif01_Gen* gen, sstring_Res3* res,
               const std::tuple<long, long, int, int>& xargs );
    static Pvals AutoCor( unif01_Gen* gen, sres_Basic* res,
               const std::tuple<long, long, int, int, int>& xargs );

    //! Assign statistical tests to selected RNGs
    template< class Suite >
    void assignTests( Suite& suite )
    {
      auto size = m_rngprops.size();
       for (std::size_t r=0; r<size; ++r) {
         auto& rng = m_rngprops[r];
         if (rng.id != tk::ctr::RNGType::NO_RNG) {
           ++m_numRNGs;
           suite.addTests( r, rng.id, rng.ptr );
         }
       }
      total();  // Count up total number of statistics expected
    }

    //! Add statistical test to battery
    template< class TestType, class Result, typename... Ts >
    void add( std::size_t id,
              const Gen01Ptr& gen,
              const tk::ctr::RNGType& rng,
              StatTest::Names&& names,
              Pvals (*runner)(unif01_Gen*, Result*, const std::tuple<Ts...>&),
              Ts&&... xargs )
    {
      m_tests.push_back( 
        tk::make_unique< TestType >( id,
                                     gen.get(),
                                     std::move(rng),
                                     std::move(names),
                                     runner,
                                     std::forward<Ts>(xargs)... ) );
    }

    static const long THOUSAND = 1000;
    static const long MILLION = THOUSAND * THOUSAND;
    static const long BILLION = THOUSAND * MILLION;

  private:
    //! Don't permit copy constructor
    TestU01Suite(const TestU01Suite&) = delete;
    //! Don't permit copy assigment
    TestU01Suite& operator=(const TestU01Suite&) = delete;
    //! Don't permit move constructor
    TestU01Suite(TestU01Suite&&) = delete;
    //! Don't permit move assigment
    TestU01Suite& operator=(TestU01Suite&&) = delete;

    //! RNG properties
    struct RNGProps {
       using TimerIds = std::vector< tk::TimerId >;
       using Times = std::vector< tk::real >;

      // The first two fields, id and ptr, represent TestU01 external RNGs
      // associated to quinoa's RNGs. This could also be a single std::map,
      // however, using two equal-sized std::vectors instead. This is
      // due to the constraints on passing the RNGs as global wrappers to
      // TestU01. Only those entries of the std::vectors are initialized and
      // used that are requested by the user (know only at runtime). RNGs are
      // always assigned to the same position (hence std::vector instead of
      // std::map), regardless of requested or not. See also TestU01Suite
      // constructor.
      tk::ctr::RNGType id;
      Gen01Ptr ptr;

      // Timer ids and measured time as tk::real in seconds
      TimerIds timer;           // size: number of threads
      Times time;               // size: number of threads

      // Fill constructor, all arguments optional
      RNGProps( tk::ctr::RNGType i = tk::ctr::RNGType::NO_RNG,
                Gen01Ptr&& p = nullptr,
                TimerIds&& r = { 0 },
                Times&& t = { 0.0 } ) :
        id( i ),
        ptr( std::forward< Gen01Ptr >( p ) ),
        timer( std::forward< TimerIds >( r ) ),
        time( std::forward< Times >( t ) ) {}
    };

    //! Properties of all RNGs tested
    std::vector< RNGProps > m_rngprops;

    //! Number of RNGs tested
    std::size_t m_numRNGs;

    //! Register a RNG for testing
    template< int id >
    void addRNG( tk::ctr::RNGType r,
                 double (*wrap)(void*,void*),
                 unsigned long (*wrap_bits)(void*,void*) )
    {
      // Create new RNG and store its pointer in global scope
      g_rng[id] = std::unique_ptr< tk::RNG >( m_base.rng[r]() );
      // Get RNG's name
      std::string name( tk::ctr::RNG().name(r) );
      // Create per-thread timers named '<rng_name> <thread_id>' and zero times
      using tk::operator+;
      auto nthreads = m_base.paradigm.nthreads();
      RNGProps::TimerIds ids( nthreads );
      RNGProps::Times times( nthreads );
      for ( int i=0; i<nthreads; ++i ) {
        ids[i] = m_base.timer.create( name + std::string(" ") + i );
        times[i] = 0.0;
      }
      // Fill RNG properties for given id
      m_rngprops[id] = {
        r,
        Gen01Ptr( unif01_CreateExternGen01( const_cast<char*>(name.c_str()),
                                            wrap, wrap_bits ) ),
        std::move( ids ),
        std::move( times ) };
   }

   //! Count up total number of statistics expected
   void total();

   //! Count up number of failed tests
   std::size_t failed();

   const std::string m_name;             //!< Test suite name

   std::size_t m_npval;                  //!< Total number of stats
   TestContainer m_tests;                //!< Statistical tests
   TestChareContainer m_charetests;      //!< Statistical tests
};

} // rngtest::

#endif // TestU01Suite_h
