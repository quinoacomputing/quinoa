//******************************************************************************
/*!
  \file      src/RNG/TestU01Suite.h
  \author    J. Bakosi
  \date      Mon 02 Dec 2013 09:58:28 PM MST
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
}

#include <Battery.h>
#include <TestU01Util.h>

namespace rngtest {

//! Global pointers to RNGs (in TestU01Suite.C)
extern std::vector< std::unique_ptr< tk::RNG > > g_rng;

//! TestU01 random number generator test suite
class TestU01Suite : public Battery {

  protected:
    //! Constructor
    explicit TestU01Suite(const Base& base);

    //! Destructor
    ~TestU01Suite() noexcept override = default;

    //! TestU01 external generator type with a custom deleter by TestU01
    using Gen01Ptr = TestU01Ptr< unif01_Gen, unif01_DeleteExternGen01 >;
    //! TestU01 external RNGs (associated to all of quinoa's RNGs)
    std::vector< Gen01Ptr > m_gen;
    //! TestU01 external RNGs associated to an integer id
    std::map< quinoa::ctr::RNGType, int > m_rng;

    //! Ids of RNGs tested
    std::vector< int > m_testRNGs;

    template< int id>
    void addRNG( quinoa::ctr::RNGType rng,
                 double (*wrap)(void*,void*),
                 unsigned long (*wrap_bits)(void*,void*)) {
      // Create new RNG and store its pointer in global scope
      g_rng.push_back( std::unique_ptr< tk::RNG >( m_base.rng[rng]() ) );
      // Create new TestU01 external RNG and associate global-scope wrapper
      char* const name = const_cast<char*>( quinoa::ctr::RNG().name(rng).c_str() );
      m_gen.push_back(
        Gen01Ptr( unif01_CreateExternGen01( name, wrap, wrap_bits ) ) );
      // Associate rng to id
      m_rng[ rng ] = id;
    }

    //! Statistical tests wrappers
    using Pvals = StatTest::Pvals;
    static Pvals BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res );
    static Pvals Collision( unif01_Gen* gen, sknuth_Res2* res );
    static Pvals Gap( unif01_Gen* gen, sres_Chi2* res );
    static Pvals SimpPoker( unif01_Gen* gen, sres_Chi2* res );
    static Pvals CouponCollector( unif01_Gen* gen, sres_Chi2* res );
    static Pvals MaxOft( unif01_Gen* gen, sknuth_Res1* res );
    static Pvals WeightDistrib( unif01_Gen* gen, sres_Chi2* res );
    static Pvals MatrixRank( unif01_Gen* gen, sres_Chi2* res );
    static Pvals HammingIndep( unif01_Gen* gen, sstring_Res* res );
    static Pvals RandomWalk1( unif01_Gen* gen, swalk_Res* res );

  private:
    //! Don't permit copy constructor
    TestU01Suite(const TestU01Suite&) = delete;
    //! Don't permit copy assigment
    TestU01Suite& operator=(const TestU01Suite&) = delete;
    //! Don't permit move constructor
    TestU01Suite(TestU01Suite&&) = delete;
    //! Don't permit move assigment
    TestU01Suite& operator=(TestU01Suite&&) = delete;

    static const long THOUSAND = 1000;
    static const long MILLION = THOUSAND * THOUSAND;
    static const long BILLION = THOUSAND * MILLION;
};

} // rngtest::

#endif // TestU01Suite_h
