//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Stack.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 05:28:57 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Stack of TestU01 tests
  \details   Stack of TestU01 tests
*/
//******************************************************************************
#ifndef TestU01Stack_h
#define TestU01Stack_h

#include <vector>

#include <boost/functional/value_factory.hpp>

extern "C" {
  #include <sres.h>
  #include <sstring.h>
  #include <sknuth.h>
  #include <swalk.h>
  #include <smarsa.h>
  #include <snpair.h>
  #include <scomp.h>
  #include <sspectral.h>
}

#include <StatTest.h>
#include <TestU01Util.h>
#include <RNGTest/Tags.h>

namespace rngtest {

//! TestU01Stack
class TestU01Stack {

  public:
    //! Constructor
    explicit TestU01Stack();

    //! Add statistical test to battery
    template< class TestType, class Proxy, typename... Ts >
    void add( Proxy& proxy,
              std::vector< std::function< StatTest() > >& tests,
              tk::ctr::RNGType r,
              unif01_Gen* const gen,
              std::vector< std::string >&& names,
              Ts&&... xargs ) const
    {
      using Model = TestType;
      using Host = StatTest;
      using Props = typename TestType::Props;
      tests.emplace_back(
        std::bind( boost::value_factory< Host >(),
                   std::function< Model() >(),
                   std::forward< Props >(
                     Props( proxy, r, std::move(names), gen,
                            std::forward<Ts>(xargs)... ) ) ) );
    }

    //! Stack of TestU01 statistical tests wrappers

    static std::vector< double >
    BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res,
      const std::tuple<long, long, int, long, int, int>& xargs );

    static std::vector< double >
    Collision( unif01_Gen* gen, sknuth_Res2* res,
      const std::tuple<long, long, int, long, int>& xargs );

    static std::vector< double >
    Gap( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, double, double>& xargs );

    static std::vector< double >
    SimpPoker( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, int, int>& xargs );

    static std::vector< double >
    CouponCollector( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    MaxOft( unif01_Gen* gen, sknuth_Res1* res,
      const std::tuple<long, long, int, int, int, gofw_TestType,
                       gofw_TestType>& xargs );

    static std::vector< double >
    WeightDistrib( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, long, double, double>& xargs );

    static std::vector< double >
    MatrixRank( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, int, int, int>& xargs );

    static std::vector< double >
    HammingIndep( unif01_Gen* gen, sstring_Res* res,
      const std::tuple<long, long, int, int, int, int>& xargs );

    static std::vector< double >
    RandomWalk1( unif01_Gen* gen, swalk_Res* res,
      const std::tuple<long, long, int, int, long, long>& xargs );

    static std::vector< double >
    SerialOver( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, long, int, long, int>& xargs );

    static std::vector< double >
    CollisionOver( unif01_Gen* gen, smarsa_Res* res,
      const std::tuple<long, long, int, long, int>& xargs );

    static std::vector< double >
    ClosePairs( unif01_Gen* gen, snpair_Res* res,
      const std::tuple<long, long, int, int, int, int, int>& xargs );

    static std::vector< double >
    ClosePairsBitMatch( unif01_Gen* gen, snpair_Res* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    Run( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    Permutation( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    CollisionPermut( unif01_Gen* gen, sknuth_Res2* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    SampleProd( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    SampleMean( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, long, int>& xargs );

    static std::vector< double >
    SampleCorr( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    AppearanceSpacings( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, long, long, int, int, int>& xargs );

    static std::vector< double >
    SumCollector( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, double>& xargs );

    static std::vector< double >
    Savir2( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, long, int>& xargs );

    static std::vector< double >
    GCD( unif01_Gen* gen, smarsa_Res2* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    LinearComp( unif01_Gen* gen, scomp_Res* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    LempelZiv( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, int, int, int>& xargs );

    static std::vector< double >
    Fourier3( unif01_Gen* gen, sspectral_Res* res,
      const std::tuple<long, int, int, int>& xargs );

    static std::vector< double >
    LongestHeadRun( unif01_Gen* gen, sstring_Res2* res,
      const std::tuple<long, long, int, int, long>& xargs );

    static std::vector< double >
    PeriodsInStrings( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    HammingWeight2( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, long, int, int, long>& xargs );

    static std::vector< double >
    HammingCorr( unif01_Gen* gen, sstring_Res* res,
      const std::tuple<long, long, int, int, int>& xargs );

    static std::vector< double >
    StringRun( unif01_Gen* gen, sstring_Res3* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    AutoCor( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, long, int, int, int>& xargs );

    //! RNG properties
    struct RNGProps {
      Gen01Ptr ptr;     //!< TestU01 RNG wrapper
      RNGProps( unif01_Gen* p ) : ptr( Gen01Ptr(p) ) {}
    };

    //! Find RNG properties based on RNG id
    const RNGProps& rngprops( tk::ctr::RNGType r ) const;

    //! Compile-time tag-based access to individual wrappers
    //! Associate tags (empty structs) to test wrappers
    tk::tuple::tagged_tuple<

      tag::BirthdaySpacings,
      std::vector< double > (*)( unif01_Gen*, sres_Poisson*,
        const std::tuple<long, long, int, long, int, int>& ),

      tag::Collision,
      std::vector< double > (*)( unif01_Gen*, sknuth_Res2*,
        const std::tuple<long, long, int, long, int>& ),

      tag::SerialOver,
      std::vector< double > (*)( unif01_Gen*, sres_Basic*,
        const std::tuple<long, long, int, long, int>& ),

      tag::RandomWalk1,
      std::vector< double > (*)( unif01_Gen*, swalk_Res*,
        const std::tuple<long, long, int, int, long, long>& )

    > runner {

      this->BirthdaySpacings,
      this->Collision,
      this->SerialOver,
      this->RandomWalk1

    };

  private:
    //! Create TestU01 RNG wrapper
    template< tk::ctr::RawRNGType id > void addRNG( tk::ctr::RNGType r );

    std::map< tk::ctr::RNGType, RNGProps > m_rngprops; //!< Selected RNG props
};

} // rngtest::

#endif // TestU01Stack_h
