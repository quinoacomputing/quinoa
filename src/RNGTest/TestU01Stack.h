// *****************************************************************************
/*!
  \file      src/RNGTest/TestU01Stack.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Stack of TestU01 RNG statistical tests
  \details   Stack of TestU01 RNG statistical tests
*/
// *****************************************************************************
#ifndef TestU01Stack_h
#define TestU01Stack_h

#include <vector>
#include <functional>
#include <iosfwd>
#include <map>
#include <tuple>
#include <type_traits>

#include "NoWarning/value_factory.h"

extern "C" {
  #include <sres.h>
  #include <sstring.h>
  #include <sknuth.h>
  #include <swalk.h>
  #include <smarsa.h>
  #include <scomp.h>
  #include <sspectral.h>
  #include <unif01.h>
  #include <snpair.h>
}

#include "StatTest.h"
#include "TestU01Util.h"
#include "Tags.h"
#include "TaggedTuple.h"
#include "Options/RNG.h"

namespace rngtest {

//! Stack of TestU01 RNG statistical tests
class TestU01Stack {

  public:
    //! Constructor
    explicit TestU01Stack();

    //! \brief Add a statistical test to battery
    //! \details Note that adding a test to the battery does not invoke the test
    //!   constructor, it only records the information on how to call the test
    //!   constructor in the future. That is it binds the constructor arguments
    //!   to the constructor call and records the the information so only a
    //!   function call "()" is necessary to instantiate it.
    //! \param[in] proxy Charm++ host proxy to which the test calls back to
    //! \param[in] tests Vector of test constructors to add tests to
    //! \param[in] r RNG ID enum
    //! \param[in] gen Raw function pointer to TestU01 statistical test
    //! \param[in] names Vector of statisical test names (can be more than one
    //!   associated with a given test, since a test can contain more than one
    //!   statistical test evaluation, yielding multiple p-values)
    //! \param[in] xargs Extra arguments to test-run
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

    /** @name Stack of TestU01 statistical tests wrappers
      * */
    ///@{
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
    SimplePoker( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, int, int>& xargs );

    static std::vector< double >
    CouponCollector( unif01_Gen* gen, sres_Chi2* res,
      const std::tuple<long, long, int, int>& xargs );

    static std::vector< double >
    MaxOft( unif01_Gen* gen, sknuth_Res1* res,
      const std::tuple<long, long, int, int, int, int, int>& xargs );

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
    AutoCorr( unif01_Gen* gen, sres_Basic* res,
      const std::tuple<long, long, int, int, int>& xargs );
    ///@}

    //! Function pointer type to use to define pointers to test runner functions
    //! \note Abstracting the function pointer away as this, also pleases
    //!    doxygen in client code, otherwise doxygen issues warnings, such as
    //!    `no matching class member found for` as it fails to parse the
    //!    function pointer correctly and fails to find a function declaration,
    //!    which, of course, does not exist. Thanks to Vladimír Vondruš,
    //!    author of m.css (http://mcss.mosra.cz) for the help with this.
    template< class... Args > using FnPtr = std::vector<double>(*)( Args... );

    //! \brief Compile-time tag-based access to individual test wrappers.
    //! \details This tagged_tuple is practically a compile-time map that
    //!   associates tags (empty structs) to test wrappers. This is used to find
    //!   the test wrapper function pointers after migration over the network.
    //! \see See also TestU01Props::pup().
    tk::tuple::tagged_tuple<

      tag::BirthdaySpacings,                                // tag
      FnPtr< unif01_Gen*, sres_Poisson*,
             const std::tuple<long, long, int, long, int, int>& >,

      tag::Collision,
      FnPtr< unif01_Gen*, sknuth_Res2*,
             const std::tuple<long, long, int, long, int>& >,

      tag::RandomWalk1,
      FnPtr< unif01_Gen*, swalk_Res*,
             const std::tuple<long, long, int, int, long, long>& >,

      tag::Gap,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, double, double>& >,

      tag::SimplePoker,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, int, int>& >,

      tag::CouponCollector,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, int>& >,

      tag::MaxOft,
      FnPtr< unif01_Gen*, sknuth_Res1*,
             const std::tuple<long, long, int, int, int, int, int>& >,

      tag::WeightDistrib,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, long, double, double>& >,

      tag::MatrixRank,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, int, int, int>& >,

      tag::HammingIndep,
      FnPtr< unif01_Gen*, sstring_Res*,
             const std::tuple<long, long, int, int, int, int>& >,

      tag::SerialOver,
      FnPtr< unif01_Gen*, sres_Basic*,
             const std::tuple<long, long, int, long, int>& >,

      tag::CollisionOver,
      FnPtr< unif01_Gen*, smarsa_Res*,
             const std::tuple<long, long, int, long, int>& >,

      tag::ClosePairs,
      FnPtr< unif01_Gen*, snpair_Res*,
             const std::tuple<long, long, int, int, int, int, int>& >,

      tag::ClosePairsBitMatch,
      FnPtr< unif01_Gen*, snpair_Res*,
             const std::tuple<long, long, int, int>& >,

      tag::Run,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, int>& >,

      tag::Permutation,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, int>& >,

      tag::CollisionPermut,
      FnPtr< unif01_Gen*, sknuth_Res2*,
             const std::tuple<long, long, int, int>& >,

      tag::SampleProd,
      FnPtr< unif01_Gen*, sres_Basic*,
             const std::tuple<long, long, int, int>& >,

      tag::SampleMean,
      FnPtr< unif01_Gen*, sres_Basic*,
             const std::tuple<long, long, int>& >,

      tag::SampleCorr,
      FnPtr< unif01_Gen*, sres_Basic*,
             const std::tuple<long, long, int, int>& >,

      tag::AppearanceSpacings,
      FnPtr< unif01_Gen*, sres_Basic*,
             const std::tuple<long, long, long, int, int, int>& >,

      tag::SumCollector,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, double>& >,

      tag::Savir2,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, long, int>& >,

      tag::GCD,
      FnPtr< unif01_Gen*, smarsa_Res2*,
             const std::tuple<long, long, int, int>& >,

      tag::LinearComp,
      FnPtr< unif01_Gen*, scomp_Res*,
             const std::tuple<long, long, int, int>& >,

      tag::LempelZiv,
      FnPtr< unif01_Gen*, sres_Basic*,
             const std::tuple<long, int, int, int>& >,

      tag::Fourier3,
      FnPtr< unif01_Gen*, sspectral_Res*,
             const std::tuple<long, int, int, int>& >,

      tag::LongestHeadRun,
      FnPtr< unif01_Gen*, sstring_Res2*,
             const std::tuple<long, long, int, int, long>& >,

      tag::PeriodsInStrings,
      FnPtr< unif01_Gen*, sres_Chi2*,
             const std::tuple<long, long, int, int>& >,

      tag::HammingWeight2,
      FnPtr< unif01_Gen*, sres_Basic*,
             const std::tuple<long, long, int, int, long>& >,

      tag::HammingCorr,
      FnPtr< unif01_Gen*, sstring_Res*,
             const std::tuple<long, long, int, int, int>& >,

      tag::StringRun,
      FnPtr< unif01_Gen*, sstring_Res3*,
             const std::tuple<long, long, int, int>& >,

      tag::AutoCorr,
      FnPtr< unif01_Gen*, sres_Basic*,
             const std::tuple<long, long, int, int, int>& >

    > runner {

      BirthdaySpacings,   // Initialize by binding to member function wrappers.
      Collision,          // Obviously the order here is important.
      RandomWalk1,
      Gap,
      SimplePoker,
      CouponCollector,
      MaxOft,
      WeightDistrib,
      MatrixRank,
      HammingIndep,
      SerialOver,
      CollisionOver,
      ClosePairs,
      ClosePairsBitMatch,
      Run,
      Permutation,
      CollisionPermut,
      SampleProd,
      SampleMean,
      SampleCorr,
      AppearanceSpacings,
      SumCollector,
      Savir2,
      GCD,
      LinearComp,
      LempelZiv,
      Fourier3,
      LongestHeadRun,
      PeriodsInStrings,
      HammingWeight2,
      HammingCorr,
      StringRun,
      AutoCorr

    };

    //! Find RNG properties based on RNG id
    unif01_Gen* generator( tk::ctr::RNGType r ) const;

  private:
    //! Create TestU01 RNG wrapper
    template< tk::ctr::RawRNGType id > void addRNG( tk::ctr::RNGType r );

    std::map< tk::ctr::RNGType, Gen01Ptr > m_generator; //!< RNG wrappers
};

} // rngtest::

#endif // TestU01Stack_h
