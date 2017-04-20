// *****************************************************************************
/*!
  \file      src/Control/Options/RNG.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator options and associations
  \details   Random number generator options and associations
*/
// *****************************************************************************
#ifndef RNGOptions_h
#define RNGOptions_h

#include <map>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/joint_view.hpp>

#include "QuinoaConfig.h"

#include "NoWarning/mkl_vsl.h"

#include "Toggle.h"
#include "Keywords.h"
#include "Options/RNGSSESeqLen.h"
#include "PUPUtil.h"
#include "StrConvUtil.h"

namespace tk {
namespace ctr {

//! Random number generator types
//! \author J. Bakosi
enum class RNGType : uint8_t { NO_RNG=0
                             #ifdef HAS_RNGSSE2
                             , RNGSSE_GM19
                             , RNGSSE_GM29
                             , RNGSSE_GM31
                             , RNGSSE_GM55
                             , RNGSSE_GM61
                             , RNGSSE_GQ581
                             , RNGSSE_GQ583
                             , RNGSSE_GQ584
                             , RNGSSE_MT19937
                             , RNGSSE_LFSR113
                             , RNGSSE_MRG32K3A
                             #endif
                             #ifdef HAS_MKL
                             , MKL_MCG31
                             , MKL_R250
                             , MKL_MRG32K3A
                             , MKL_MCG59
                             , MKL_WH
                             , MKL_MT19937
                             , MKL_MT2203
                             , MKL_SFMT19937
                             , MKL_SOBOL
                             , MKL_NIEDERR
                             // MKL VSL's abstract RNGs are for use with
                             // external (i.e., external to MKL or user-defined
                             // generator via call-back functions), disabled for
                             // now
                             //, MKL_IABSTRACT
                             //, MKL_DABSTRACT
                             //, MKL_SABSTRACT
                             , MKL_NONDETERM
                             #endif
                             , R123_THREEFRY
                             , R123_PHILOX
};

//! \brief Pack/Unpack RNGType: forward overload to generic enum class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, RNGType& e ) { PUP::pup( p, e ); }

//! Enum class underlying type shortcut
//! \author J. Bakosi
using RawRNGType = std::underlying_type< RNGType >::type;

//! Return underlying type
//! \param[in] r RNG enum class value
//! \return Enum class underlying value using static_cast
//! \author J. Bakosi
constexpr RawRNGType raw( RNGType r ) { return static_cast< RawRNGType >( r ); }

//! Random number generator library types
//! \author J. Bakosi
enum class RNGLibType : uint8_t { NO_LIB=0,
                                  MKL,
                                  RNGSSE,
                                  PRAND,
                                  R123 };

//! \brief RNG options: outsource searches to base templated on enum type
//! \author J. Bakosi
class RNG : public tk::Toggle< RNGType > {

  private:
    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywordsMKL = boost::mpl::vector<
                                          #ifdef HAS_MKL
                                            kw::mkl_mcg31
                                          , kw::mkl_r250
                                          , kw::mkl_mrg32k3a
                                          , kw::mkl_mcg59
                                          , kw::mkl_wh
                                          , kw::mkl_mt19937
                                          , kw::mkl_mt2203
                                          , kw::mkl_sfmt19937
                                          , kw::mkl_sobol
                                          , kw::mkl_niederr
                                          //, kw::mkl_iabstract
                                          //, kw::mkl_dabstract
                                          //, kw::mkl_sabstract
                                          , kw::mkl_nondeterm
                                          #endif
                                          >;
    using keywordsRNGSSE2 = boost::mpl::vector<
                                              #ifdef HAS_RNGSSE2
                                                kw::rngsse_gm19
                                              , kw::rngsse_gm29
                                              , kw::rngsse_gm31
                                              , kw::rngsse_gm55
                                              , kw::rngsse_gm61
                                              , kw::rngsse_gq581
                                              , kw::rngsse_gq583
                                              , kw::rngsse_gq584
                                              , kw::rngsse_mt19937
                                              , kw::rngsse_lfsr113
                                              , kw::rngsse_mrg32k3a
                                              #endif
                                              >;
    using keywordsR123 = boost::mpl::vector< kw::r123_threefry
                                           , kw::r123_philox >;

  public:
    using ParamType = int;
    using LibType = RNGLibType;

    using keywords = boost::mpl::joint_view< keywordsMKL,
                       boost::mpl::joint_view< keywordsRNGSSE2,
                                               keywordsR123 > >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit RNG() :
      tk::Toggle< RNGType >(
        //! Group, i.e., options, name
        "Random number generator",
        //! Enums -> names
        { { RNGType::NO_RNG, "n/a" }
        , { RNGType::R123_THREEFRY, kw::r123_threefry::name() }
        , { RNGType::R123_PHILOX, kw::r123_philox::name() }
        #ifdef HAS_RNGSSE2
        , { RNGType::RNGSSE_GM19, kw::rngsse_gm19::name() }
        , { RNGType::RNGSSE_GM29, kw::rngsse_gm29::name() }
        , { RNGType::RNGSSE_GM31, kw::rngsse_gm31::name() }
        , { RNGType::RNGSSE_GM55, kw::rngsse_gm55::name() }
        , { RNGType::RNGSSE_GM61, kw::rngsse_gm61::name() }
        , { RNGType::RNGSSE_GQ581, kw::rngsse_gq581::name() }
        , { RNGType::RNGSSE_GQ583, kw::rngsse_gq583::name() }
        , { RNGType::RNGSSE_GQ584, kw::rngsse_gq584::name() }
        , { RNGType::RNGSSE_MT19937, kw::rngsse_mt19937::name() }
        , { RNGType::RNGSSE_LFSR113, kw::rngsse_lfsr113::name() }
        , { RNGType::RNGSSE_MRG32K3A, kw::rngsse_mrg32k3a::name() }
        #endif
        #ifdef HAS_MKL
        , { RNGType::MKL_MCG31, kw::mkl_mcg31::name() }
        , { RNGType::MKL_R250, kw::mkl_r250::name() }
        , { RNGType::MKL_MRG32K3A, kw::mkl_mrg32k3a::name() }
        , { RNGType::MKL_MCG59, kw::mkl_mcg59::name() }
        , { RNGType::MKL_WH, kw::mkl_wh::name() }
        , { RNGType::MKL_MT19937, kw::mkl_mt19937::name() }
        , { RNGType::MKL_MT2203, kw::mkl_mt2203::name() }
        , { RNGType::MKL_SFMT19937, kw::mkl_sfmt19937::name() }
        , { RNGType::MKL_SOBOL, kw::mkl_sobol::name() }
        , { RNGType::MKL_NIEDERR, kw::mkl_niederr::name() }
        //, { RNGType::MKL_IABSTRACT, kw::mkl_iabstract::name() }
        //, { RNGType::MKL_DABSTRACT, kw::mkl_dabstract::name() }
        //, { RNGType::MKL_SABSTRACT, kw::mkl_sabstract::name() }
        , { RNGType::MKL_NONDETERM, kw::mkl_nondeterm::name() }
        #endif
        },
        //! keywords -> Enums
        { { "no_rng", RNGType::NO_RNG }
        , { kw::r123_threefry::string(), RNGType::R123_THREEFRY }
        , { kw::r123_philox::string(), RNGType::R123_PHILOX }
        #ifdef HAS_RNGSSE2
        , { kw::rngsse_gm19::string(), RNGType::RNGSSE_GM19 }
        , { kw::rngsse_gm29::string(), RNGType::RNGSSE_GM29 }
        , { kw::rngsse_gm31::string(), RNGType::RNGSSE_GM31 }
        , { kw::rngsse_gm55::string(), RNGType::RNGSSE_GM55 }
        , { kw::rngsse_gm61::string(), RNGType::RNGSSE_GM61 }
        , { kw::rngsse_gq581::string(), RNGType::RNGSSE_GQ581 }
        , { kw::rngsse_gq583::string(), RNGType::RNGSSE_GQ583 }
        , { kw::rngsse_gq584::string(), RNGType::RNGSSE_GQ584 }
        , { kw::rngsse_mt19937::string(), RNGType::RNGSSE_MT19937 }
        , { kw::rngsse_lfsr113::string(), RNGType::RNGSSE_LFSR113 }
        , { kw::rngsse_mrg32k3a::string(), RNGType::RNGSSE_MRG32K3A }
        #endif
        #ifdef HAS_MKL
        , { kw::mkl_mcg31::string(), RNGType::MKL_MCG31 }
        , { kw::mkl_r250::string(), RNGType::MKL_R250 }
        , { kw::mkl_mrg32k3a::string(), RNGType::MKL_MRG32K3A }
        , { kw::mkl_mcg59::string(), RNGType::MKL_MCG59 }
        , { kw::mkl_wh::string(), RNGType::MKL_WH }
        , { kw::mkl_mt19937::string(), RNGType::MKL_MT19937 }
        , { kw::mkl_mt2203::string(), RNGType::MKL_MT2203 }
        , { kw::mkl_sfmt19937::string(), RNGType::MKL_SFMT19937 }
        , { kw::mkl_sobol::string(), RNGType::MKL_SOBOL }
        , { kw::mkl_niederr::string(), RNGType::MKL_NIEDERR }
        //, { kw::mkl_iabstract::string(), RNGType::MKL_IABSTRACT }
        //, { kw::mkl_dabstract::string(), RNGType::MKL_DABSTRACT }
        //, { kw::mkl_sabstract::string(), RNGType::MKL_SABSTRACT }
        , { kw::mkl_nondeterm::string(), RNGType::MKL_NONDETERM }
        #endif
        } ) {}

    //! \brief Return parameter based on Enum
    //! \details Here 'parameter' is the library-specific identifier of the
    //!    option, i.e., as the library identifies the given option
    //! \param[in] rng Enum value of the option requested
    //! \return Library-specific parameter of the option
    //! \author J. Bakosi
    const ParamType& param( RNGType rng ) const {
      using tk::operator<<;
      auto it = brng.find( rng );
      Assert( it != end(brng),
              std::string("Cannot find parameter for RNG \"") << rng << "\"" );
      return it->second;
    }
 
    //! \brief Return field from RNG parameters bundle
    //! \details If user has specified it, return field. If user did not
    //!   specify it, return default.
    //! \param[in] rng Enum value of the option
    //! \param[in] def Default value of requested field if field is not found
    //! \param[in] bundle Parameter bundle to search in
    //! \return Field for requested enum if found, default if not found
    //! \author J. Bakosi
    template< class tag, class Param, class Field >
    Field param( RNGType rng, const Field& def, const Param& bundle ) const {
      auto it = bundle.find( rng );
      if ( it != bundle.end() ) return it->second.template get< tag >();
      else return def;
    }

    //! \brief Return RNG library type based on RNG options enum
    //! \param[in] rng Enum value of the option
    //! \return Library type enum
    //! \see tk::ctr::RNGLibType
    //! \author J. Bakosi
    RNGLibType lib( RNGType rng ) const {
      const auto& n = name( rng );
      if ( found( "MKL", n ) ) return RNGLibType::MKL;
      else if ( found( "RNGSSE", n ) ) return RNGLibType::RNGSSE;
      else if ( found( "PRAND", n ) ) return RNGLibType::PRAND;
      else if ( found( "Random123", n ) ) return RNGLibType::R123;
      else return RNGLibType::NO_LIB;
    }

    //! \brief Return whether RNG supports sequence option
    //! \param[in] rng Enum value of the option
    //! \return True if RNG supports sequence option
    //! \author J. Bakosi
    bool supportsSeq( RNGType rng ) const
    { return support.find( rng ) != end( support ) ? true : false; }

    //! \brief Return whether RNG supports sequence option given
    //! \param[in] rng Enum value of the option
    //! \param[in] option Option type
    //! \return True if RNG supports sequence option
    //! \author J. Bakosi
    template< class OptionType >
    bool supportsOpt( RNGType rng, const OptionType& option ) const {
      auto it = support.find( rng );
      if ( it != end( support ) ) {
        for (auto& o : it->second)
          if (o == option) return true;
      }
      return false;
    }

  private:
    //! \brief Search for 'kw' in 'str'
    //! \param[in] kw Keyword to search for
    //! \param[in] str String to search in
    //! \return True if found, false if not
    //! \author J. Bakosi
    bool found( const std::string& kw, const std::string& str ) const
    { return str.find( kw ) != std::string::npos ? true : false; }

    //! Enums -> MKL VSL BRNG parameters
    std::map< RNGType, ParamType > brng {
        { RNGType::NO_RNG, -1 }
      , { RNGType::R123_THREEFRY, 0 }
      , { RNGType::R123_PHILOX, 1 }
      #ifdef HAS_RNGSSE2
      , { RNGType::RNGSSE_GM19, 2 }
      , { RNGType::RNGSSE_GM29, 3 }
      , { RNGType::RNGSSE_GM31, 4 }
      , { RNGType::RNGSSE_GM55, 5 }
      , { RNGType::RNGSSE_GM61, 6 }
      , { RNGType::RNGSSE_GQ581, 7 }
      , { RNGType::RNGSSE_GQ583, 8 }
      , { RNGType::RNGSSE_GQ584, 9 }
      , { RNGType::RNGSSE_MT19937, 10 }
      , { RNGType::RNGSSE_LFSR113, 11 }
      , { RNGType::RNGSSE_MRG32K3A, 12 }
      #endif
      #ifdef HAS_MKL
      , { RNGType::MKL_MCG31, VSL_BRNG_MCG31 }
      , { RNGType::MKL_R250, VSL_BRNG_R250 }
      , { RNGType::MKL_MRG32K3A, VSL_BRNG_MRG32K3A }
      , { RNGType::MKL_MCG59, VSL_BRNG_MCG59 }
      , { RNGType::MKL_WH, VSL_BRNG_WH }
      , { RNGType::MKL_MT19937, VSL_BRNG_MT19937 }
      , { RNGType::MKL_MT2203, VSL_BRNG_MT2203 }
      , { RNGType::MKL_SFMT19937, VSL_BRNG_SFMT19937 }
      , { RNGType::MKL_SOBOL, VSL_BRNG_SOBOL }
      , { RNGType::MKL_NIEDERR, VSL_BRNG_NIEDERR }
      //, { RNGType::MKL_IABSTRACT, VSL_BRNG_IABSTRACT }
      //, { RNGType::MKL_DABSTRACT, VSL_BRNG_DABSTRACT }
      //, { RNGType::MKL_SABSTRACT, VSL_BRNG_SABSTRACT }
      , { RNGType::MKL_NONDETERM, VSL_BRNG_NONDETERM }
      #endif
    };

    //! Enums -> sequence length options supported
    std::map< RNGType, std::vector< RNGSSESeqLenType > > support {
      #ifdef HAS_RNGSSE2
      { RNGType::RNGSSE_GM29,    { RNGSSESeqLenType::SHORT,
                                   RNGSSESeqLenType::MEDIUM,
                                   RNGSSESeqLenType::LONG } },
      { RNGType::RNGSSE_GM31,    { RNGSSESeqLenType::SHORT,
                                   RNGSSESeqLenType::MEDIUM,
                                   RNGSSESeqLenType::LONG } },
      { RNGType::RNGSSE_GM55,    { RNGSSESeqLenType::SHORT,
                                   RNGSSESeqLenType::LONG } },
      { RNGType::RNGSSE_GQ581,   { RNGSSESeqLenType::SHORT,
                                   RNGSSESeqLenType::MEDIUM,
                                   RNGSSESeqLenType::LONG } },
      { RNGType::RNGSSE_GQ583,   { RNGSSESeqLenType::SHORT,
                                   RNGSSESeqLenType::MEDIUM,
                                   RNGSSESeqLenType::LONG } },
      { RNGType::RNGSSE_GQ584,   { RNGSSESeqLenType::SHORT,
                                   RNGSSESeqLenType::MEDIUM,
                                   RNGSSESeqLenType::LONG } },
      { RNGType::RNGSSE_GM61,    { RNGSSESeqLenType::SHORT,
                                   RNGSSESeqLenType::LONG } },
      { RNGType::RNGSSE_LFSR113, { RNGSSESeqLenType::SHORT,
                                   RNGSSESeqLenType::LONG } }
      #endif
    };
};

} // ctr::
} // tk::

#endif // RNGOptions_h
