//******************************************************************************
/*!
  \file      src/Control/Options/RNG.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 03:07:36 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Quinoa's random number generator options and associations
  \details   Quinoa's random number generator options and associations
*/
//******************************************************************************
#ifndef RNGOptions_h
#define RNGOptions_h

#include <map>
#include <list>

#include <Config.h>

#ifdef HAS_MKL
#include <mkl_vsl.h>
#endif

#include <TaggedTuple.h>
#include <Toggle.h>
#include <Keywords.h>
#include <Options/RNGSSESeqLen.h>

namespace tk {
namespace ctr {

//! Random number generator types
enum class RNGType : uint8_t { NO_RNG=0
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
};

//! Pack/Unpack: delegate to tk::
inline void operator|( PUP::er& p, RNGType& e ) { PUP::pup( p, e ); }

//! Underlying type shortcut
using RawRNGType = std::underlying_type< RNGType >::type;

//! Return underlying type
constexpr RawRNGType raw( RNGType r ) { return static_cast< RawRNGType >( r ); }

//! Random number generator library types
enum class RNGLibType : uint8_t { NO_LIB=0,
                                  MKL,
                                  RNGSSE,
                                  PRAND };

//! Class with base templated on the above enum class with associations
class RNG : public tk::Toggle< RNGType > {

  public:
    using ParamType = int;
    using LibType = RNGLibType;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit RNG() :
      Toggle< RNGType >( "Random number generator",
        //! Enums -> names
        { { RNGType::NO_RNG, "n/a" }
        , { RNGType::RNGSSE_GM19, kw::rngsse_gm19().name() }
        , { RNGType::RNGSSE_GM29, kw::rngsse_gm29().name() }
        , { RNGType::RNGSSE_GM31, kw::rngsse_gm31().name() }
        , { RNGType::RNGSSE_GM55, kw::rngsse_gm55().name() }
        , { RNGType::RNGSSE_GM61, kw::rngsse_gm61().name() }
        , { RNGType::RNGSSE_GQ581, kw::rngsse_gq581().name() }
        , { RNGType::RNGSSE_GQ583, kw::rngsse_gq583().name() }
        , { RNGType::RNGSSE_GQ584, kw::rngsse_gq584().name() }
        , { RNGType::RNGSSE_MT19937, kw::rngsse_mt19937().name() }
        , { RNGType::RNGSSE_LFSR113, kw::rngsse_lfsr113().name() }
        , { RNGType::RNGSSE_MRG32K3A, kw::rngsse_mrg32k3a().name() }
        #ifdef HAS_MKL
        , { RNGType::MKL_MCG31, kw::mkl_mcg31().name() }
        , { RNGType::MKL_R250, kw::mkl_r250().name() }
        , { RNGType::MKL_MRG32K3A, kw::mkl_mrg32k3a().name() }
        , { RNGType::MKL_MCG59, kw::mkl_mcg59().name() }
        , { RNGType::MKL_WH, kw::mkl_wh().name() }
        , { RNGType::MKL_MT19937, kw::mkl_mt19937().name() }
        , { RNGType::MKL_MT2203, kw::mkl_mt2203().name() }
        , { RNGType::MKL_SFMT19937, kw::mkl_sfmt19937().name() }
        , { RNGType::MKL_SOBOL, kw::mkl_sobol().name() }
        , { RNGType::MKL_NIEDERR, kw::mkl_niederr().name() }
        //, { RNGType::MKL_IABSTRACT, kw::mkl_iabstract().name() }
        //, { RNGType::MKL_DABSTRACT, kw::mkl_dabstract().name() }
        //, { RNGType::MKL_SABSTRACT, kw::mkl_sabstract().name() }
        , { RNGType::MKL_NONDETERM, kw::mkl_nondeterm().name() }
        #endif
        },
        //! keywords -> Enums
        { { "no_rng", RNGType::NO_RNG }
        , { kw::rngsse_gm19().string(), RNGType::RNGSSE_GM19 }
        , { kw::rngsse_gm29().string(), RNGType::RNGSSE_GM29 }
        , { kw::rngsse_gm31().string(), RNGType::RNGSSE_GM31 }
        , { kw::rngsse_gm55().string(), RNGType::RNGSSE_GM55 }
        , { kw::rngsse_gm61().string(), RNGType::RNGSSE_GM61 }
        , { kw::rngsse_gq581().string(), RNGType::RNGSSE_GQ581 }
        , { kw::rngsse_gq583().string(), RNGType::RNGSSE_GQ583 }
        , { kw::rngsse_gq584().string(), RNGType::RNGSSE_GQ584 }
        , { kw::rngsse_mt19937().string(), RNGType::RNGSSE_MT19937 }
        , { kw::rngsse_lfsr113().string(), RNGType::RNGSSE_LFSR113 }
        , { kw::rngsse_mrg32k3a().string(), RNGType::RNGSSE_MRG32K3A }
        #ifdef HAS_MKL
        , { kw::mkl_mcg31().string(), RNGType::MKL_MCG31 }
        , { kw::mkl_r250().string(), RNGType::MKL_R250 }
        , { kw::mkl_mrg32k3a().string(), RNGType::MKL_MRG32K3A }
        , { kw::mkl_mcg59().string(), RNGType::MKL_MCG59 }
        , { kw::mkl_wh().string(), RNGType::MKL_WH }
        , { kw::mkl_mt19937().string(), RNGType::MKL_MT19937 }
        , { kw::mkl_mt2203().string(), RNGType::MKL_MT2203 }
        , { kw::mkl_sfmt19937().string(), RNGType::MKL_SFMT19937 }
        , { kw::mkl_sobol().string(), RNGType::MKL_SOBOL }
        , { kw::mkl_niederr().string(), RNGType::MKL_NIEDERR }
        //, { kw::mkl_iabstract().string(), RNGType::MKL_IABSTRACT }
        //, { kw::mkl_dabstract().string(), RNGType::MKL_DABSTRACT }
        //, { kw::mkl_sabstract().string(), RNGType::MKL_SABSTRACT }
        , { kw::mkl_nondeterm().string(), RNGType::MKL_NONDETERM }
        #endif
        } ) {}

    //! Return parameter based on Enum
    const ParamType& param( RNGType rng ) const {
      auto it = brng.find( rng );
      Assert( it != end(brng),
              std::string("Cannot find parameter for RNG \"") << rng << "\"" );
      return it->second;
    }
 
    //! Return field from RNG parameters bundle: if user has specified it,
    //! return it, if user did not specify it, return default
    template< class tag, class Param, class Field >
    Field param( RNGType rng, const Field& def, const Param& bundle ) const {
      auto it = bundle.find( rng );
      if ( it != bundle.end() ) return it->second.template get< tag >();
      else return def;
    }

    //! Return RNG library type based on Enum
    RNGLibType lib( RNGType rng ) const {
      const auto& n = name( rng );
      if ( found( "MKL", n ) ) return RNGLibType::MKL;
      else if ( found( "RNGSSE", n ) ) return RNGLibType::RNGSSE;
      else if ( found( "PRAND", n) ) return RNGLibType::PRAND;
      else return RNGLibType::NO_LIB;
    }

    //! Return whether RNG supports sequence option
    bool supportsSeq( RNGType rng ) const
    { return support.find( rng ) != end( support ) ? true : false; }

    //! Return whether RNG supports sequence option given
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
    //! Search for 'kw' in 'str'
    //! \param[in]  kw   Keyword to search for
    //! \param[in]  str  String to search in
    //! \return     True if found, false if not
    bool found(const std::string& kw, const std::string& str) const
    { return str.find( kw ) != std::string::npos ? true : false; }

    //! Enums -> MKL VSL BRNG parameters
    std::map< RNGType, ParamType > brng {
        { RNGType::NO_RNG, -1 }
      , { RNGType::RNGSSE_GM19, 0 }
      , { RNGType::RNGSSE_GM29, 1 }
      , { RNGType::RNGSSE_GM31, 2 }
      , { RNGType::RNGSSE_GM55, 3 }
      , { RNGType::RNGSSE_GM61, 4 }
      , { RNGType::RNGSSE_GQ581, 5 }
      , { RNGType::RNGSSE_GQ583, 6 }
      , { RNGType::RNGSSE_GQ584, 7 }
      , { RNGType::RNGSSE_MT19937, 8 }
      , { RNGType::RNGSSE_LFSR113, 9 }
      , { RNGType::RNGSSE_MRG32K3A, 10 }
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
    };
};

} // ctr::
} // tk::

#endif // RNGOptions_h
