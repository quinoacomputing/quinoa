//******************************************************************************
/*!
  \file      src/Control/Options/RNG.h
  \author    J. Bakosi
  \date      Tue 03 Jun 2014 07:15:46 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
#include <Quinoa/Tags.h>
#include <Quinoa/InputDeck/Keywords.h>
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
  // MKL VSL's abstract RNGs are for use with external (i.e., external to MKL or
  // user-defined generator via call-back functions), disabled for now
                             //, MKL_IABSTRACT
                             //, MKL_DABSTRACT
                             //, MKL_SABSTRACT
                             , MKL_NONDETERM
                             #endif
};

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
      Toggle< RNGType >( "Random number generator", names, values ) {}

    //! Return parameter based on Enum
    const ParamType& param(RNGType rng) const;
 
    //! Return field from RNG parameters bundle: if user has specified it,
    //! return it, if user did not specify it, return default
    template< class tag, class Param, class Field >
    Field param( RNGType rng, const Field& def, const Param& bundle ) const {
      auto it = bundle.find( rng );
      if ( it != bundle.end() ) return it->second.template get< tag >();
      else return def;
    }

    //! Return RNG library type based on Enum
    RNGLibType lib( RNGType rng ) const;

    //! Return whether RNG supports sequence option
    bool supportsSeq( RNGType rng ) const {
      auto it = support.find( rng );
      if ( it != support.end() ) return true;
      else return false;
    }

    //! Return whether RNG supports sequence option given
    template< class OptionType >
    bool supportsOpt( RNGType rng, const OptionType& option ) const {
      auto it = support.find( rng );
      if ( it != support.end() )
        for (auto& o : it->second)
          if (o == option) return true;
      return false;
    }

  private:
    //! Don't permit copy constructor
    RNG(const RNG&) = delete;
    //! Don't permit copy assigment
    RNG& operator=(const RNG&) = delete;
    //! Don't permit move constructor
    RNG(RNG&&) = delete;
    //! Don't permit move assigment
    RNG& operator=(RNG&&) = delete;

    //! Search for 'kw' in 'str'
    //! \param[in]  kw   Keyword to search for
    //! \param[in]  str  String to search in
    //! \return     True if found, false if not
    bool found(const std::string& kw, const std::string& str) const;

    //! Get access to RNG keywords
    const tk::kw::rngsse_gm19 rngsse_gm19 {};
    const tk::kw::rngsse_gm29 rngsse_gm29 {};
    const tk::kw::rngsse_gm31 rngsse_gm31 {};
    const tk::kw::rngsse_gm55 rngsse_gm55 {};
    const tk::kw::rngsse_gm61 rngsse_gm61 {};
    const tk::kw::rngsse_gq581 rngsse_gq581 {};
    const tk::kw::rngsse_gq583 rngsse_gq583 {};
    const tk::kw::rngsse_gq584 rngsse_gq584 {};
    const tk::kw::rngsse_mt19937 rngsse_mt19937 {};
    const tk::kw::rngsse_lfsr113 rngsse_lfsr113 {};
    const tk::kw::rngsse_mrg32k3a rngsse_mrg32k3a {};
    #ifdef HAS_MKL
    const tk::kw::mkl_mcg31 mkl_mcg31 {};
    const tk::kw::mkl_r250 mkl_r250 {};
    const tk::kw::mkl_mrg32k3a mkl_mrg32k3a {};
    const tk::kw::mkl_mcg59 mkl_mcg59 {};
    const tk::kw::mkl_wh mkl_wh {};
    const tk::kw::mkl_mt19937 mkl_mt19937 {};
    const tk::kw::mkl_mt2203 mkl_mt2203 {};
    const tk::kw::mkl_sfmt19937 mkl_sfmt19937 {};
    const tk::kw::mkl_sobol mkl_sobol {};
    const tk::kw::mkl_niederr mkl_niederr {};
    //const tk::kw::mkl_iabstract mkl_iabstract {};
    //const tk::kw::mkl_dabstract mkl_dabstract {};
    //const tk::kw::mkl_sabstract mkl_sabstract {};
    const tk::kw::mkl_nondeterm mkl_nondeterm {};
    #endif

    //! Enums -> names
    const std::map<RNGType, std::string> names {
        { RNGType::NO_RNG, "n/a" }
      , { RNGType::RNGSSE_GM19, rngsse_gm19.name() }
      , { RNGType::RNGSSE_GM29, rngsse_gm29.name() }
      , { RNGType::RNGSSE_GM31, rngsse_gm31.name() }
      , { RNGType::RNGSSE_GM55, rngsse_gm55.name() }
      , { RNGType::RNGSSE_GM61, rngsse_gm61.name() }
      , { RNGType::RNGSSE_GQ581, rngsse_gq581.name() }
      , { RNGType::RNGSSE_GQ583, rngsse_gq583.name() }
      , { RNGType::RNGSSE_GQ584, rngsse_gq584.name() }
      , { RNGType::RNGSSE_MT19937, rngsse_mt19937.name() }
      , { RNGType::RNGSSE_LFSR113, rngsse_lfsr113.name() }
      , { RNGType::RNGSSE_MRG32K3A, rngsse_mrg32k3a.name() }
      #ifdef HAS_MKL
      , { RNGType::MKL_MCG31, mkl_mcg31.name() }
      , { RNGType::MKL_R250, mkl_r250.name() }
      , { RNGType::MKL_MRG32K3A, mkl_mrg32k3a.name() }
      , { RNGType::MKL_MCG59, mkl_mcg59.name() }
      , { RNGType::MKL_WH, mkl_wh.name() }
      , { RNGType::MKL_MT19937, mkl_mt19937.name() }
      , { RNGType::MKL_MT2203, mkl_mt2203.name() }
      , { RNGType::MKL_SFMT19937, mkl_sfmt19937.name() }
      , { RNGType::MKL_SOBOL, mkl_sobol.name() }
      , { RNGType::MKL_NIEDERR, mkl_niederr.name() }
      //, { RNGType::MKL_IABSTRACT, mkl_iabstract.name() }
      //, { RNGType::MKL_DABSTRACT, mkl_dabstract.name() }
      //, { RNGType::MKL_SABSTRACT, mkl_sabstract.name() }
      , { RNGType::MKL_NONDETERM, mkl_nondeterm.name() }
      #endif
    };

    //! keywords -> Enums
    const std::map<std::string, RNGType> values {
        { "no_rng", RNGType::NO_RNG }
      , { rngsse_gm19.string(), RNGType::RNGSSE_GM19 }
      , { rngsse_gm29.string(), RNGType::RNGSSE_GM29 }
      , { rngsse_gm31.string(), RNGType::RNGSSE_GM31 }
      , { rngsse_gm55.string(), RNGType::RNGSSE_GM55 }
      , { rngsse_gm61.string(), RNGType::RNGSSE_GM61 }
      , { rngsse_gq581.string(), RNGType::RNGSSE_GQ581 }
      , { rngsse_gq583.string(), RNGType::RNGSSE_GQ583 }
      , { rngsse_gq584.string(), RNGType::RNGSSE_GQ584 }
      , { rngsse_mt19937.string(), RNGType::RNGSSE_MT19937 }
      , { rngsse_lfsr113.string(), RNGType::RNGSSE_LFSR113 }
      , { rngsse_mrg32k3a.string(), RNGType::RNGSSE_MRG32K3A }
      #ifdef HAS_MKL
      , { mkl_mcg31.string(), RNGType::MKL_MCG31 }
      , { mkl_r250.string(), RNGType::MKL_R250 }
      , { mkl_mrg32k3a.string(), RNGType::MKL_MRG32K3A }
      , { mkl_mcg59.string(), RNGType::MKL_MCG59 }
      , { mkl_wh.string(), RNGType::MKL_WH }
      , { mkl_mt19937.string(), RNGType::MKL_MT19937 }
      , { mkl_mt2203.string(), RNGType::MKL_MT2203 }
      , { mkl_sfmt19937.string(), RNGType::MKL_SFMT19937 }
      , { mkl_sobol.string(), RNGType::MKL_SOBOL }
      , { mkl_niederr.string(), RNGType::MKL_NIEDERR }
      //, { mkl_iabstract.string(), RNGType::MKL_IABSTRACT }
      //, { mkl_dabstract.string(), RNGType::MKL_DABSTRACT }
      //, { mkl_sabstract.string(), RNGType::MKL_SABSTRACT }
      , { mkl_nondeterm.string(), RNGType::MKL_NONDETERM }
      #endif
    };

    //! Enums -> MKL VSL BRNG parameters
    const std::map<RNGType, ParamType> brng {
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
    const std::map< RNGType, std::vector< RNGSSESeqLenType > > support {
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
