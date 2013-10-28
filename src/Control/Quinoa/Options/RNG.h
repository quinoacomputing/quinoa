//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/RNG.h
  \author    J. Bakosi
  \date      Mon Oct 28 07:55:40 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's random number generator options and associations
  \details   Quinoa's random number generator options and associations
*/
//******************************************************************************
#ifndef QuinoaRNGOptions_h
#define QuinoaRNGOptions_h

#include <map>

#include <Config.h>

#ifdef HAS_MKL
#include <mkl_vsl.h>
#endif

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>
#include <RNG.h>

namespace quinoa {
namespace ctr {

//! Random number generator types
enum class RNGType : uint8_t { NO_RNG=0
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
                             , MKL_IABSTRACT
                             , MKL_DABSTRACT
                             , MKL_SABSTRACT
                             , MKL_NONDETERM
                             #endif
};

//! Random number generator library types
enum class RNGLibType : uint8_t { NO_LIB=0,
                                  MKL,
                                  RNGSSELIB,
                                  PRAND };

//! Random number generator factory type
using RNGFactory = std::map< RNGType, std::function<tk::RNG*()> >;

//! Class with base templated on the above enum class with associations
class RNG : public tk::Toggle<RNGType> {

  public:
    using ParamType = int;
    using LibType = RNGLibType;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit RNG() :
      Toggle<RNGType>("Random number generator", names, values) {}

    //! Return parameter based on Enum
    const ParamType& param(RNGType rng) const;

    //! Register random number generators into factory
    void initFactory(RNGFactory& f, int nthreads, unsigned int seed) const;
 
  private:
    //! Don't permit copy constructor
    RNG(const RNG&) = delete;
    //! Don't permit copy assigment
    RNG& operator=(const RNG&) = delete;
    //! Don't permit move constructor
    RNG(RNG&&) = delete;
    //! Don't permit move assigment
    RNG& operator=(RNG&&) = delete;

    //! Return RNG library type based on Enum
    RNGLibType lib(RNGType rng) const;

    //! Search for 'kw' in 'str'
    //! \param[in]  kw   Keyword to search for
    //! \param[in]  str  String to search in
    //! \return     True if found, false if not
    bool found(const std::string& kw, const std::string& str) const;

    //! Get access to RNG keywords
    const kw::mkl_mcg31 mkl_mcg31 {};
    const kw::mkl_r250 mkl_r250 {};
    const kw::mkl_mrg32k3a mkl_mrg32k3a {};
    const kw::mkl_mcg59 mkl_mcg59 {};
    const kw::mkl_wh mkl_wh {};
    const kw::mkl_mt19937 mkl_mt19937 {};
    const kw::mkl_mt2203 mkl_mt2203 {};
    const kw::mkl_sfmt19937 mkl_sfmt19937 {};
    const kw::mkl_sobol mkl_sobol {};
    const kw::mkl_niederr mkl_niederr {};
    const kw::mkl_iabstract mkl_iabstract {};
    const kw::mkl_dabstract mkl_dabstract {};
    const kw::mkl_sabstract mkl_sabstract {};
    const kw::mkl_nondeterm mkl_nondeterm {};

    //! Enums -> names
    const std::map<RNGType, std::string> names {
        { RNGType::NO_RNG, "n/a" }
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
      , { RNGType::MKL_IABSTRACT, mkl_iabstract.name() }
      , { RNGType::MKL_DABSTRACT, mkl_dabstract.name() }
      , { RNGType::MKL_SABSTRACT, mkl_sabstract.name() }
      , { RNGType::MKL_NONDETERM, mkl_nondeterm.name() }
      #endif
    };

    //! keywords -> Enums
    const std::map<std::string, RNGType> values {
        { "no_rng", RNGType::NO_RNG }
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
      , { mkl_iabstract.string(), RNGType::MKL_IABSTRACT }
      , { mkl_dabstract.string(), RNGType::MKL_DABSTRACT }
      , { mkl_sabstract.string(), RNGType::MKL_SABSTRACT }
      , { mkl_nondeterm.string(), RNGType::MKL_NONDETERM }
      #endif
    };

    //! Enums -> MKL VSL BRNG parameters
    const std::map<RNGType, ParamType> brng {
        { RNGType::NO_RNG, -1 }
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
      , { RNGType::MKL_IABSTRACT, VSL_BRNG_IABSTRACT }
      , { RNGType::MKL_DABSTRACT, VSL_BRNG_DABSTRACT }
      , { RNGType::MKL_SABSTRACT, VSL_BRNG_SABSTRACT }
      , { RNGType::MKL_NONDETERM, VSL_BRNG_NONDETERM }
      #endif
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaRNGOptions_h
