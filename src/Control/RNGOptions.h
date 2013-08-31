//******************************************************************************
/*!
  \file      src/Control/RNGOptions.h
  \author    J. Bakosi
  \date      Fri 30 Aug 2013 06:23:53 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator options and associations
  \details   Random number generator options and associations
*/
//******************************************************************************
#ifndef RNGOptions_h
#define RNGOptions_h

#include <map>

#include <mkl_vsl.h>

#include <Exception.h>
#include <Toggle.h>

namespace rngtest {

namespace select {

//! Random number generator test types
enum class RNGType : uint8_t { NO_RNG=0,
                               MKL_MCG31,
                               MKL_R250,
                               MKL_MRG32K3A,
                               MKL_MCG59,
                               MKL_WH,
                               MKL_MT19937,
                               MKL_MT2203,
                               MKL_SFMT19937,
                               MKL_SOBOL,
                               MKL_NIEDERR,
                               MKL_IABSTRACT,
                               MKL_DABSTRACT,
                               MKL_SABSTRACT,
                               MKL_NONDETERM };

//! Random number generator library types
enum class RNGLibType : uint8_t { NO_LIB=0,
                                  MKL,
                                  RNGSSELIB,
                                  PRAND };

using quinoa::select::operator+;
using quinoa::select::operator<<;

//! Class with base templated on the above enum class with associations
class RNG : public quinoa::select::Toggle<RNGType> {

  public:
    using ParamType = int;
    using LibType = RNGLibType;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit RNG() : quinoa::select::Toggle<RNGType>(names, values) {}

  private:
    //! Don't permit copy constructor
    RNG(const RNG&) = delete;
    //! Don't permit copy assigment
    RNG& operator=(const RNG&) = delete;
    //! Don't permit move constructor
    RNG(RNG&&) = delete;
    //! Don't permit move assigment
    RNG& operator=(RNG&&) = delete;

    //! Return parameter based on Enum
    const ParamType& param(RNGType rng) const {
      auto it = brng.find(rng);
      Assert(it != brng.end(), quinoa::ExceptType::FATAL,
             std::string("Cannot find parameter for RNG \"") + rng + "\"");
      return it->second;
    }

    //! Return RNG library type based on Enum
    RNGLibType lib(RNGType rng) const {
      auto it = names.find(rng);
      Assert(it != names.end(), quinoa::ExceptType::FATAL,
             std::string("Cannot find name for RNG \"") + rng + "\"");
      if (found("MKL", it->second)) return RNGLibType::MKL;
      else if (found("RNGSSELIB", it->second)) return RNGLibType::RNGSSELIB;
      else if (found("PRAND", it->second)) return RNGLibType::PRAND;
      else return RNGLibType::NO_LIB;
    }

    //! Search for 'kw' in 'str'
    //! \param[in]  kw   Keyword to search for
    //! \param[in]  str  String to search in
    //! \return     True if found, false if not
    bool found(const std::string& kw, const std::string& str) const {
      std::size_t f = str.find(kw);
      if (f != std::string::npos) return true; else return false;
    }

    //! Enums -> names
    const std::map<RNGType, std::string> names {
      { RNGType::NO_RNG, "No RNG" },
      { RNGType::MKL_MCG31, "MKL_VSL_MCG31" },
      { RNGType::MKL_R250, "MKL_VSL_R250" },
      { RNGType::MKL_MRG32K3A, "MKL_VSL_MRG32K3A" },
      { RNGType::MKL_MCG59, "MKL_VSL_MCG59" },
      { RNGType::MKL_WH, "MKL_VSL_WH" },
      { RNGType::MKL_MT19937, "MKL_VSL_MT19937" },
      { RNGType::MKL_MT2203, "MKL_VSL_MT2203" },
      { RNGType::MKL_SFMT19937, "MKL_VSL_SFMT19937" },
      { RNGType::MKL_SOBOL, "MKL_VSL_SOBOL" },
      { RNGType::MKL_NIEDERR, "MKL_VSL_NIEDERR" },
      { RNGType::MKL_IABSTRACT, "MKL_VSL_IABSTRACT" },
      { RNGType::MKL_DABSTRACT, "MKL_VSL_DABSTRACT" },
      { RNGType::MKL_SABSTRACT, "MKL_VSL_SABSTRACT" },
      { RNGType::MKL_NONDETERM, "MKL_VSL_NONDETERM" }
    };

    //! keywords -> Enums
    const std::map<std::string, RNGType> values {
      { "no_rng", RNGType::NO_RNG },
      { "mkl_mcg31", RNGType::MKL_MCG31 },
      { "mkl_r250", RNGType::MKL_R250 },
      { "mkl_mrg32k3a", RNGType::MKL_MRG32K3A },
      { "mkl_mcg59", RNGType::MKL_MCG59 },
      { "mkl_wh", RNGType::MKL_WH },
      { "mkl_mt19937", RNGType::MKL_MT19937 },
      { "mkl_mt2203", RNGType::MKL_MT2203 },
      { "mkl_sfmt19937", RNGType::MKL_SFMT19937 },
      { "mkl_sobol", RNGType::MKL_SOBOL },
      { "mkl_niederr", RNGType::MKL_NIEDERR },
      { "mkl_iabstract", RNGType::MKL_IABSTRACT },
      { "mkl_dabstract", RNGType::MKL_DABSTRACT },
      { "mkl_sabstract", RNGType::MKL_SABSTRACT },
      { "mkl_nondeterm", RNGType::MKL_NONDETERM }
    };

    //! Enums -> MKL VSL BRNG parameters
    const std::map<RNGType, ParamType> brng {
      { RNGType::NO_RNG, -1 },
      { RNGType::MKL_MCG31, VSL_BRNG_MCG31 },
      { RNGType::MKL_R250, VSL_BRNG_R250 },
      { RNGType::MKL_MRG32K3A, VSL_BRNG_MRG32K3A },
      { RNGType::MKL_MCG59, VSL_BRNG_MCG59 },
      { RNGType::MKL_WH, VSL_BRNG_WH },
      { RNGType::MKL_MT19937, VSL_BRNG_MT19937 },
      { RNGType::MKL_MT2203, VSL_BRNG_MT2203 },
      { RNGType::MKL_SFMT19937, VSL_BRNG_SFMT19937 },
      { RNGType::MKL_SOBOL, VSL_BRNG_SOBOL },
      { RNGType::MKL_NIEDERR, VSL_BRNG_NIEDERR },
      { RNGType::MKL_IABSTRACT, VSL_BRNG_IABSTRACT },
      { RNGType::MKL_DABSTRACT, VSL_BRNG_DABSTRACT },
      { RNGType::MKL_SABSTRACT, VSL_BRNG_SABSTRACT },
      { RNGType::MKL_NONDETERM, VSL_BRNG_NONDETERM }
    };
};

} // namespace select

} // namespace rngtest

#endif // RNGOptions_h
