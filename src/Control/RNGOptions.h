//******************************************************************************
/*!
  \file      src/Control/RNGOptions.h
  \author    J. Bakosi
  \date      Fri Sep 27 09:00:46 2013
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
#include <QuinoaKeywords.h>

namespace quinoa {
namespace sel {

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

//! Class with base templated on the above enum class with associations
class RNG : public quinoa::sel::Toggle<RNGType> {

  public:
    using ParamType = int;
    using LibType = RNGLibType;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit RNG() :
      quinoa::sel::Toggle<RNGType>("Random number generator", names, values) {}

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

    //! Get access to RNG keywords
    const grm::kw::mkl_mcg31 mkl_mcg31 {};
    const grm::kw::mkl_r250 mkl_r250 {};
    const grm::kw::mkl_mrg32k3a mkl_mrg32k3a {};
    const grm::kw::mkl_mcg59 mkl_mcg59 {};
    const grm::kw::mkl_wh mkl_wh {};
    const grm::kw::mkl_mt19937 mkl_mt19937 {};
    const grm::kw::mkl_mt2203 mkl_mt2203 {};
    const grm::kw::mkl_sfmt19937 mkl_sfmt19937 {};
    const grm::kw::mkl_sobol mkl_sobol {};
    const grm::kw::mkl_niederr mkl_niederr {};
    const grm::kw::mkl_iabstract mkl_iabstract {};
    const grm::kw::mkl_dabstract mkl_dabstract {};
    const grm::kw::mkl_sabstract mkl_sabstract {};
    const grm::kw::mkl_nondeterm mkl_nondeterm {};

    //! Enums -> names
    const std::map<RNGType, std::string> names {
      { RNGType::NO_RNG, "n/a" },
      { RNGType::MKL_MCG31, mkl_mcg31.name() },
      { RNGType::MKL_R250, mkl_r250.name() },
      { RNGType::MKL_MRG32K3A, mkl_mrg32k3a.name() },
      { RNGType::MKL_MCG59, mkl_mcg59.name() },
      { RNGType::MKL_WH, mkl_wh.name() },
      { RNGType::MKL_MT19937, mkl_mt19937.name() },
      { RNGType::MKL_MT2203, mkl_mt2203.name() },
      { RNGType::MKL_SFMT19937, mkl_sfmt19937.name() },
      { RNGType::MKL_SOBOL, mkl_sobol.name() },
      { RNGType::MKL_NIEDERR, mkl_niederr.name() },
      { RNGType::MKL_IABSTRACT, mkl_iabstract.name() },
      { RNGType::MKL_DABSTRACT, mkl_dabstract.name() },
      { RNGType::MKL_SABSTRACT, mkl_sabstract.name() },
      { RNGType::MKL_NONDETERM, mkl_nondeterm.name() }
    };

    //! keywords -> Enums
    const std::map<std::string, RNGType> values {
      { "no_rng", RNGType::NO_RNG },
      { mkl_mcg31.string(), RNGType::MKL_MCG31 },
      { mkl_r250.string(), RNGType::MKL_R250 },
      { mkl_mrg32k3a.string(), RNGType::MKL_MRG32K3A },
      { mkl_mcg59.string(), RNGType::MKL_MCG59 },
      { mkl_wh.string(), RNGType::MKL_WH },
      { mkl_mt19937.string(), RNGType::MKL_MT19937 },
      { mkl_mt2203.string(), RNGType::MKL_MT2203 },
      { mkl_sfmt19937.string(), RNGType::MKL_SFMT19937 },
      { mkl_sobol.string(), RNGType::MKL_SOBOL },
      { mkl_niederr.string(), RNGType::MKL_NIEDERR },
      { mkl_iabstract.string(), RNGType::MKL_IABSTRACT },
      { mkl_dabstract.string(), RNGType::MKL_DABSTRACT },
      { mkl_sabstract.string(), RNGType::MKL_SABSTRACT },
      { mkl_nondeterm.string(), RNGType::MKL_NONDETERM }
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

} // sel::
} // quinoa::

#endif // RNGOptions_h
