//******************************************************************************
/*!
  \file      src/Control/RNGOptions.h
  \author    J. Bakosi
  \date      Thu Aug  1 14:40:41 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator options and associations
  \details   Random number generator options and associations
*/
//******************************************************************************
#ifndef RNGOptions_h
#define RNGOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace Quinoa {

namespace select {

//! Random number generator test types
enum class RNGTypes : uint8_t { NO_RNG=0,
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

//! Class with base templated on the above enum class with associations
class RNG : public Toggle<RNGTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    RNG() : Toggle<RNGTypes>(names, values) {
      //! Enums -> names
      names[RNGTypes::NO_RNG] = "No RNG";
      names[RNGTypes::MKL_MCG31] = "MKL_VSL_MCG31";
      names[RNGTypes::MKL_R250] = "MKL_VSL_R250";
      names[RNGTypes::MKL_MRG32K3A] = "MKL_VSL_MRG32K3A";
      names[RNGTypes::MKL_MCG59] = "MKL_VSL_MCG59";
      names[RNGTypes::MKL_WH] = "MKL_VSL_WH";
      names[RNGTypes::MKL_MT19937] = "MKL_VSL_MT19937";
      names[RNGTypes::MKL_MT2203] = "MKL_VSL_MT2203";
      names[RNGTypes::MKL_SFMT19937] = "MKL_VSL_SFMT19937";
      names[RNGTypes::MKL_SOBOL] = "MKL_VSL_SOBOL";
      names[RNGTypes::MKL_NIEDERR] = "MKL_VSL_NIEDERR";
      names[RNGTypes::MKL_IABSTRACT] = "MKL_VSL_IABSTRACT";
      names[RNGTypes::MKL_DABSTRACT] = "MKL_VSL_DABSTRACT";
      names[RNGTypes::MKL_SABSTRACT] = "MKL_VSL_SABSTRACT";
      names[RNGTypes::MKL_NONDETERM] = "MKL_VSL_NONDETERM";
      //! keywords -> Enums
      values["no_rng"] = RNGTypes::NO_RNG;
      values["mkl_mcg31"] = RNGTypes::MKL_MCG31;
      values["mkl_r250"] = RNGTypes::MKL_R250;
      values["mkl_mrg32k3a"] = RNGTypes::MKL_MRG32K3A;
      values["mkl_mcg59"] = RNGTypes::MKL_MCG59;
      values["mkl_wh"] = RNGTypes::MKL_WH;
      values["mkl_mt19937"] = RNGTypes::MKL_MT19937;
      values["mkl_mt2203"] = RNGTypes::MKL_MT2203;
      values["mkl_sfmt19937"] = RNGTypes::MKL_SFMT19937;
      values["mkl_sobol"] = RNGTypes::MKL_SOBOL;
      values["mkl_niederr"] = RNGTypes::MKL_NIEDERR;
      values["mkl_iabstract"] = RNGTypes::MKL_IABSTRACT;
      values["mkl_dabstract"] = RNGTypes::MKL_DABSTRACT;
      values["mkl_sabstract"] = RNGTypes::MKL_SABSTRACT;
      values["mkl_nondeterm"] = RNGTypes::MKL_NONDETERM;
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

    std::map<RNGTypes, std::string> names;
    std::map<std::string, RNGTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // RNGOptions_h
