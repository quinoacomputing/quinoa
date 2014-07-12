//******************************************************************************
/*!
  \file      src/Control/Options/MKLGaussianMethod.h
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 08:58:39 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Intel MKL Gaussian RNG method options
  \details   Intel MKL Gaussian RNG method options
*/
//******************************************************************************
#ifndef MKLGaussianMethodOptions_h
#define MKLGaussianMethodOptions_h

#include <map>

#include <mkl_vsl.h>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace tk {
namespace ctr {

//! MKL Gaussian random number generator methods
enum class MKLGaussianMethodType : uint8_t { BOXMULLER,
                                             BOXMULLER2,
                                             ICDF };

//! Class with base templated on the above enum class with associations
class MKLGaussianMethod : public tk::Toggle< MKLGaussianMethodType > {

  public:
    using ParamType = int;

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MKLGaussianMethod() :
      Toggle< MKLGaussianMethodType >("Gaussian method", names, values) {}

    //! Return parameter based on Enum
    const ParamType& param(MKLGaussianMethodType rng) const;

  private:
    //! Don't permit copy constructor
    MKLGaussianMethod(const MKLGaussianMethod&) = delete;
    //! Don't permit copy assigment
    MKLGaussianMethod& operator=(const MKLGaussianMethod&) = delete;
    //! Don't permit move constructor
    MKLGaussianMethod(MKLGaussianMethod&&) = delete;
    //! Don't permit move assigment
    MKLGaussianMethod& operator=(MKLGaussianMethod&&) = delete;

    //! Get access to MKL Gaussian method keywords
    const tk::kw::boxmuller boxmuller{};
    const tk::kw::boxmuller2 boxmuller2 {};
    const tk::kw::icdf icdf {};

    //! Enums -> names
    const std::map< MKLGaussianMethodType, std::string > names {
      { MKLGaussianMethodType::BOXMULLER, boxmuller.name() },
      { MKLGaussianMethodType::BOXMULLER2, boxmuller2.name() },
      { MKLGaussianMethodType::ICDF, icdf.name() }
    };

    //! keywords -> Enums
    const std::map< std::string, MKLGaussianMethodType > values {
      { boxmuller.string(), MKLGaussianMethodType::BOXMULLER },
      { boxmuller2.string(), MKLGaussianMethodType::BOXMULLER2 },
      { icdf.string(), MKLGaussianMethodType::ICDF }
    };

    //! Enums -> MKL VSL RNG GAUSSIAN METHOD parameters
    const std::map< MKLGaussianMethodType, ParamType > method {
      { MKLGaussianMethodType::BOXMULLER, VSL_RNG_METHOD_GAUSSIAN_BOXMULLER },
      { MKLGaussianMethodType::BOXMULLER2, VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2 },
      { MKLGaussianMethodType::ICDF, VSL_RNG_METHOD_GAUSSIAN_ICDF }
    };
};

} // ctr::
} // tk::

#endif // MKLGaussianMethodOptions_h
