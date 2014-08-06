//******************************************************************************
/*!
  \file      src/Control/Options/MKLGaussianMethod.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 10:42:41 AM MDT
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

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MKLGaussianMethod() :
      Toggle< MKLGaussianMethodType >( "Gaussian method",
        //! Enums -> names
        { { MKLGaussianMethodType::BOXMULLER, tk::kw::boxmuller().name() },
          { MKLGaussianMethodType::BOXMULLER2, tk::kw::boxmuller2().name() },
          { MKLGaussianMethodType::ICDF, tk::kw::icdf().name() } },
        //! keywords -> Enums
        { { tk::kw::boxmuller().string(), MKLGaussianMethodType::BOXMULLER },
          { tk::kw::boxmuller2().string(), MKLGaussianMethodType::BOXMULLER2 },
          { tk::kw::icdf().string(), MKLGaussianMethodType::ICDF } } ) {}

    //! Return parameter based on Enum
    const ParamType& param( MKLGaussianMethodType m ) const {
      auto it = method.find( m );
      Assert( it != end(method),
              std::string("Cannot find parameter for MKLGaussianMethod \"")
              << m << "\"" );
      return it->second;
    }

  private:
    //! Enums -> MKL VSL RNG GAUSSIAN METHOD parameters
    std::map< MKLGaussianMethodType, ParamType > method {
      { MKLGaussianMethodType::BOXMULLER, VSL_RNG_METHOD_GAUSSIAN_BOXMULLER },
      { MKLGaussianMethodType::BOXMULLER2, VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2 },
      { MKLGaussianMethodType::ICDF, VSL_RNG_METHOD_GAUSSIAN_ICDF }
    };
};

} // ctr::
} // tk::

#endif // MKLGaussianMethodOptions_h
