//******************************************************************************
/*!
  \file      src/Control/Options/MKLGaussianMethod.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 03:08:21 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Intel MKL Gaussian RNG method options
  \details   Intel MKL Gaussian RNG method options
*/
//******************************************************************************
#ifndef MKLGaussianMethodOptions_h
#define MKLGaussianMethodOptions_h

#include <map>

#include <mkl_vsl.h>

#include <Toggle.h>
#include <Keywords.h>

namespace tk {
namespace ctr {

//! MKL Gaussian random number generator methods
enum class MKLGaussianMethodType : uint8_t { BOXMULLER,
                                             BOXMULLER2,
                                             ICDF };

//! Pack/Unpack: delegate to tk::
inline void operator|( PUP::er& p, MKLGaussianMethodType& e )
{ PUP::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class MKLGaussianMethod : public tk::Toggle< MKLGaussianMethodType > {

  public:
    using ParamType = int;

    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MKLGaussianMethod() :
      Toggle< MKLGaussianMethodType >( "Gaussian method",
        //! Enums -> names
        { { MKLGaussianMethodType::BOXMULLER, kw::boxmuller().name() },
          { MKLGaussianMethodType::BOXMULLER2, kw::boxmuller2().name() },
          { MKLGaussianMethodType::ICDF, kw::icdf().name() } },
        //! keywords -> Enums
        { { kw::boxmuller().string(), MKLGaussianMethodType::BOXMULLER },
          { kw::boxmuller2().string(), MKLGaussianMethodType::BOXMULLER2 },
          { kw::icdf().string(), MKLGaussianMethodType::ICDF } } ) {}

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
