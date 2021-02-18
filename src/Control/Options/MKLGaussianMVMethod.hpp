// *****************************************************************************
/*!
  \file      src/Control/Options/MKLGaussianMVMethod.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Intel MKL multi-variate Gaussian RNG method options
  \details   Intel MKL multi-variate Gaussian RNG method options
*/
// *****************************************************************************
#ifndef MKLGaussianMVMethodOptions_h
#define MKLGaussianMVMethodOptions_h

#include <map>

#include <mkl_vsl_defines.h>

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! MKL Gaussian random number generator method types
enum class MKLGaussianMVMethodType : uint8_t { BOXMULLER,
                                               BOXMULLER2,
                                               ICDF };

//! \brief Pack/Unpack MKLGaussianMVMethodType: forward overload to generic enum
//!   class packer
inline void operator|( PUP::er& p, MKLGaussianMVMethodType& e )
{ PUP::pup( p, e ); }

//! \brief MKLGaussianMVMethod options: outsource searches to base templated on
//!   enum type
class MKLGaussianMVMethod : public tk::Toggle< MKLGaussianMVMethodType > {

  public:
    using ParamType = int;

    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::boxmuller
                                  , kw::boxmuller2
                                  , kw::icdf
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit MKLGaussianMVMethod() :
      tk::Toggle< MKLGaussianMVMethodType >(
        //! Group, i.e., options, name
        "Gaussian method",
        //! Enums -> names
        { { MKLGaussianMVMethodType::BOXMULLER, kw::boxmuller::name() },
          { MKLGaussianMVMethodType::BOXMULLER2, kw::boxmuller2::name() },
          { MKLGaussianMVMethodType::ICDF, kw::icdf::name() } },
        //! keywords -> Enums
        { { kw::boxmuller::string(), MKLGaussianMVMethodType::BOXMULLER },
          { kw::boxmuller2::string(), MKLGaussianMVMethodType::BOXMULLER2 },
          { kw::icdf::string(), MKLGaussianMVMethodType::ICDF } } ) {}

    //! \brief Return parameter based on Enum
    //! \details Here 'parameter' is the library-specific identifier of the
    //!    option, i.e., as the library identifies the given option
    //! \param[in] m Enum value of the option requested
    //! \return Library-specific parameter of the option
    const ParamType& param( MKLGaussianMVMethodType m ) const {
      using tk::operator<<;
      auto it = method.find( m );
      Assert( it != end(method),
              std::string("Cannot find parameter for MKLGaussianMVMethod \"")
              << m << "\"" );
      return it->second;
    }

  private:
    //! Enums -> MKL VSL RNG GAUSSIAN METHOD parameters
    std::map< MKLGaussianMVMethodType, ParamType > method {
      { MKLGaussianMVMethodType::BOXMULLER,
        VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER },
      { MKLGaussianMVMethodType::BOXMULLER2,
        VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2 },
      { MKLGaussianMVMethodType::ICDF, VSL_RNG_METHOD_GAUSSIANMV_ICDF }
    };
};

} // ctr::
} // tk::

#endif // MKLGaussianMVMethodOptions_h
