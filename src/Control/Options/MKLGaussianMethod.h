// *****************************************************************************
/*!
  \file      src/Control/Options/MKLGaussianMethod.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Intel MKL Gaussian RNG method options
  \details   Intel MKL Gaussian RNG method options
*/
// *****************************************************************************
#ifndef MKLGaussianMethodOptions_h
#define MKLGaussianMethodOptions_h

#include <map>

#include <boost/mpl/vector.hpp>

#include "NoWarning/mkl_vsl.h"

#include "Toggle.h"
#include "Keywords.h"
#include "PUPUtil.h"

namespace tk {
namespace ctr {

//! MKL Gaussian random number generator method types
//! \author J. Bakosi
enum class MKLGaussianMethodType : uint8_t { BOXMULLER,
                                             BOXMULLER2,
                                             ICDF };

//! \brief Pack/Unpack MKLGaussianMethodType: forward overload to generic enum
//!   class packer
//! \author J. Bakosi
inline void operator|( PUP::er& p, MKLGaussianMethodType& e )
{ PUP::pup( p, e ); }

//! \brief MKLGaussianMethod options: outsource searches to base templated on
//!   enum type
//! \author J. Bakosi
class MKLGaussianMethod : public tk::Toggle< MKLGaussianMethodType > {

  public:
    using ParamType = int;

    //! Valid expected choices to make them also available at compile-time
    //! \author J. Bakosi
    using keywords = boost::mpl::vector< kw::boxmuller
                                       , kw::boxmuller2
                                       , kw::icdf
                                       >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    //! \author J. Bakosi
    explicit MKLGaussianMethod() :
      tk::Toggle< MKLGaussianMethodType >(
        //! Group, i.e., options, name
        "Gaussian method",
        //! Enums -> names
        { { MKLGaussianMethodType::BOXMULLER, kw::boxmuller::name() },
          { MKLGaussianMethodType::BOXMULLER2, kw::boxmuller2::name() },
          { MKLGaussianMethodType::ICDF, kw::icdf::name() } },
        //! keywords -> Enums
        { { kw::boxmuller::string(), MKLGaussianMethodType::BOXMULLER },
          { kw::boxmuller2::string(), MKLGaussianMethodType::BOXMULLER2 },
          { kw::icdf::string(), MKLGaussianMethodType::ICDF } } ) {}

    //! \brief Return parameter based on Enum
    //! \details Here 'parameter' is the library-specific identifier of the
    //!    option, i.e., as the library identifies the given option
    //! \param[in] m Enum value of the option requested
    //! \return Library-specific parameter of the option
    //! \author J. Bakosi
    const ParamType& param( MKLGaussianMethodType m ) const {
      using tk::operator<<;
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
