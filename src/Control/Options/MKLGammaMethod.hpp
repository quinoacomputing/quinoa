// *****************************************************************************
/*!
  \file      src/Control/Options/MKLGammaMethod.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Intel MKL Gamma RNG method options
  \details   Intel MKL Gamma RNG method options
*/
// *****************************************************************************
#ifndef MKLGammaMethodOptions_h
#define MKLGammaMethodOptions_h

#include <map>

#include <mkl_vsl_defines.h>

#include <brigand/sequences/list.hpp>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! MKL gamma random number generator method types
enum class MKLGammaMethodType : uint8_t { GNORM,
                                          GNORM_ACCURATE };

//! \brief Pack/Unpack MKLGammaMethodType: forward overload to generic enum
//!   class packer
inline void operator|( PUP::er& p, MKLGammaMethodType& e ) { PUP::pup( p, e ); }

//! \brief MKLGammaMethod options: outsource searches to base templated on
//!   enum type
class MKLGammaMethod : public tk::Toggle< MKLGammaMethodType > {

  public:
    using ParamType = int;

    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::gnorm
                                  , kw::gnorm_accurate
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit MKLGammaMethod() :
      tk::Toggle< MKLGammaMethodType >(
        //! Group, i.e., options, name
        "Gamma method",
        //! Enums -> names
        { { MKLGammaMethodType::GNORM, kw::gnorm::name() },
          { MKLGammaMethodType::GNORM_ACCURATE, kw::gnorm_accurate::name() } },
        //! keywords -> Enums
        { { kw::gnorm::string(), MKLGammaMethodType::GNORM },
          { kw::gnorm_accurate::string(),
            MKLGammaMethodType::GNORM_ACCURATE } } ) {}

    //! \brief Return parameter based on Enum
    //! \details Here 'parameter' is the library-specific identifier of the
    //!    option, i.e., as the library identifies the given option
    //! \param[in] m Enum value of the option requested
    //! \return Library-specific parameter of the option
    const ParamType& param( MKLGammaMethodType m ) const {
      using tk::operator<<;
      auto it = method.find( m );
      Assert( it != end(method),
              std::string("Cannot find parameter for MKLGammaMethod \"")
              << m << "\"" );
      return it->second;
    }

  private:
    //! Enums -> MKL VSL RNG GAMMA METHOD parameters
    std::map< MKLGammaMethodType, ParamType > method {
      { MKLGammaMethodType::GNORM, VSL_RNG_METHOD_GAMMA_GNORM },
      { MKLGammaMethodType::GNORM_ACCURATE,
        VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE }
    };
};

} // ctr::
} // tk::

#endif // MKLGammaMethodOptions_h
