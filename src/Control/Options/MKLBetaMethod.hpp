// *****************************************************************************
/*!
  \file      src/Control/Options/MKLBetaMethod.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Intel MKL Beta RNG method options
  \details   Intel MKL Beta RNG method options
*/
// *****************************************************************************
#ifndef MKLBetaMethodOptions_h
#define MKLBetaMethodOptions_h

#include <map>

#include <brigand/sequences/list.hpp>

#include <mkl_vsl_defines.h>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! MKL beta random number generator method types
enum class MKLBetaMethodType : uint8_t { CJA,
                                         CJA_ACCURATE };

//! \brief Pack/Unpack MKLBetaMethodType: forward overload to generic enum
//!   class packer
inline void operator|( PUP::er& p, MKLBetaMethodType& e ) { PUP::pup( p, e ); }

//! \brief MKLBetaMethod options: outsource searches to base templated on
//!   enum type
class MKLBetaMethod : public tk::Toggle< MKLBetaMethodType > {

  public:
    using ParamType = int;

    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::cja
                                  , kw::cja_accurate
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit MKLBetaMethod() :
      tk::Toggle< MKLBetaMethodType >(
        //! Group, i.e., options, name
        "Beta method",
        //! Enums -> names
        { { MKLBetaMethodType::CJA, kw::cja::name() },
          { MKLBetaMethodType::CJA_ACCURATE, kw::cja_accurate::name() } },
        //! keywords -> Enums
        { { kw::cja::string(), MKLBetaMethodType::CJA },
          { kw::cja_accurate::string(), MKLBetaMethodType::CJA_ACCURATE } } ) {}

    //! \brief Return parameter based on Enum
    //! \details Here 'parameter' is the library-specific identifier of the
    //!    option, i.e., as the library identifies the given option
    //! \param[in] m Enum value of the option requested
    //! \return Library-specific parameter of the option
    const ParamType& param( MKLBetaMethodType m ) const {
      using tk::operator<<;
      auto it = method.find( m );
      Assert( it != end(method),
              std::string("Cannot find parameter for MKLBetaMethod \"")
              << m << "\"" );
      return it->second;
    }

  private:
    //! Enums -> MKL VSL RNG BETA METHOD parameters
    std::map< MKLBetaMethodType, ParamType > method {
      { MKLBetaMethodType::CJA, VSL_RNG_METHOD_BETA_CJA },
      { MKLBetaMethodType::CJA_ACCURATE, VSL_RNG_METHOD_BETA_CJA_ACCURATE }
    };
};

} // ctr::
} // tk::

#endif // MKLBetaMethodOptions_h
