// *****************************************************************************
/*!
  \file      src/Control/Options/MKLUniformMethod.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Intel MKL uniform RNG method options
  \details   Intel MKL uniform RNG method options
*/
// *****************************************************************************
#ifndef MKLUniformMethodOptions_h
#define MKLUniformMethodOptions_h

#include <map>

#include <brigand/sequences/list.hpp>

#include <mkl_vsl_defines.h>

#include "Toggle.hpp"
#include "Keywords.hpp"
#include "PUPUtil.hpp"

namespace tk {
namespace ctr {

//! MKL uniform random number generator method types
enum class MKLUniformMethodType : uint8_t { STANDARD,
                                            ACCURATE };

//! \brief Pack/Unpack MKLGUniformMethodType: forward overload to generic enum
//!   class packer
inline void operator|( PUP::er& p, MKLUniformMethodType& e )
{ PUP::pup( p, e ); }

//! \brief MKLUniformMethod options: outsource searches to base templated on
//!   enum type
class MKLUniformMethod : public tk::Toggle< MKLUniformMethodType > {

  public:
    using ParamType = int;

    //! Valid expected choices to make them also available at compile-time
    using keywords = brigand::list< kw::standard
                                  , kw::accurate
                                  >;

    //! \brief Options constructor
    //! \details Simply initialize in-line and pass associations to base, which
    //!    will handle client interactions
    explicit MKLUniformMethod() :
      tk::Toggle< MKLUniformMethodType >(
        //! Group, i.e., options, name
        "uniform method",
        //! Enums -> names
        { { MKLUniformMethodType::STANDARD, kw::standard::name() },
          { MKLUniformMethodType::ACCURATE, kw::accurate::name() }
        },
        //! keywords -> Enums
        { { kw::standard::string(), MKLUniformMethodType::STANDARD },
          { kw::accurate::string(), MKLUniformMethodType::ACCURATE } } ) {}

    //! \brief Return parameter based on Enum
    //! \details Here 'parameter' is the library-specific identifier of the
    //!    option, i.e., as the library identifies the given option
    //! \param[in] m Enum value of the option requested
    //! \return Library-specific parameter of the option
    const ParamType& param( MKLUniformMethodType m ) const {
      using tk::operator<<;
      auto it = method.find( m );
      Assert( it != end(method),
              std::string("Cannot find parameter for MKLUniformMethod \"")
              << m << "\"" );
      return it->second;
    }

  private:
    //! Enums -> MKL VSL RNG UNIFORM METHOD parameters
    std::map< MKLUniformMethodType, ParamType > method {
      { MKLUniformMethodType::STANDARD, VSL_RNG_METHOD_UNIFORM_STD },
      { MKLUniformMethodType::ACCURATE, VSL_RNG_METHOD_UNIFORM_STD_ACCURATE }
    };
};

} // ctr::
} // tk::

#endif // MKLUniformMethodOptions_h
