//******************************************************************************
/*!
  \file      src/Control/Options/MKLUniformMethod.h
  \author    J. Bakosi
  \date      Mon 15 Sep 2014 07:59:18 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Intel MKL uniform RNG method options
  \details   Intel MKL uniform RNG method options
*/
//******************************************************************************
#ifndef MKLUniformMethodOptions_h
#define MKLUniformMethodOptions_h

#include <map>

#include <mkl_vsl.h>

#include <Toggle.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace tk {
namespace ctr {

//! MKL uniform random number generator methods
enum class MKLUniformMethodType : uint8_t { STANDARD,
                                            ACCURATE };

//! Pack/Unpack: delegate to tk::
inline void operator|( PUP::er& p, MKLUniformMethodType& e )
{ tk::pup( p, e ); }

//! Class with base templated on the above enum class with associations
class MKLUniformMethod : public tk::Toggle< MKLUniformMethodType > {

  public:
    using ParamType = int;

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MKLUniformMethod() :
      Toggle< MKLUniformMethodType >( "uniform method",
        //! Enums -> names
        { { MKLUniformMethodType::STANDARD, tk::kw::standard().name() },
          { MKLUniformMethodType::ACCURATE, tk::kw::accurate().name() } },
        //! keywords -> Enums
        { { tk::kw::standard().string(), MKLUniformMethodType::STANDARD },
          { tk::kw::accurate().string(), MKLUniformMethodType::ACCURATE } } ) {}

    //! Return parameter based on Enum
    const ParamType& param( MKLUniformMethodType m ) const {
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
