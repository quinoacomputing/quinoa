//******************************************************************************
/*!
  \file      src/Control/Options/MKLUniformMethod.h
  \author    J. Bakosi
  \date      Tue 05 Aug 2014 03:43:38 PM MDT
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
    const ParamType& param( MKLUniformMethodType rng ) const;

  private:
    //! Don't permit copy constructor
    MKLUniformMethod(const MKLUniformMethod&) = delete;
    //! Don't permit copy assigment
    MKLUniformMethod& operator=(const MKLUniformMethod&) = delete;
    //! Don't permit move constructor
    MKLUniformMethod(MKLUniformMethod&&) = delete;
    //! Don't permit move assigment
    MKLUniformMethod& operator=(MKLUniformMethod&&) = delete;

    //! Enums -> MKL VSL RNG UNIFORM METHOD parameters
    const std::map< MKLUniformMethodType, ParamType > method {
      { MKLUniformMethodType::STANDARD, VSL_RNG_METHOD_UNIFORM_STD },
      { MKLUniformMethodType::ACCURATE, VSL_RNG_METHOD_UNIFORM_STD_ACCURATE }
    };
};

} // ctr::
} // tk::

#endif // MKLUniformMethodOptions_h
