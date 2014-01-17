//******************************************************************************
/*!
  \file      src/Control/Options/MKLUniformMethod.h
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 08:36:28 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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

namespace quinoa {
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
      Toggle< MKLUniformMethodType >("uniform method", names, values) {}

    //! Return parameter based on Enum
    const ParamType& param(MKLUniformMethodType rng) const;

  private:
    //! Don't permit copy constructor
    MKLUniformMethod(const MKLUniformMethod&) = delete;
    //! Don't permit copy assigment
    MKLUniformMethod& operator=(const MKLUniformMethod&) = delete;
    //! Don't permit move constructor
    MKLUniformMethod(MKLUniformMethod&&) = delete;
    //! Don't permit move assigment
    MKLUniformMethod& operator=(MKLUniformMethod&&) = delete;

    //! Get access to MKL uniform method keywords
    const tk::kw::standard standard {};
    const tk::kw::accurate accurate {};

    //! Enums -> names
    const std::map< MKLUniformMethodType, std::string > names {
      { MKLUniformMethodType::STANDARD, standard.name() },
      { MKLUniformMethodType::ACCURATE, accurate.name() }
    };

    //! keywords -> Enums
    const std::map< std::string, MKLUniformMethodType > values {
      { standard.string(), MKLUniformMethodType::STANDARD },
      { accurate.string(), MKLUniformMethodType::ACCURATE }
    };

    //! Enums -> MKL VSL RNG UNIFORM METHOD parameters
    const std::map< MKLUniformMethodType, ParamType > method {
      { MKLUniformMethodType::STANDARD, VSL_RNG_METHOD_UNIFORM_STD },
      { MKLUniformMethodType::ACCURATE, VSL_RNG_METHOD_UNIFORM_STD_ACCURATE }
    };
};

} // ctr::
} // quinoa::

#endif // MKLUniformMethodOptions_h
