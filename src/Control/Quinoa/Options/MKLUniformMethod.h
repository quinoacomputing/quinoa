//******************************************************************************
/*!
  \file      src/Control/Quinoa/Options/MKLUniformMethod.h
  \author    J. Bakosi
  \date      Wed 06 Nov 2013 10:21:27 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Intel MKL uniform RNG method options
  \details   Intel MKL uniform RNG method options
*/
//******************************************************************************
#ifndef QuinoaMKLUniformMethodOptions_h
#define QuinoaMKLUniformMethodOptions_h

#include <map>

#include <Toggle.h>

namespace quinoa {
namespace ctr {

//! MKL uniform random number generator methods
enum class MKLUniformMethodType : uint8_t { NO_METHOD=0,
                                            STANDARD,
                                            ACCURATE };

//! Class with base templated on the above enum class with associations
class MKLUniformMethod : public tk::Toggle< MKLUniformMethodType > {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit MKLUniformMethod() :
      Toggle< MKLUniformMethodType >("MKL uniform method", names, values) {}

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
    const kw::standard standard {};
    const kw::accurate accurate {};

    //! Enums -> names
    const std::map< MKLUniformMethodType, std::string > names {
      { MKLUniformMethodType::NO_METHOD, "n/a" },
      { MKLUniformMethodType::STANDARD, standard.name() },
      { MKLUniformMethodType::ACCURATE, accurate.name() }
    };

    //! keywords -> Enums
    const std::map< std::string, MKLUniformMethodType > values {
      { "no_method", MKLUniformMethodType::NO_METHOD },
      { standard.string(), MKLUniformMethodType::STANDARD },
      { accurate.string(), MKLUniformMethodType::ACCURATE }
    };
};

} // ctr::
} // quinoa::

#endif // QuinoaMKLUniformMethodOptions_h
