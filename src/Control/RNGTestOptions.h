//******************************************************************************
/*!
  \file      src/Control/RNGTestOptions.h
  \author    J. Bakosi
  \date      Tue 30 Jul 2013 08:07:13 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Turbulence frequency model options and associations
  \details   Turbulence frequency model options and associations
*/
//******************************************************************************
#ifndef RNGTestOptions_h
#define RNGTestOptions_h

#include <map>

#include <Exception.h>
#include <Toggle.h>

namespace Quinoa {

namespace select {

//! Random number generator test types
enum class RNGTestTypes : uint8_t { NO_RNGTEST=0,
                                    SMALLCRUSH,
                                    CRUSH,
                                    BIGCRUSH };

//! Class with base templated on the above enum class with associations
class RNGTest : public Toggle<RNGTestTypes> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    RNGTest() : Toggle<RNGTestTypes>(names, values) {
      //! Enums -> names
      names[RNGTestTypes::NO_RNGTEST] = "No RNG test suite";
      names[RNGTestTypes::SMALLCRUSH] = "SmallCrush";
      names[RNGTestTypes::CRUSH] = "Crush";
      names[RNGTestTypes::BIGCRUSH] = "BigCrush";
      //! keywords -> Enums
      values["no_rngtest"] = RNGTestTypes::NO_RNGTEST;
      values["smallcrush"] = RNGTestTypes::SMALLCRUSH;
      values["crush"] = RNGTestTypes::CRUSH;
      values["bigcrush"] = RNGTestTypes::BIGCRUSH;
    }

  private:
    //! Don't permit copy constructor
    RNGTest(const RNGTest&) = delete;
    //! Don't permit copy assigment
    RNGTest& operator=(const RNGTest&) = delete;
    //! Don't permit move constructor
    RNGTest(RNGTest&&) = delete;
    //! Don't permit move assigment
    RNGTest& operator=(RNGTest&&) = delete;

    std::map<RNGTestTypes, std::string> names;
    std::map<std::string, RNGTestTypes> values;
};

} // namespace select

} // namespace Quinoa

#endif // RNGTestOptions_h
