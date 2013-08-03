//******************************************************************************
/*!
  \file      src/Control/RNGTestOptions.h
  \author    J. Bakosi
  \date      Fri Aug  2 15:43:06 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite options and associations
  \details   Random number generator test suite options and associations
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
enum class RNGTestType : uint8_t { NO_RNGTEST=0,
                                   SMALLCRUSH,
                                   CRUSH,
                                   BIGCRUSH };

//! Class with base templated on the above enum class with associations
class RNGTest : public Toggle<RNGTestType> {

  public:
    //! Constructor initializing associations
    // ICC: use initializer lists
    RNGTest() : Toggle<RNGTestType>(names, values) {
      //! Enums -> names
      names[RNGTestType::NO_RNGTEST] = "No RNG test suite";
      names[RNGTestType::SMALLCRUSH] = "SmallCrush";
      names[RNGTestType::CRUSH] = "Crush";
      names[RNGTestType::BIGCRUSH] = "BigCrush";
      //! keywords -> Enums
      values["no_rngtest"] = RNGTestType::NO_RNGTEST;
      values["smallcrush"] = RNGTestType::SMALLCRUSH;
      values["crush"] = RNGTestType::CRUSH;
      values["bigcrush"] = RNGTestType::BIGCRUSH;
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

    std::map<RNGTestType, std::string> names;
    std::map<std::string, RNGTestType> values;
};

} // namespace select

} // namespace Quinoa

#endif // RNGTestOptions_h
