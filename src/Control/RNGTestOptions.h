//******************************************************************************
/*!
  \file      src/Control/RNGTestOptions.h
  \author    J. Bakosi
  \date      Thu Aug 29 16:41:43 2013
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

namespace rngtest {

namespace select {

//! Random number generator test types
enum class RNGTestType : uint8_t { NO_RNGTEST=0,
                                   SMALLCRUSH,
                                   CRUSH,
                                   BIGCRUSH };

using quinoa::select::Toggle;

//! Class with base templated on the above enum class with associations
class RNGTest : public Toggle<RNGTestType> {

  public:
    //! Constructor: pass associations references to base, which will handle
    //! class-user interactions
    explicit RNGTest() : Toggle<RNGTestType>(names, values) {}

  private:
    //! Don't permit copy constructor
    RNGTest(const RNGTest&) = delete;
    //! Don't permit copy assigment
    RNGTest& operator=(const RNGTest&) = delete;
    //! Don't permit move constructor
    RNGTest(RNGTest&&) = delete;
    //! Don't permit move assigment
    RNGTest& operator=(RNGTest&&) = delete;

    //! Enums -> names
    const std::map<RNGTestType, std::string> names {
      { RNGTestType::NO_RNGTEST, "No RNG test suite" },
      { RNGTestType::SMALLCRUSH, "SmallCrush" },
      { RNGTestType::CRUSH, "Crush" },
      { RNGTestType::BIGCRUSH, "BigCrush" }
    };

    //! keywords -> Enums
    const std::map<std::string, RNGTestType> values {
      { "no_rngtest", RNGTestType::NO_RNGTEST },
      { "smallcrush", RNGTestType::SMALLCRUSH },
      { "crush", RNGTestType::CRUSH },
      { "bigcrush", RNGTestType::BIGCRUSH }
    };
};

} // namespace select

} // namespace rngtest

#endif // RNGTestOptions_h
