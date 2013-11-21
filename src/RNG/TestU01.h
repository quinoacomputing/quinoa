//******************************************************************************
/*!
  \file      src/RNG/TestU01.h
  \author    J. Bakosi
  \date      Thu 21 Nov 2013 02:49:12 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 random number generator test suite
  \details   TestU01 random number generator test suite
*/
//******************************************************************************
#ifndef TestU01_h
#define TestU01_h

#include <Battery.h>

namespace rngtest {

//! TestU01 random number generator test suite
class TestU01 : public Battery {

  protected:
    //! Constructor
    explicit TestU01(const Base& base) : Battery(base) {}

    //! Destructor
    ~TestU01() noexcept override = default;

  private:
    //! Don't permit copy constructor
    TestU01(const TestU01&) = delete;
    //! Don't permit copy assigment
    TestU01& operator=(const TestU01&) = delete;
    //! Don't permit move constructor
    TestU01(TestU01&&) = delete;
    //! Don't permit move assigment
    TestU01& operator=(TestU01&&) = delete;
};

} // rngtest::

#endif // TestU01_h
