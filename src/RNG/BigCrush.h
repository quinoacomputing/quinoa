//******************************************************************************
/*!
  \file      src/RNG/BigCrush.h
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 06:24:42 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************
#ifndef BigCrush_h
#define BigCrush_h

#include <TestU01Suite.h>

namespace rngtest {

//! BigCrush : TestU01Suite
class BigCrush : public TestU01Suite {

  public:
    //! Constructor
    explicit BigCrush(const Base& base) : TestU01Suite(base) {}

    //! Destructor
    ~BigCrush() noexcept override = default;

    //! Run battery of RNG tests
    void run() override { std::cout << "BigCrush::run() unimplemented!" << std::endl; }

    //! Print list of registered statistical tests
    void print() const override {  std::cout << "BigCrush::print() unimplemented!" << std::endl; }

  private:
    //! Don't permit copy constructor
    BigCrush(const BigCrush&) = delete;
    //! Don't permit copy assigment
    BigCrush& operator=(const BigCrush&) = delete;
    //! Don't permit move constructor
    BigCrush(BigCrush&&) = delete;
    //! Don't permit move assigment
    BigCrush& operator=(BigCrush&&) = delete;
};

} // rngtest::

#endif // BigCrush_h
