//******************************************************************************
/*!
  \file      src/RNG/Crush.h
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 06:25:04 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Crush battery
  \details   Crush battery
*/
//******************************************************************************
#ifndef Crush_h
#define Crush_h

#include <TestU01Suite.h>

namespace rngtest {

//! Crush : TestU01Suite
class Crush : public TestU01Suite {

  public:
    //! Constructor
    explicit Crush(const Base& base) : TestU01Suite(base) {}

    //! Destructor
    ~Crush() noexcept override = default;

    //! Run battery of RNG tests
    void run() override { std::cout << "Crush::run() unimplemented!" << std::endl; }

    //! Print list of registered statistical tests
    void print() const override {  std::cout << "Crush::print() unimplemented!" << std::endl; }

  private:
    //! Don't permit copy constructor
    Crush(const Crush&) = delete;
    //! Don't permit copy assigment
    Crush& operator=(const Crush&) = delete;
    //! Don't permit move constructor
    Crush(Crush&&) = delete;
    //! Don't permit move assigment
    Crush& operator=(Crush&&) = delete;
};

} // rngtest::

#endif // Crush_h
