//******************************************************************************
/*!
  \file      src/RNG/TestU01Suite.h
  \author    J. Bakosi
  \date      Mon 25 Nov 2013 11:12:42 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 random number generator test suite
  \details   TestU01 random number generator test suite
*/
//******************************************************************************
#ifndef TestU01Suite_h
#define TestU01Suite_h

extern "C" {
  #include <unif01.h>
  #include <sres.h>
}

#include <Battery.h>
#include <TestU01Wrap.h>

namespace rngtest {

//! TestU01 random number generator test suite
class TestU01Suite : public Battery {

  protected:
    //! Constructor
    explicit TestU01Suite(const Base& base) : Battery(base) {}

    //! Destructor
    ~TestU01Suite() noexcept override = default;

    //! TestU01 external generator type with a custom deleter by TestU01
    using Gen01Ptr = TestU01Ptr< unif01_Gen, unif01_DeleteExternGen01 >;
    //! TestU01 external generator
    Gen01Ptr m_gen;

    //! Marsaglia's BirthdaySpacing test
    static double BirthdaySpacings( unif01_Gen* gen, sres_Poisson* res );

  private:
    //! Don't permit copy constructor
    TestU01Suite(const TestU01Suite&) = delete;
    //! Don't permit copy assigment
    TestU01Suite& operator=(const TestU01Suite&) = delete;
    //! Don't permit move constructor
    TestU01Suite(TestU01Suite&&) = delete;
    //! Don't permit move assigment
    TestU01Suite& operator=(TestU01Suite&&) = delete;

    static const long THOUSAND = 1000;
    static const long MILLION = THOUSAND * THOUSAND;
    static const long BILLION = THOUSAND * MILLION;
};

} // rngtest::

#endif // TestU01Suite_h
