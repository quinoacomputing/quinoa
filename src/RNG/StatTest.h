//******************************************************************************
/*!
  \file      src/RNG/StatTest.h
  \author    J. Bakosi
  \date      Sat 23 Nov 2013 04:44:29 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical test base
  \details   Statistical test base
*/
//******************************************************************************
#ifndef StatTest_h
#define StatTest_h

namespace rngtest {

//! StatTest base
class StatTest {

  public:
    //! Constructor
    explicit StatTest() = default;

    //! Destructor
    virtual ~StatTest() noexcept = default;

    //! Run
    virtual double run() = 0;

    //! Test name accessor
    virtual const char* name() const = 0;

  protected:
    static const long THOUSAND = 1000;
    static const long MILLION = THOUSAND * THOUSAND;
    static const long BILLION = THOUSAND * MILLION;

  private:
    //! Don't permit copy constructor
    StatTest(const StatTest&) = delete;
    //! Don't permit copy assigment
    StatTest& operator=(const StatTest&) = delete;
    //! Don't permit move constructor
    StatTest(StatTest&&) = delete;
    //! Don't permit move assigment
    StatTest& operator=(StatTest&&) = delete;
};

} // rngtest::

#endif // StatTest_h
