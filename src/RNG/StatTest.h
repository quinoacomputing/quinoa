//******************************************************************************
/*!
  \file      src/RNG/StatTest.h
  \author    J. Bakosi
  \date      Fri 22 Nov 2013 05:36:53 PM MST
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
