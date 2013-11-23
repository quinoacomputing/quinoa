//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.h
  \author    J. Bakosi
  \date      Sat 23 Nov 2013 12:25:08 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

extern "C" {
  #include <unif01.h>
}

#include <TestU01Wrap.h>
#include <TestU01.h>

namespace rngtest {

//! SmallCrush : TestU01
class SmallCrush : public TestU01 {

  public:
    //! Constructor
    explicit SmallCrush(const Base& base);

    //! Destructor
    ~SmallCrush() noexcept override = default;

    //! Run battery of RNG tests
    void run() override;

  private:
    //! Don't permit copy constructor
    SmallCrush(const SmallCrush&) = delete;
    //! Don't permit copy assigment
    SmallCrush& operator=(const SmallCrush&) = delete;
    //! Don't permit move constructor
    SmallCrush(SmallCrush&&) = delete;
    //! Don't permit move assigment
    SmallCrush& operator=(SmallCrush&&) = delete;

    //! Name of RNGs
    std::string m_rngname;

    //! Statistical tests
    TestContainer m_tests;
};

} // rngtest::

#endif // SmallCrush_h
