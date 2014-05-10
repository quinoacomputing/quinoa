//******************************************************************************
/*!
  \file      src/RNGTest/SmallCrush.h
  \author    J. Bakosi
  \date      Sat 10 May 2014 10:19:10 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

#include <TestU01Util.h>
#include <TestU01Suite.h>

namespace rngtest {

//! SmallCrush : TestU01Suite
class SmallCrush : public TestU01Suite {

  public:
    //! Constructor
    explicit SmallCrush(const Base& base);

    //! Destructor
    ~SmallCrush() override = default;

    //! Add statistical tests to battery
    void addTests( std::size_t id,
                   tk::ctr::RNGType& rng,
                   Gen01Ptr& gen );

  private:
    //! Don't permit copy constructor
    SmallCrush(const SmallCrush&) = delete;
    //! Don't permit copy assigment
    SmallCrush& operator=(const SmallCrush&) = delete;
    //! Don't permit move constructor
    SmallCrush(SmallCrush&&) = delete;
    //! Don't permit move assigment
    SmallCrush& operator=(SmallCrush&&) = delete;
};

} // rngtest::

#endif // SmallCrush_h
