//******************************************************************************
/*!
  \file      src/RNGTest/BigCrush.h
  \author    J. Bakosi
  \date      Wed 09 Apr 2014 04:08:13 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************
#ifndef BigCrush_h
#define BigCrush_h

#include <TestU01Util.h>
#include <TestU01Suite.h>

namespace rngtest {

//! BigCrush : TestU01Suite
class BigCrush : public TestU01Suite {

  public:
    //! Constructor
    explicit BigCrush( const Base& base );

    //! Destructor
    ~BigCrush() noexcept override = default;

    //! Add statistical tests to battery
    void addTests( const StatTest::Rsize& id,
                   const tk::ctr::RNGType& rng,
                   const Gen01Ptr& gen );

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
