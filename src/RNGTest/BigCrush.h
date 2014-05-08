//******************************************************************************
/*!
  \file      src/RNGTest/BigCrush.h
  \author    J. Bakosi
  \date      Wed 07 May 2014 10:43:24 PM MDT
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
    ~BigCrush() override = default;

    //! Add statistical tests to battery
    void addTests( std::size_t id,
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
