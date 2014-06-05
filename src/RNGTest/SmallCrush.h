//******************************************************************************
/*!
  \file      src/RNGTest/SmallCrush.h
  \author    J. Bakosi
  \date      Sat 24 May 2014 09:11:40 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

#include <TestU01SuitePolicy.h>

namespace rngtest {

//! SmallCrush : TestU01SuitePolicy
class SmallCrush : public TestU01SuitePolicy {

  public:
    //! Constructor
    explicit SmallCrush() = default;

    //! Return string identifying policy
    const std::string& policy() const;

    //! Add statistical tests to battery
    void addTests( std::vector< std::unique_ptr< StatTest > >& tests );

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
