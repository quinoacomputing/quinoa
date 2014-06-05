//******************************************************************************
/*!
  \file      src/RNGTest/BigCrush.h
  \author    J. Bakosi
  \date      Sat 24 May 2014 09:11:24 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************
#ifndef BigCrush_h
#define BigCrush_h

#include <TestU01SuitePolicy.h>

namespace rngtest {

//! BigCrush : TestU01Suite
class BigCrush : public TestU01SuitePolicy {

  public:
    //! Constructor
    explicit BigCrush() = default;

    //! Return string identifying policy
    const std::string& policy() const;

    //! Add statistical tests to battery
    void addTests( std::vector< std::unique_ptr< StatTest > >& tests );

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
