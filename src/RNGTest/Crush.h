//******************************************************************************
/*!
  \file      src/RNGTest/Crush.h
  \author    J. Bakosi
  \date      Sat 24 May 2014 09:12:00 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Crush battery
  \details   Crush battery
*/
//******************************************************************************
#ifndef Crush_h
#define Crush_h

#include <TestU01SuitePolicy.h>

namespace rngtest {

//! Crush : TestU01Suite
class Crush : public TestU01SuitePolicy {

  public:
    //! Constructor
    explicit Crush() = default;

    //! Return string identifying policy
    const std::string& policy() const;

    //! Add statistical tests to battery
    void addTests( std::vector< std::unique_ptr< StatTest > >& tests );

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
