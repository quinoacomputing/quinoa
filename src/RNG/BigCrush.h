//******************************************************************************
/*!
  \file      src/RNG/BigCrush.h
  \author    J. Bakosi
  \date      Sat 02 Nov 2013 10:49:06 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************
#ifndef BigCrush_h
#define BigCrush_h

#include <Battery.h>

namespace rngtest {

//! BigCrush : Battery
class BigCrush : public Battery {

  public:
    //! Constructor
    explicit BigCrush() = default;

    //! Destructor
    virtual ~BigCrush() noexcept = default;

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
