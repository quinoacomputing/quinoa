//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.h
  \author    J. Bakosi
  \date      Sat 09 Nov 2013 02:45:10 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************
#ifndef SmallCrush_h
#define SmallCrush_h

#include <Battery.h>

namespace rngtest {

//! SmallCrush : Battery
class SmallCrush : public Battery {

  public:
    //! Constructor
    explicit SmallCrush(const Base& base) : Battery(base) {};

    //! Destructor
    virtual ~SmallCrush() noexcept = default;

    //! Run battery of RNG tests
    virtual void run() override;

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
