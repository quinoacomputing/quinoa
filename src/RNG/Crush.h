//******************************************************************************
/*!
  \file      src/RNG/Crush.h
  \author    J. Bakosi
  \date      Sat 02 Nov 2013 10:48:49 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Crush battery
  \details   Crush battery
*/
//******************************************************************************
#ifndef Crush_h
#define Crush_h

#include <Battery.h>

namespace rngtest {

//! Crush : Battery
class Crush : public Battery {

  public:
    //! Constructor
    explicit Crush() = default;

    //! Destructor
    virtual ~Crush() noexcept = default;

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
