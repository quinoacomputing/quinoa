//******************************************************************************
/*!
  \file      src/RNG/Battery.h
  \author    J. Bakosi
  \date      Sat 09 Nov 2013 03:28:37 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Battery base
  \details   Battery base
*/
//******************************************************************************
#ifndef Battery_h
#define Battery_h

#include <Base.h>

namespace rngtest {

//! Battery
class Battery {

  public:
    //! Constructor
    explicit Battery(const Base& base) : m_base(base) {}

    //! Destructor
    virtual ~Battery() noexcept = default;

    //! Run battery of RNG tests
    virtual void run() = 0;

  private:
    //! Don't permit copy constructor
    Battery(const Battery&) = delete;
    //! Don't permit copy assigment
    Battery& operator=(const Battery&) = delete;
    //! Don't permit move constructor
    Battery(Battery&&) = delete;
    //! Don't permit move assigment
    Battery& operator=(Battery&&) = delete;

    const Base& m_base;                   //!< Essentials
};

} // rngtest::

#endif // Battery_h
