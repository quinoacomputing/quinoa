//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Wed Jul 31 13:45:09 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver random number test suite driver
  \details   Driver random number test suite driver
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <Driver.h>

namespace Quinoa {

void MKLErrChk(int vslerr);

class Timer;

//! RNGTestDriver base class
class RNGTestDriver : public Driver {

  public:
    //! Constructor
    RNGTestDriver(int argc, char** argv);

    //! Destructor
    virtual ~RNGTestDriver() noexcept;

    //! Solve
    virtual void execute() const;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept;

  private:
    //! Don't permit copy constructor
    RNGTestDriver(const RNGTestDriver&) = delete;
    //! Don't permit assigment constructor
    RNGTestDriver& operator=(const RNGTestDriver&) = delete;
    //! Don't permit move constructor
    RNGTestDriver(RNGTestDriver&&) = delete;
    //! Don't permit move assignment
    RNGTestDriver& operator=(RNGTestDriver&&) = delete;

    Timer* m_timer;                   //!< Timer object
};

} // namespace Quinoa

#endif // RNGTestDriver_h
