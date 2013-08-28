//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Wed Aug 28 15:11:03 2013
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

//! RNGTestDriver base class
class RNGTestDriver : public Driver {

  public:
    //! Constructor
    explicit RNGTestDriver(int argc, char** argv);

    //! Destructor
    virtual ~RNGTestDriver() noexcept;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    virtual void finalize() noexcept;

    //! Solve
    virtual void execute() const;

  private:
    //! Don't permit copy constructor
    RNGTestDriver(const RNGTestDriver&) = delete;
    //! Don't permit assigment constructor
    RNGTestDriver& operator=(const RNGTestDriver&) = delete;
    //! Don't permit move constructor
    RNGTestDriver(RNGTestDriver&&) = delete;
    //! Don't permit move assignment
    RNGTestDriver& operator=(RNGTestDriver&&) = delete;
};

} // namespace Quinoa

#endif // RNGTestDriver_h
