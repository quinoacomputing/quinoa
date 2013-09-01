//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Sun 01 Sep 2013 02:01:31 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver random number test suite driver
  \details   Driver random number test suite driver
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <Driver.h>

namespace rngtest {

void MKLErrChk(int vslerr);

//! RNGTestDriver base class
class RNGTestDriver : public quinoa::Driver {

  public:
    //! Constructor
    explicit RNGTestDriver(int argc, char** argv);

    //! Destructor
    ~RNGTestDriver() noexcept override;

    //! Finalize, single exit point, called implicitly from destructor or
    //! explicitly from anywhere else
    void finalize() noexcept override;

    //! Solve
    void execute() const override;

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

} // namespace rngtest

#endif // RNGTestDriver_h
