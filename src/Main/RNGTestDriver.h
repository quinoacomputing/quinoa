//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 03:31:37 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver random number test suite driver
  \details   Driver random number test suite driver
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <Driver.h>
#include <Base.h>

//! Everything that contributes to the rngtest executable
namespace rngtest {

void MKLErrChk(int vslerr);

//! RNGTestDriver base class
class RNGTestDriver : public quinoa::Driver {

  public:
    //! Constructor
    explicit RNGTestDriver(int argc, char** argv, Base& base);

    //! Destructor
    ~RNGTestDriver() noexcept override = default;

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

    Base& m_base;                           //!< Essentials
};

} // namespace rngtest

#endif // RNGTestDriver_h
