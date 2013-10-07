//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Mon Oct  7 08:33:39 2013
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
class RNGTestDriver : public tk::Driver {

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

} // rngtest::

#endif // RNGTestDriver_h
