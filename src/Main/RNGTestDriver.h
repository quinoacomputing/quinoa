//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Fri Oct 18 12:21:25 2013
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
    explicit RNGTestDriver(int argc, char** argv, const tk::Print& print);

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

    std::unique_ptr< ctr::InputDeck > m_control; //!< Control
    std::unique_ptr< RNGTestPrint > m_print;     //!< Pretty printer
    std::unique_ptr< tk::Paradigm > m_paradigm;  //!< Parallel compute env.
    std::unique_ptr< tk::Timer > m_timer;        //!< Timer
    std::unique_ptr< Base > m_base;              //!< Essentials bundle
};

} // rngtest::

#endif // RNGTestDriver_h
