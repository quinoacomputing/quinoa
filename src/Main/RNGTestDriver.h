//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Sat 02 Nov 2013 12:30:38 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver random number test suite driver
  \details   Driver random number test suite driver
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <Driver.h>
#include <Base.h>
#include <Battery.h>
#include <RNG.h>

//! Everything that contributes to the rngtest executable
namespace rngtest {

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

    //! Initialize factories
    void initFactories(const tk::Print& print);

    //! Echo information on random number test suite
    void echo();

    //! Factories
    quinoa::ctr::RNGFactory m_RNGFactory;        //!< RNG factory
    ctr::BatteryFactory m_batteryFactory;        //!< Battery factory

    //! Pointers to selected options
    std::unique_ptr< tk::RNG > m_rng;            //!< Random number generator
    std::unique_ptr< Battery > m_battery;        //!< Battery

    //! Pointers to essentials to be created
    std::unique_ptr< ctr::InputDeck > m_control; //!< Control
    std::unique_ptr< RNGTestPrint > m_print;     //!< Pretty printer
    std::unique_ptr< tk::Paradigm > m_paradigm;  //!< Parallel compute env.
    std::unique_ptr< tk::Timer > m_timer;        //!< Timer
    std::unique_ptr< Base > m_base;              //!< Essentials bundle
};

} // rngtest::

#endif // RNGTestDriver_h
