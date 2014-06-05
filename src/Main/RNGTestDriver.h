//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Mon 26 May 2014 04:52:53 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Driver random number test suite driver
  \details   Driver random number test suite driver
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <RNGDriver.h>
#include <Base.h>
#include <Battery.h>

//! Everything that contributes to the rngtest executable
namespace rngtest {

//! Battery factory type
using BatteryFactory = std::map< ctr::BatteryType, std::function< Battery*() > >;

//! RNGTestDriver
class RNGTestDriver : public tk::RNGDriver {

  public:
    //! Constructor
    explicit RNGTestDriver( int argc, char** argv, const tk::Print& print );

    //! Execute driver
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

//     //! Create all global-scope RNG wrappers based on user selection
//     void initRNGs();
// 
//     //! Register a single global-scope RNG wrapper
//     template< std::size_t gid >
//     void registerRNG( tk::ctr::RNGType r ) {
//       auto it = m_base->rng.find( r );
//       if ( it != m_base->rng.end() ) {
//         g_rng[gid] = std::unique_ptr< tk::RNG >( it->second() );
//       } else {
//         Throw( tk::ExceptType::FATAL, "RNG not found in factory" );
//       }
//     }

    //! Echo information on random number test suite
    void echo();

    //! Number of test to run
    std::size_t m_ntest;

    //! Factories
    tk::RNGFactory m_RNGFactory;                 //!< RNG factory
    BatteryFactory m_batteryFactory;             //!< Battery factory

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
