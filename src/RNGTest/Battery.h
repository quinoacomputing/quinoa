//******************************************************************************
/*!
  \file      src/RNGTest/Battery.h
  \author    J. Bakosi
  \date      Sat 24 May 2014 01:28:08 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Battery base
  \details   Battery base
*/
//******************************************************************************
#ifndef Battery_h
#define Battery_h

#include <charm++.h>

#include <StatTest.h>
#include <RNGTestPrint.h>

namespace rngtest {

//! Battery
class Battery : public PUP::able {

  public:
    //! Constructor
    explicit Battery() = default;

    //! Destructor
    virtual ~Battery() = default;

    //! Run battery of RNG tests
    virtual void run() = 0;

    //! Print list of registered statistical tests
    virtual void print( const RNGTestPrint& print ) const = 0;

    //! Return number of statistical tests in battery
    virtual std::size_t ntest() const = 0;

    //! Return number of statistics produced by battery
    virtual std::size_t nstat() const = 0;

    //! Evaluate a statistical test
    virtual void evaluate( std::size_t id ) = 0;

    //! Enable PUP for virtual function pointers,
    //! see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html,
    //! section "Subclass allocation via PUP::able"
    PUPable_abstract( Battery );  

  private:
    //! Don't permit copy constructor
    Battery(const Battery&) = delete;
    //! Don't permit copy assigment
    Battery& operator=(const Battery&) = delete;
    //! Don't permit move constructor
    Battery(Battery&&) = delete;
    //! Don't permit move assigment
    Battery& operator=(Battery&&) = delete;
};

//! Runner chare for statistical test batteries (a full suite of statistical
//! tests) deriving from Battery. This wrapper is used for virtual dispatch of
//! batteries with base Battery. See also section "Subclass allocation via
//! PUP::able" in http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
class BatteryRunner : public Chare {
  public:
    BatteryRunner( CkMigrateMessage* ) {}
    BatteryRunner( Battery& battery ) {
      std::cout << "BatteryRunner::BatteryRunner()" << std::endl;
      //battery.run();
      delete this;
    }
};

} // rngtest::

#endif // Battery_h
