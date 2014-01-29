//******************************************************************************
/*!
  \file      src/MonteCarlo/TestSDE.h
  \author    J. Bakosi
  \date      Tue 28 Jan 2014 05:04:12 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SDE testbed
  \details   SDE testbed
*/
//******************************************************************************
#ifndef TestSDE_h
#define TestSDE_h

#include <MonteCarlo.h>
#include <Timer.h>

namespace quinoa {

//! TestSDE : MonteCarlo
class TestSDE : public MonteCarlo {

  public:
    //! Constructor
    explicit TestSDE( const Base& base ) : MonteCarlo( base ) {}

    //! Destructor
    ~TestSDE() override = default;

    //! Run
    void run() override;

  private:
    //! Don't permit copy constructor
    TestSDE(const TestSDE&) = delete;
    //! Don't permit copy assigment
    TestSDE& operator=(const TestSDE&) = delete;
    //! Don't permit move constructor
    TestSDE(TestSDE&&) = delete;
    //! Don't permit move assigment
    TestSDE& operator=(TestSDE&&) = delete;

    //! Advance
    void advance( tk::real dt );

    //! Output joint scalar PDF
    void outJpdf( tk::real t );
};

} // quinoa::

#endif // TestSDE_h
