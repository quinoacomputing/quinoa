//******************************************************************************
/*!
  \file      src/MonteCarlo/TestSDE.h
  \author    J. Bakosi
  \date      Tue 31 Dec 2013 12:28:22 PM MST
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

    //! Initialize
    void init() override;

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

    //! One-liner report
    void reportHeader() const;
    void report( const uint64_t it,
                 const uint64_t nstep,
                 const tk::real t,
                 const tk::real dt,
                 const bool wroteJpdf,
                 const bool wroteGlob,
                 const bool wrotePlot );

    //! Advance
    void advance( tk::real dt );

    //! Output joint scalar PDF
    void outJpdf( tk::real t );
};

} // quinoa::

#endif // TestSDE_h
