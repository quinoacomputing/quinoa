//******************************************************************************
/*!
  \file      src/RNGTest/StatTest.h
  \author    J. Bakosi
  \date      Wed 14 May 2014 07:40:38 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical test base
  \details   Statistical test base
*/
//******************************************************************************
#ifndef StatTest_h
#define StatTest_h

#include <charm++.h>

namespace rngtest {

//! StatTest base
class StatTest : public PUP::able {

  public:
    //! Constructor
    explicit StatTest() = default;

    //! Destructor
    virtual ~StatTest() = default;

    //! Run
    virtual void run( std::size_t id ) = 0;

    //! Test name accessor
    virtual const std::string& name( std::size_t i ) const = 0;

    //! Number of results/test accessor
    virtual std::size_t nstat() const = 0;

    //! RNG enum accessor
    virtual const tk::ctr::RNGType& rng() const = 0;

    //! RNG id accessor
    virtual std::size_t id() const = 0;

    //! Query whether test is failed
    virtual bool fail( std::size_t p ) const = 0;

    //! p-value accessors
    virtual double pval( std::size_t p ) const = 0;
    virtual std::string pvalstr( std::size_t p ) const = 0;

    //! Return number of failed tests
    virtual std::size_t nfail() const = 0;

    PUPable_abstract( StatTest );

  private:
    //! Don't permit copy constructor
    StatTest(const StatTest&) = delete;
    //! Don't permit copy assigment
    StatTest& operator=(const StatTest&) = delete;
    //! Don't permit move constructor
    StatTest(StatTest&&) = delete;
    //! Don't permit move assigment
    StatTest& operator=(StatTest&&) = delete;
};

} // rngtest::

#endif // StatTest_h
