//******************************************************************************
/*!
  \file      src/RNG/StatTest.h
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 08:17:27 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistical test base
  \details   Statistical test base
*/
//******************************************************************************
#ifndef StatTest_h
#define StatTest_h

namespace rngtest {

//! StatTest base
class StatTest {

  public:
    //! Constructor
    explicit StatTest() = default;

    //! Destructor
    virtual ~StatTest() noexcept = default;

    //! Container types for the p-values and names of statistical tests
    using Pvals = std::vector< double >;
    using Names = std::vector< std::string >;

    //! Run
    virtual Pvals run() = 0;

    //! Test name accessor
    virtual const Names::value_type&
    name( const Names::size_type& i ) const = 0;

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
