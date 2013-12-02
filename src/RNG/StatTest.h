//******************************************************************************
/*!
  \file      src/RNG/StatTest.h
  \author    J. Bakosi
  \date      Mon 02 Dec 2013 05:43:10 PM MST
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

    //! Number of results/test accessor
    virtual const Names::size_type& nresult() const = 0;

    //! RNG enum accessor
    virtual const quinoa::ctr::RNGType& rng() const = 0;

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
