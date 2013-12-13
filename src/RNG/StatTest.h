//******************************************************************************
/*!
  \file      src/RNG/StatTest.h
  \author    J. Bakosi
  \date      Thu 12 Dec 2013 08:17:01 PM MST
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

    using Psize = Pvals::size_type;
    using Nsize = Names::size_type;

    //! Run
    virtual void run() = 0;

    //! Test name accessor
    virtual const Names::value_type& name( const Nsize& i ) const = 0;

    //! Number of results/test accessor
    virtual const Nsize& nresult() const = 0;

    //! RNG enum accessor
    virtual const quinoa::ctr::RNGType& rng() const = 0;

    //! Query whether test is failed
    virtual bool fail( const Nsize& p ) const = 0;

    //! p-value accessors
    virtual double pval( const Nsize& p ) const = 0;
    virtual std::string pvalstr( const Nsize& p ) const = 0;

    //! Return number of failed tests
    virtual Psize nfail() const = 0;

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
