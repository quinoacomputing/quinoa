//******************************************************************************
/*!
  \file      src/RNGTest/TestU01.h
  \author    J. Bakosi
  \date      Wed Apr 23 13:36:08 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 statistical tests
  \details   TestU01 statistical tests
*/
//******************************************************************************
#ifndef TestU01_h
#define TestU01_h

#include <StatTest.h>

namespace rngtest {

//! TestU01 : StatTest
template< class Result,                             //!< Results type
          Result* (*Creator)(void),                 //!< Results creator
          void (*Deleter)(Result *),                //!< Results deleter
          typename... Ts >                          //!< Extra runner args types
class TestU01 : public StatTest {

    using Xargs = std::tuple< Ts... >;

    //! Test runner function pointer type
    using RunFn = StatTest::Pvals (*)( unif01_Gen*, Result*, const Xargs& );

  public:
    //! Constructor
    explicit TestU01( std::size_t i,
                      const unif01_Gen* const gen,
                      const tk::ctr::RNGType& r,
                      Names&& names,
                      RunFn runner,
                      Ts&&... xargs ) :
      m_id( i ),
      m_gen( gen ),
      m_rng( r ),
      m_npval( names.size() ),
      m_names( std::move(names) ),
      m_runner( runner ),
      m_xargs( std::forward<Ts>(xargs)... ),
      m_pvals( names.size(), -1.0 ),
      m_res( ResultPtr( Creator() ) ) {}

    //! Destructor
    ~TestU01() noexcept override = default;

    //! Run test, awful that TestU01 does not guarantee the constness of gen
    void run() override {
      m_pvals = m_runner( const_cast<unif01_Gen*>(m_gen), m_res.get(), m_xargs );
    }

    //! Test name accessor
    const Names::value_type& name( std::size_t i ) const override {
      Assert( i < m_names.size(), tk::ExceptType::FATAL,
              "Indexing outside of container bounds in TestU01::names()" );
      return m_names[i];
    }

    //! Number of results/test (i.e., p-values) accessor
    std::size_t nstat() const override { return m_npval; }

    //! RNG enum accessor
    const tk::ctr::RNGType& rng() const override { return m_rng; }

    //! RNG id accessor
    std::size_t id() const override { return m_id; }

    //! Query whether test is failed
    bool fail( std::size_t p ) const override {
      if ((m_pvals[p] <= gofw_Suspectp) || (m_pvals[p] >= 1.0-gofw_Suspectp))
        return true;
      else
        return false;
    }

    //! Return number of failed tests
    std::size_t nfail() const override {
      std::size_t f = 0;
      for (std::size_t p = 0; p<m_npval; ++p) {
        if (fail(p)) ++f;
      }
      return f;
    }

    //! p-value accessor
    double pval( std::size_t p ) const override {
      return m_pvals[p];
    }

    //! Return humand-readable p-value (ala TestU01::bbattery.c::WritePval)
    std::string pvalstr( std::size_t p ) const override {
      std::stringstream ss;
      double val = m_pvals[p];
      if (val < gofw_Suspectp) {

        if ((val >= 0.01) && (val <= 0.99))
          ss << val;
        else if (val < gofw_Epsilonp)
          ss << "eps";
        else if (val < 0.01)
          ss << val;
        else if (val >= 1.0 - gofw_Epsilonp1)
          ss << "1 - eps1";
        else if (val < 1.0 - 1.0e-4)
          ss << val;
        else
          ss << 1.0 - val;

      } else if (val > 1.0 - gofw_Suspectp) {

        if (val >= 1.0 - gofw_Epsilonp1)
          ss << "1 - eps1";
        else if (val >= 1.0 - 1.0e-4)
          ss << "1 - " << 1.0 - val;
        else
          ss << val;

      }
      return ss.str();
    }

  private:
    //! Don't permit copy constructor
    TestU01(const TestU01&) = delete;
    //! Don't permit copy assigment
    TestU01& operator=(const TestU01&) = delete;
    //! Don't permit move constructor
    TestU01(TestU01&&) = delete;
    //! Don't permit move assigment
    TestU01& operator=(TestU01&&) = delete;

    const std::size_t m_id;                    //!< RNG id
    const unif01_Gen* const m_gen;             //!< Raw ptr to TestU01 generator
    const tk::ctr::RNGType m_rng;              //!< RNG selected
    const std::size_t m_npval;                 //!< Number of p-values produced
    const Names m_names;                       //!< Name(s) of tests
    const RunFn m_runner;                      //!< Test runner function
    const Xargs m_xargs;                       //!< Extra args for run()

    StatTest::Pvals m_pvals;                   //!< p-value(s)

    //! TestU01 results type with a custom deleter by TestU01
    using ResultPtr = TestU01Ptr< Result, Deleter >;
    //! TestU01 results struct (wrapped to std::unique_ptr)
    ResultPtr m_res;
};

} // rngtest::

#endif // TestU01_h
