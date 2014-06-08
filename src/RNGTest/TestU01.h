//******************************************************************************
/*!
  \file      src/RNGTest/TestU01.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 01:37:15 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 statistical tests
  \details   TestU01 statistical tests
*/
//******************************************************************************
#ifndef TestU01_h
#define TestU01_h

extern "C" {
  #include <sres.h>
  #include <sknuth.h>
}

#include <TestU01Util.h>

namespace rngtest {

//! TestU01 statistical test used polymorphically with StatTest
template< class Result,                 //!< Results type
          Result* (*Creator)(void),     //!< Results creator
          void (*Deleter)(Result *),    //!< Results deleter
          typename... Ts >              //!< Extra runner args types
class TestU01 {

    //! Test extra arguments type
    using Xargs = std::tuple< Ts... >;

    //! Test runner function pointer type
    using RunFn = std::vector<double> (*)( unif01_Gen*, Result*, const Xargs& );

    //! TestU01 results type with a custom deleter by TestU01
    using ResultPtr = TestU01Ptr< Result, Deleter >;

  public:
    //! Constructor
    explicit TestU01( tk::ctr::RNGType r,
                      //std::size_t gid,
                      unif01_Gen* gen,
                      std::vector< std::string >&& names,
                      RunFn runner,
                      Ts&&... xargs ) :
      m_id( r ),
      m_gid( 0 ),
      m_gen( gen ),
      m_names( std::move(names) ),
      m_runner( runner ),
      m_xargs( std::forward<Ts>(xargs)... ),
      m_npval( m_names.size() ),
      m_pvals( m_names.size(), -1.0 ),
      m_res( ResultPtr(Creator()) ) {}

    //! Run test
    void run( std::size_t id ) {
      m_pvals = m_runner( m_gen, m_res.get(), m_xargs );
      //m_host->evaluate( id );
    }

    //! Test name accessor
    const std::string& name( std::size_t i ) const {
      Assert( i < m_names.size(), tk::ExceptType::FATAL,
              "Indexing outside of container bounds in TestU01::names()" );
      return m_names[i];
    }

    //! Number of results/test (i.e., p-values) accessor
    std::size_t nstat() const { return m_npval; }

    //! RNG enum accessor
    const tk::ctr::RNGType& rng() const { return m_id; }

    //! RNG id accessor
    std::size_t id() const { return m_gid; }

    //! Query whether test is failed
    bool fail( std::size_t p ) const {
      if ( (m_pvals[p] <= gofw_Suspectp) || (m_pvals[p] >= 1.0-gofw_Suspectp) )
        return true;
      else
        return false;
    }

    //! Return number of failed tests
    std::size_t nfail() const {
      std::size_t f = 0;
      for (std::size_t p = 0; p<m_npval; ++p) if (fail(p)) ++f;
      return f;
    }

    //! p-value accessor
    double pval( std::size_t p ) const { return m_pvals[p]; }

    //! Return human-readable p-value (based on TestU01::bbattery.c::WritePval)
    std::string pvalstr( std::size_t p ) const {
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

    //! Copy assignment
    TestU01& operator=( const TestU01& x) {
      m_id = x.m_id;
      m_gid = x.m_gid;
      m_gen = x.m_gen;
      m_names = x.m_names;
      m_runner = x.m_runner;
      m_xargs = x.m_xargs;
      m_npval = x.m_npval;
      m_pvals = x.m_pvals;
      m_res = ResultPtr(Creator());
      return *this;
    }
    //! Copy constructor: in terms of copy assignment
    TestU01( const TestU01& x ) { operator=(x); }
    //! Move assignment
    TestU01& operator=( TestU01&& ) noexcept = default;
    //! Move constructor
    TestU01( TestU01&& ) noexcept = default;

  private:
    tk::ctr::RNGType m_id;              //!< RNG id
    size_t m_gid;                       //!< global-scope RNG id
    unif01_Gen* m_gen;                  //!< pointer to generator
    std::vector< std::string > m_names; //!< name(s) of tests
    RunFn m_runner;                     //!< test runner function
    Xargs m_xargs;                      //!< extra args for run()
    std::size_t m_npval;                //!< number of p-values produced
    std::vector< double > m_pvals;      //!< p-value(s)
    ResultPtr m_res;                    //!< testU01 results
};

} // rngtest::

#endif // TestU01_h
