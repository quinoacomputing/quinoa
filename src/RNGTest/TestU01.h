//******************************************************************************
/*!
  \file      src/RNGTest/TestU01.h
  \author    J. Bakosi
  \date      Sun 25 May 2014 05:41:42 PM MDT
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

#include <Base.h>
#include <RNG.h>
#include <TestU01Util.h>
#include <StatTest.h>
#include <PUPUtil.h>
#include <TestU01PUP.h>
#include <testu01.decl.h>

//! PUP operator for TestU01's unif01_Gen struct as POD
PUPbytes( unif01_Gen )

namespace rngtest {

//! TestU01 : StatTest
template< class Result,                 //!< Results type
          Result* (*Creator)(void),     //!< Results creator
          void (*Deleter)(Result *),    //!< Results deleter
          typename... Ts >              //!< Extra runner args types
class TestU01 : public StatTest {

    //! Test extra arguments type
    using Xargs = std::tuple< Ts... >;

    //! Test runner function pointer type
    using RunFn = std::vector<double> (*)( unif01_Gen*, Result*, const Xargs& );

    //! TestU01 results type with a custom deleter by TestU01
    using ResultPtr = TestU01Ptr< Result, Deleter >;

  public:
    PUPable_decl_template( TestU01 );

    //! Constructor
    explicit TestU01( std::vector< std::string >&& names,
                      RunFn runner,
                      Ts&&... xargs ) :
      m_npval( names.size() ),
      m_names( std::move(names) ),
      m_runner( runner ),
      m_xargs( std::forward<Ts>(xargs)... ),
      m_pvals( names.size(), -1.0 ),
      m_res( ResultPtr(Creator()) )
    {
      //initRNG( base ),
      PUPable_reg( TestU01 );
    }

    //! Migrator
    TestU01( CkMigrateMessage* m = 0 ) {}

    //! Destructor
    ~TestU01() override = default;

//     //! Pack/Unpack RunFn (function pointer)
//     void PUP_RunFn( PUP::er& p, RunFn r ) {
//     }
// 
    //! Pack/Unpack
    void pup( PUP::er& p ) override {
//       //p | *m_gen;       //! Test if this is okay as POD
//       p | m_rng;
//       p | m_npval;
//       PUP_vector( p, m_names );
//       PUP_RunFn( p, m_runner );
//       p | m_xargs;
//       PUP_vector( p, m_pvals );
//       PUP_unique_ptr( p, m_res, ResultPtr(Creator()) );
    }

    //! Run test
    void run( std::size_t id ) override {
      std::cout << "TestU01::run( " << id << " )" << std::endl;
      //std::cout << m_gen << std::endl;
      // TestU01 does not guarantee the constness of gen
      //m_pvals = m_runner( const_cast<unif01_Gen*>(m_gen), m_res.get(), m_xargs );
      //m_host->evaluate( id );
    }

    //! Test name accessor
    const std::string& name( std::size_t i ) const override {
      Assert( i < m_names.size(), tk::ExceptType::FATAL,
              "Indexing outside of container bounds in TestU01::names()" );
      return m_names[i];
    }

    //! Number of results/test (i.e., p-values) accessor
    std::size_t nstat() const override { return m_npval; }

    //! RNG enum accessor
    const tk::ctr::RNGType& rng() const override { return m_id; }

    //! RNG id accessor
    std::size_t id() const override { return m_gid; }

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

    std::size_t m_npval;                //!< Number of p-values produced
    std::vector< std::string > m_names; //!< Name(s) of tests
    RunFn m_runner;                     //!< Test runner function
    Xargs m_xargs;                      //!< Extra args for run()
    std::vector< double > m_pvals;      //!< p-value(s)
    ResultPtr m_res;                    //!< TestU01 results

    std::size_t m_gid;                  //!< RNG global id
    tk::ctr::RNGType m_id;              //!< RNG selected
    Gen01Ptr m_gen;                     //!< Pointer to generator
};

} // rngtest::

#endif // TestU01_h
