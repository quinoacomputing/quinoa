//******************************************************************************
/*!
  \file      src/Integrator/Distributor.h
  \author    J. Bakosi
  \date      Tue 30 Sep 2014 01:27:48 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations
*/
//******************************************************************************
#ifndef Distributor_h
#define Distributor_h

#include <TaggedTuple.h>
#include <QuinoaPrint.h>
#include <Quinoa/CmdLine/CmdLine.h>
#include <distributor.decl.h>

namespace quinoa {

//! Distributor drives the time integration of differential equations
class Distributor : public CBase_Distributor {

  public:
    //! Constructor
    explicit Distributor( const ctr::CmdLine& cmdline );

    //! Finish initialization
    void init();

    //! Finish estimation of ordinary moments
    void estimateOrd( const std::vector< tk::real >& ord );

    //! Finish estimation of central moments
    void estimateCen( const std::vector< tk::real >& ctr );

    //! Finish estimation of PDFs
    void estimatePDF( const std::vector< UniPDF >& updf,
                      const std::vector< BiPDF >& bpdf,
                      const std::vector< TriPDF >& tpdf );

  private:
    using CProxyInt = CProxy_Integrator< CProxy_Distributor >;

    //! Compute load distribution for given total work and virtualization
    void computeLoadDistribution( uint64_t& chunksize, uint64_t& remainder );

    //! compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header() const;

    //! Evaluate time step
    void evaluateTime();

    //! Print out one-liner report on time step
    void report();

    //! Output univariate PDFs to file(s)
    void outUniPDF();

    //! Output bivariate PDFs to file(s)
    void outBiPDF();

    //! Output trivariate PDFs to file(s)
    void outTriPDF();

    //! Pretty printer
    QuinoaPrint m_print;

    //! Counters of integrator chares completing a function
    tk::tuple::tagged_tuple< tag::init,     uint64_t,
                             tag::ordinary, uint64_t,
                             tag::central,  uint64_t,
                             tag::pdf,      uint64_t,
                             tag::chare,    uint64_t > m_count;

    //! Output indicators
    tk::tuple::tagged_tuple< tag::stat, bool,
                             tag::pdf,  bool > m_output;

    uint64_t m_it;                              //!< Iteration count
    tk::real m_t;                               //!< Physical time
    std::vector< CProxyInt > m_proxy;           //!< Integrator proxies
    std::vector< tk::Timer > m_timer;           //!< Timers
    std::vector< std::string > m_nameOrdinary;  //!< Ordinary moment names
    std::vector< std::string > m_nameCentral;   //!< Central moment names
    std::vector< tk::real > m_ordinary;         //!< Ordinary moments
    std::vector< tk::real > m_central;          //!< Central moments
    std::vector< UniPDF > m_updf;               //!< Univariate PDFs
    std::vector< BiPDF > m_bpdf;                //!< Bivariate PDFs
    std::vector< TriPDF > m_tpdf;               //!< Trivariate PDFs
};

} // quinoa::

#endif // Distributor_h
