//******************************************************************************
/*!
  \file      src/Walker/Distributor.h
  \author    J. Bakosi
  \date      Thu 29 Jan 2015 09:20:18 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations.
    The implementation uses the Charm++ runtime system and is fully asynchronous,
    overlapping computation, communication as well I/O. The algorithm utilizes
    the structured dagger (SDAG) Charm++ functionality. The high-level overview
    of the algorithm structure and how it interfaces with Charm++ is discussed
    in the Charm++ interface file src/Walker/distributor.ci.
*/
//******************************************************************************
#ifndef Distributor_h
#define Distributor_h

#include <TaggedTuple.h>
#include <WalkerPrint.h>
#include <Walker/CmdLine/CmdLine.h>
#include <distributor.decl.h>

namespace walker {

//! Distributor drives the time integration of differential equations
class Distributor : public CBase_Distributor {

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Distributor_SDAG_CODE

  public:
    //! Constructor
    explicit Distributor( const ctr::CmdLine& cmdline );

    //! Finish initialization
    void init();

    //! Finish estimation of ordinary moments
    void estimateOrd( const std::vector< tk::real >& ord );

    //! Finish estimation of central moments
    void estimateCen( const std::vector< tk::real >& ctr );

    //! Finish estimation of ordinary PDFs
    void estimateOrdPDF( const std::vector< tk::UniPDF >& updf,
                         const std::vector< tk::BiPDF >& bpdf,
                         const std::vector< tk::TriPDF >& tpdf );

    //! Finish estimation of central PDFs
    void estimateCenPDF( const std::vector< tk::UniPDF >& updf,
                         const std::vector< tk::BiPDF >& bpdf,
                         const std::vector< tk::TriPDF >& tpdf );

  private:
    using CProxyInt = CProxy_Integrator< CProxy_Distributor >;

    //! Print information at startup
    void info( uint64_t chunksize, uint64_t remainder ) const;

    //! Compute load distribution for given total work and virtualization
    void computeLoadDistribution( uint64_t& chunksize, uint64_t& remainder );

    //! Compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header() const;

    //! Print out one-liner report on time step
    void report();

    //! Output statistics to file
    void outStat();

    //! Write univariate PDF to file
    void writeUniPDF( const tk::UniPDF& p, std::size_t& cnt );

    //! Write bivariate PDF to file
    void writeBiPDF( const tk::BiPDF& p, std::size_t& cnt );

    //! Write trivariate PDF to file
    void writeTriPDF( const tk::TriPDF& p, std::size_t& cnt );

    //! Output PDFs to file
    void outPDF();

    //! Output all requested univariate PDFs to file(s)
    std::size_t outUniPDF();

    //! Output all requested bivariate PDFs to file(s)
    std::size_t outBiPDF();

    //! Output all requested trivariate PDFs to file(s)
    std::size_t outTriPDF();

    //! Evaluate time step, compute new time step size
    void evaluateTime();

    //! Pretty printer
    WalkerPrint m_print;

    //! Counters of integrator chares completing a function
    tk::tuple::tagged_tuple< tag::init,     uint64_t,
                             tag::ordinary, uint64_t,
                             tag::central,  uint64_t,
                             tag::ordpdf,   uint64_t,
                             tag::cenpdf,   uint64_t,
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
    std::vector< tk::UniPDF > m_ordupdf;        //!< Ordinary univariate PDFs
    std::vector< tk::BiPDF > m_ordbpdf;         //!< Ordinary bivariate PDFs
    std::vector< tk::TriPDF > m_ordtpdf;        //!< Ordinary trivariate PDFs
    std::vector< tk::UniPDF > m_cenupdf;        //!< Central univariate PDFs
    std::vector< tk::BiPDF > m_cenbpdf;         //!< Central bivariate PDFs
    std::vector< tk::TriPDF > m_centpdf;        //!< Central trivariate PDFs
};

} // walker::

#endif // Distributor_h
