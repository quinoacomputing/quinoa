//******************************************************************************
/*!
  \file      src/Walker/Distributor.h
  \author    J. Bakosi
  \date      Wed 15 Jul 2015 08:35:35 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
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

#include <vector>
#include <map>
#include <iosfwd>
#include <cstdint>

#include "Types.h"
#include "Timer.h"
#include "Tags.h"
#include "TaggedTuple.h"
#include "StatCtr.h"
#include "UniPDF.h"
#include "BiPDF.h"
#include "TriPDF.h"
#include "WalkerPrint.h"
#include "Walker/CmdLine/CmdLine.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <sdag.h>
#include "integrator.decl.h"
#include "distributor.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace walker {

//! Distributor drives the time integration of differential equations
class Distributor : public CBase_Distributor {

  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Distributor_SDAG_CODE

  public:
    //! Constructor
    explicit Distributor( const ctr::CmdLine& cmdline );

    //! \brief Reduction target indicating that all Integrator chares have
    //!   registered with the statistics merger (collector)
    //! \details This function is a Charm++ reduction target that is called when
    //!   all Integrator chares have registered with their local branch of the
    //!   statistics merger group, Collector. Once this is done, we issue a
    //!   broadcast to all Itegrator chares to continue with their setup.
    void registered() { m_intproxy.setup( m_dt, m_it, m_moments ); }

    //! Finish initialization
    void init() const;

    //! Estimate ordinary moments
    void estimateOrd( tk::real* ord, std::size_t n );

    //! Estimate central moments
    void estimateCen( tk::real* cen, std::size_t n );

    //! Finish estimation of ordinary PDFs
    void estimateOrdPDF( const std::vector< tk::UniPDF >& updf,
                         const std::vector< tk::BiPDF >& bpdf,
                         const std::vector< tk::TriPDF >& tpdf );

    //! Finish estimation of central PDFs
    void estimateCenPDF( const std::vector< tk::UniPDF >& updf,
                         const std::vector< tk::BiPDF >& bpdf,
                         const std::vector< tk::TriPDF >& tpdf );

    //! Charm++ reduction target enabling shortcutting sync points if no stats
    void nostats();

  private:
    //! Print information at startup
    void info( uint64_t chunksize ) const;

    //! Compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header() const;

    //! Print out one-liner report on time step
    void report();

    //! Output statistics to file
    void outStat();

    //! Write univariate PDF to file
    void writeUniPDF( const tk::UniPDF& p, int& cnt );

    //! Write bivariate PDF to file
    void writeBiPDF( const tk::BiPDF& p, int& cnt );

    //! Write trivariate PDF to file
    void writeTriPDF( const tk::TriPDF& p, int& cnt );

    //! Output PDFs to file
    void outPDF();

    //! Output all requested univariate PDFs to file(s)
    int outUniPDF();

    //! Output all requested bivariate PDFs to file(s)
    int outBiPDF();

    //! Output all requested trivariate PDFs to file(s)
    int outTriPDF();

    //! Evaluate time step, compute new time step size
    void evaluateTime();

    //! Pretty printer
    WalkerPrint m_print;

    //! Counters of integrator chares completing a function
    tk::tuple::tagged_tuple< tag::init,     uint64_t,
                             tag::ordpdf,   uint64_t,
                             tag::cenpdf,   uint64_t,
                             tag::chare,    uint64_t > m_count;

    //! Output indicators
    tk::tuple::tagged_tuple< tag::stat, bool,
                             tag::pdf,  bool > m_output;

    uint64_t m_it;                              //!< Iteration count
    tk::real m_npar;                            //!< Total number of particles
    tk::real m_t;                               //!< Physical time
    tk::real m_dt;                              //!< Physical time step size
    bool m_nostat;                              //!< Any stats to estimate?
    CProxy_Integrator m_intproxy;               //!< Integrator array proxy
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

    //! Map used to lookup moments
    std::map< tk::ctr::Product, tk::real > m_moments;
};

} // walker::

#endif // Distributor_h
