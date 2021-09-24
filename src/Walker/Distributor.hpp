// *****************************************************************************
/*!
  \file      src/Walker/Distributor.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations.
    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation, communication as well as I/O. The
    algorithm utilizes the structured dagger (SDAG) Charm++ functionality. The
    high-level overview of the algorithm structure and how it interfaces with
    Charm++ is discussed in the Charm++ interface file
    src/Walker/distributor.ci.
*/
// *****************************************************************************
#ifndef Distributor_h
#define Distributor_h

#include <vector>
#include <map>
#include <iosfwd>
#include <cstdint>

#include "Types.hpp"
#include "Timer.hpp"
#include "Tags.hpp"
#include "Table.hpp"
#include "TaggedTuple.hpp"
#include "StatCtr.hpp"
#include "UniPDF.hpp"
#include "BiPDF.hpp"
#include "TriPDF.hpp"
#include "WalkerPrint.hpp"
#include "Walker/CmdLine/CmdLine.hpp"

#include "NoWarning/integrator.decl.h"

namespace walker {

//! Distributor drives the time integration of differential equations
// cppcheck-suppress noConstructor
class Distributor : public CBase_Distributor {

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-parameter"
    #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-parameter"
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(__INTEL_COMPILER)
    #pragma warning( push )
    #pragma warning( disable: 1478 )
  #endif
  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Distributor_SDAG_CODE
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #elif defined(__INTEL_COMPILER)
    #pragma warning( pop )
  #endif

  public:
    //! Constructor
    explicit Distributor();

    //! \brief Reduction target indicating that all Integrator chares have
    //!   registered with the statistics merger (collector)
    //! \details This function is a Charm++ reduction target that is called when
    //!   all Integrator chares have registered with their local branch of the
    //!   statistics merger group, Collector. Once this is done, we issue a
    //!   broadcast to all Itegrator chares to continue with their setup.
    void registered() { m_intproxy.setup( m_dt, m_t, m_it, m_moments ); }

    //! Estimate ordinary moments
    void estimateOrd( tk::real* ord, int n );

    //! Estimate central moments
    void estimateCen( tk::real* cen, int n );

    //! Estimate ordinary PDFs
    void estimateOrdPDF( CkReductionMsg* msg );

    //! Finish estimation of central PDFs
    void estimateCenPDF( CkReductionMsg* msg );

    //! Charm++ reduction target enabling shortcutting sync points if no stats
    void nostat();

  private:
    //! Type alias for output indicators
    using OutputIndicators = tk::TaggedTuple< brigand::list<
                                 tag::stat,      bool
                               , tag::pdf,       bool
                             > >;

    OutputIndicators m_output;                  //!< Output indicators
    uint64_t m_it;                              //!< Iteration count
    tk::real m_npar;                            //!< Total number of particles
    tk::real m_t;                               //!< Physical time
    tk::real m_dt;                              //!< Physical time step size
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
    //! Names of and tables to sample and output to statistics file
    std::pair< std::vector< std::string >,
               std::vector< tk::Table<1> > > m_tables;
    //! Map used to lookup moments
    std::map< tk::ctr::Product, tk::real > m_moments;

    //! Print information at startup
    void info( const WalkerPrint& print,
               uint64_t chunksize,
               std::size_t nchare );

    //! Compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header( const WalkerPrint& print ) const;

    //! Print out one-liner report on time step
    void report();

    //! Output statistics to file
    void outStat();

    //! Write univariate PDF to file
    void writeUniPDF( std::uint64_t it,
                      tk::real t,
                      const tk::UniPDF& p,
                      tk::ctr::Moment m,
                      std::size_t idx );

    //! Write bivariate PDF to file
    void writeBiPDF( std::uint64_t it, tk::real t,
                     const tk::BiPDF& p,
                     tk::ctr::Moment m,
                     std::size_t idx );

    //! Write trivariate PDF to file
    void writeTriPDF( std::uint64_t it,
                      tk::real t,
                      const tk::TriPDF& p,
                      tk::ctr::Moment m,
                      std::size_t idx );

    //! Output PDFs to file
    void outPDF();

    //! Output all requested univariate PDFs to file(s)
    void outUniPDF( std::uint64_t it, tk::real t );

    //! Output all requested bivariate PDFs to file(s)
    void outBiPDF( std::uint64_t it, tk::real t );

    //! Output all requested trivariate PDFs to file(s)
    void outTriPDF( std::uint64_t it, tk::real t );

    //! Evaluate time step, compute new time step size
    void evaluateTime();

    //! Create pretty printer specialized to Walker
    //! \return Pretty printer
    WalkerPrint printer() const {
      const auto& def =
        g_inputdeck_defaults.get< tag::cmd, tag::io, tag::screen >();
      auto nrestart = g_inputdeck.get< tag::cmd, tag::io, tag::nrestart >();
      return WalkerPrint(
        g_inputdeck.get< tag::cmd >().logname( def, nrestart ),
        g_inputdeck.get< tag::cmd, tag::verbose >() ? std::cout : std::clog,
        std::ios_base::app );
    }

    //! Normal finish of time stepping
    void finish();
};

} // walker::

#endif // Distributor_h
