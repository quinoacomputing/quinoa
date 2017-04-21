// *****************************************************************************
/*!
  \file      src/Walker/Distributor.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Distributor drives the time integration of differential equations
  \details   Distributor drives the time integration of differential equations.
    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation, communication as well as I/O. The
    algorithm utilizes the structured dagger (SDAG) Charm++ functionality. The
    high-level overview of the algorithm structure and how it interfaces with
    Charm++ is discussed in the Charm++ interface file
    src/Walker/distributor.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file distributor.ci,
    which also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points, orange frames with white fill
    are partial synchronization points that overlap with other tasks, and dashed
    lines are potential shortcuts that allow jumping over some of the task-graph
    under some circumstances. See the detailed discussion in distributor.ci.
    \dot
    digraph "Distributor SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      AdvP [ label="AdvP" tooltip="advance particles" style=filled
             color="#c96a00" URL="\ref walker::Integrator::advance"];
      OrdM [ label="OrdM" tooltip="estimate ordinary moments" color="#e6851c"
             URL="\ref walker::Integrator::accumulateOrd"];
      CenM [ label="CenM" tooltip="estimate central moments" color="#e6851c"
             URL="\ref walker::Integrator::accumulateCen"];
      OutS [ label="OutS" tooltip="output statistical moments"
             URL="\ref walker::Distributor::outStat"];
      OrdP [ label="OrdP" tooltip="estimate ordinary PDFs"
             URL="\ref walker::Integrator::accumulateOrdPDF"];
      CenP [ label="CenP" tooltip="estimate central PDFs" color="#e6851c"
             URL="\ref walker::Integrator::accumulateCenPDF"];
      OutP [ label="OutP" tooltip="output PDFs"
             URL="\ref walker::Distributor::outPDF"];
      EvT  [ label="EvT" tooltip="evaluate time step" style=filled
             color="#c96a00" URL="\ref walker::Distributor::evaluateTime"];
      AdvP -> OrdM -> CenM -> OutS -> EvT [ style="solid" ];
      AdvP -> OrdP -> OutP -> EvT [ style="dashed" ];
      OrdM -> CenP -> OutP [ style="dashed" ];
      AdvP -> EvT [ style="dashed" label="NoSt" ];
    }
    \enddot
    \include Walker/distributor.ci
*/
// *****************************************************************************
#ifndef Distributor_h
#define Distributor_h

#include <vector>
#include <map>
#include <iosfwd>
#include <cstdint>

#include "Types.h"
#include "Timer.h"
#include "Tags.h"
#include "Table.h"
#include "TaggedTuple.h"
#include "StatCtr.h"
#include "UniPDF.h"
#include "BiPDF.h"
#include "TriPDF.h"
#include "WalkerPrint.h"
#include "Walker/CmdLine/CmdLine.h"

#include "NoWarning/integrator.decl.h"

namespace walker {

//! Distributor drives the time integration of differential equations
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
    explicit Distributor( const ctr::CmdLine& cmdline );

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
    //! Print information at startup
    void info( uint64_t chunksize, std::size_t nchare );

    //! Compute size of next time step
    tk::real computedt();

    //! Print out time integration header
    void header() const;

    //! Print out one-liner report on time step
    void report();

    //! Output statistics to file
    void outStat();

    //! Write univariate PDF to file
    void writeUniPDF( const tk::UniPDF& p, tk::ctr::Moment m, std::size_t idx );

    //! Write bivariate PDF to file
    void writeBiPDF( const tk::BiPDF& p, tk::ctr::Moment m, std::size_t idx );

    //! Write trivariate PDF to file
    void writeTriPDF( const tk::TriPDF& p, tk::ctr::Moment m, std::size_t idx );

    //! Output PDFs to file
    void outPDF();

    //! Output all requested univariate PDFs to file(s)
    void outUniPDF();

    //! Output all requested bivariate PDFs to file(s)
    void outBiPDF();

    //! Output all requested trivariate PDFs to file(s)
    void outTriPDF();

    //! Evaluate time step, compute new time step size
    void evaluateTime();

    //! Pretty printer
    WalkerPrint m_print;
    //! Output indicators
    tk::tuple::tagged_tuple< tag::stat, bool,
                             tag::pdf,  bool > m_output;
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
               std::vector< tk::Table > > m_tables;

    //! Map used to lookup moments
    std::map< tk::ctr::Product, tk::real > m_moments;

    //! Normal finish of time stepping
    void finish();
};

} // walker::

#endif // Distributor_h
