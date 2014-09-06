//******************************************************************************
/*!
  \file      src/Integrator/Distributor.h
  \author    J. Bakosi
  \date      Fri 05 Sep 2014 12:40:29 PM MDT
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
#include <TxtStatWriter.h>
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

    //! Pretty printer
    QuinoaPrint m_print;
    //! Number of integrators completed setting the initial conditions
    uint64_t m_ninit;
    //! Number of integrators completed accumulating the ordinary moments
    uint64_t m_nAccOrd;
    //! Number of integrators completed accumulating the central moments
    uint64_t m_nAccCen;
    //! Number of integrators fired up
    uint64_t m_numchares;
    //! Iteration count
    uint64_t m_it;
    //! Physical time
    tk::real m_t;
    //! Integrator proxies
    std::vector< CProxyInt > m_proxy;
    //! Timers
    std::vector< tk::Timer > m_timer;
    //! Bools indicating whether to plot ordinary moments (in stats output)
    std::vector< bool > m_plotOrdinary;
    //! Ordinary moment names
    std::vector< std::string > m_nameOrdinary;
    //! Central moment names
    std::vector< std::string > m_nameCentral;
    //! Ordinary moments
    std::vector< tk::real > m_ordinary;
    //! Central moments
    std::vector< tk::real > m_central;
    //! Statistics file writer
    TxtStatWriter m_statWriter;
    //! Output indicators
    bool m_wroteStat;
};

} // quinoa::

#endif // Distributor_h
