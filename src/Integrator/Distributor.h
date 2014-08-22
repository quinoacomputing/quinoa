//******************************************************************************
/*!
  \file      src/Integrator/Distributor.h
  \author    J. Bakosi
  \date      Fri 15 Aug 2014 10:18:47 AM MDT
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

    //! Finish a step in time
    void step();

    //! Finish time stepping
    void finish();

  private:
    using CProxyInt = CProxy_Integrator< CProxy_Distributor >;

    //! Compute load distribution for given total work and virtualization
    void computeLoadDistribution( uint64_t& chunksize, uint64_t& remainder );

    //! Print out time integration header
    void header() const;

    //! Print out one-liner report on time step
    void report();

    QuinoaPrint m_print;                //!< Pretty printer
    uint64_t m_ninit;                   //!< Number of integrators complete init
    uint64_t m_nstep;                   //!< Number of integrators complete step
    uint64_t m_nfinish;                 //!< Number of integrators finish steps
    uint64_t m_it;                      //!< Iteration count
    tk::real m_t;                       //!< Time
    uint64_t m_numchares;               //!< Number of integrators to fire up
    std::vector< CProxyInt > m_proxy;   //!< Integrator proxies
    std::vector< tk::Timer > m_timer;   //!< Timers
};

} // quinoa::

#endif // Distributor_h
