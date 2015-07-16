//******************************************************************************
/*!
  \file      src/Walker/Collector.h
  \author    J. Bakosi
  \date      Wed 15 Jul 2015 08:59:46 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Charm++ module interface file for collecting contributions from
             Integrators
  \details   Charm++ module interface file for collecting contributions from
             Integrators.
*/
//******************************************************************************
#ifndef Collector_h
#define Collector_h

#include <cstddef>

#include "Types.h"
#include "Distributor.h"
#include "Walker/InputDeck/InputDeck.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "collector.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace walker {

extern ctr::InputDeck g_inputdeck;

//! Collector Charm++ chare group class
//! \details Instantiations of Collector comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). The
//!   group's elements are used to collect information from all Integrator chare
//!   objects that happen to be on a given PE. See also the Charm++ interface
//!   file collector.ci.
//! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
//! \author J. Bakosi
class Collector : public CBase_Collector {

  public:
    //! Constructor
    Collector( CProxy_Distributor& hostproxy ) :
      m_hostproxy( hostproxy ),
      m_nchare( 0 ),
      m_nord( 0 ),
      m_ncen( 0 ),
      m_ordinary( g_inputdeck.momentNames( tk::ctr::ordinary ).size(), 0.0 ),
      m_central( g_inputdeck.momentNames( tk::ctr::central ) .size(), 0.0 )
    {}

    //! Chares register on my PE
    //! \note This function does not have to be declared as a Charm++ entry
    //!   method since it is always called by chares on the same PE.
    void checkin() { ++m_nchare; }

    //! Chares contribute ordinary moments
    void chareOrd( const std::vector< tk::real >& ord );

    //! Chares contribute central moments
    void chareCen( const std::vector< tk::real >& cen );

  private:
    CProxy_Distributor m_hostproxy;     //!< Host proxy    
    std::size_t m_nchare;  //!< Number of chares contributing to my PE
    std::size_t m_nord;    //!< Number of chares contributed ordinary moments
    std::size_t m_ncen;    //!< Number of chares contributed central moments
    std::vector< tk::real > m_ordinary; //!< Ordinary moments
    std::vector< tk::real > m_central;  //!< Central moments
};

} // tk::

#endif // Collector_h
