//******************************************************************************
/*!
  \file      src/Main/InciterDriver.C
  \author    J. Bakosi
  \date      Sun 08 Nov 2015 01:03:21 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************
#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "performer.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include "InciterDriver.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern CProxy_Conductor g_ConductorProxy;

class InciterPrint;

} // inciter::

using inciter::InciterDriver;

InciterDriver::InciterDriver( const InciterPrint& print ) : m_print( print )
//******************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)
}

void
InciterDriver::execute( std::size_t npoin, uint64_t nchare ) const
//******************************************************************************
//  Execute: Fire up Conductor, a Charm++ chare that conducts Performers
//! \param[in] npoin Total number of points in computational mesh
//! \param[in] nchare Total number of chares
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate Conductor chare which drives the time-integration of a PDE via
  // several Performer chares. Charm++ chare Conductor fires up performers.
  // Store proxy handle in global-scope to make it available to individual
  // Performers so they can call back to Conductor. Since this is called inside
  // the main chare constructor, the Charm++ runtime system distributes the
  // handle along with all other global-scope data.
  g_ConductorProxy = CProxy_Conductor::ckNew( npoin, nchare );
}
