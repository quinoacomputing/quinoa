//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Tue 07 Apr 2015 10:04:43 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Performer advances the Euler equations
  \details   Performer advances the Euler equations. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances the Euler equations in time.
*/
//******************************************************************************
#ifndef Performer_h
#define Performer_h

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <inciter.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <performer.decl.h>
#include <Inciter/InputDeck/InputDeck.h>

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Performer Charm++ chare used to advance the Euler equations in time
template< class Proxy >
class Performer : public CBase_Performer< Proxy > {

  public:
    //! Constructor
    //! \param[in] proxy Host proxy to call back to (here: Conductor)
    explicit Performer( Proxy& proxy ) : m_proxy( proxy )
    {
    }

  private:
    Proxy m_proxy;                      //!< Host proxy
};

} // inciter::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include <performer.def.h>
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // Performer_h
