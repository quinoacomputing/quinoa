// *****************************************************************************
/*!
  \file      src/Base/ChareStateCollector.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare state collector group
  \details   Charm++ chare state collectory group used for debugging.
*/
// *****************************************************************************
#ifndef ChareStateCollector_h
#define ChareStateCollector_h

#include <string>

#include "Timer.hpp"
#include "ChareState.hpp"

#include "NoWarning/charestatecollector.decl.h"

namespace tk {

//! Chare state Charm++ chare group class
//! \details Instantiations of ChareStateCollector comprise a processor aware
//!   Charm++ chare group. When instantiated, a new object is created on each
//!   PE and not more (as opposed to individual chares or chare array object
//!   elements). See also the Charm++ interface file charestatecollector.ci.
class ChareStateCollector : public CBase_ChareStateCollector {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Constructor
    //! \details Start timer when constructor is called
    explicit ChareStateCollector() : m_state(), m_timer() {}

    //! Migrate constructor
    explicit ChareStateCollector( CkMigrateMessage* m ) :
      CBase_ChareStateCollector( m ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Insert new state entry
    void insert( const std::string& ch, int id, int pe, uint64_t it,
                 const std::string& fn );

    //! Collect chare state
    void collect( bool error, CkCallback cb );

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \note This is a Charm++ group, pup() is thus only for
    //!    checkpoint/restart.
    void pup( PUP::er &p ) override {
      p | m_state;
      p | m_timer;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i ChareStateCollector object reference
    friend void operator|( PUP::er& p, ChareStateCollector& i ) { i.pup(p); }
    //@}

  private:
    std::vector< ChareState > m_state;  //!< Chare states
    Timer m_timer;                      //!< Timer for getting time stamps
};

} // tk::

#endif // ChareStateCollector_h
