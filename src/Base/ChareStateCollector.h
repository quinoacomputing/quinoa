// *****************************************************************************
/*!
  \file      src/Base/ChareStateCollector.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare state collector group
  \details   Charm++ chare state collectory group used for debugging.
*/
// *****************************************************************************
#ifndef ChareStateCollector_h
#define ChareStateCollector_h

#include <string>

#include "Timer.h"
#include "ChareState.h"

#include "NoWarning/charestatecollector.decl.h"

namespace tk {

//! Chare state Charm++ chare group class
//! \details Instantiations of ChareStateCollector comprise a processor aware
//!   Charm++ chare group. When instantiated, a new object is created on each
//!   PE and not more (as opposed to individual chares or chare array object
//!   elements). See also the Charm++ interface file charestatecollector.ci.
class ChareStateCollector : public CBase_ChareStateCollector {

  public:
    //! Constructor
    //! \details Start timer when constructor is called
    ChareStateCollector() : m_state(), m_timer() {}

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Insert new state entry
    void insert( const std::string& ch, int id, int pe, uint64_t it,
                 const std::string& fn );

    //! Collect chare state
    void collect( bool error, CkCallback cb );

  private:
    std::vector< ChareState > m_state;  //!< Chare states
    Timer m_timer;                      //!< Timer for getting time stamps
};

} // tk::

#endif // ChareStateCollector_h
