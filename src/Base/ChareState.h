// *****************************************************************************
/*!
  \file      src/Base/ChareState.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare state collector group
  \details   Charm++ chare state collectory group used for debugging.
*/
// *****************************************************************************
#ifndef ChareState_h
#define ChareState_h

#include <vector>

#include "Tags.h"
#include "TaggedTuple.h"

#include "NoWarning/charestate.decl.h"

namespace tk {

//! Chare state Charm++ chare group class
//! \details Instantiations of ChareState comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). See
//!   also the Charm++ interface file charestate.ci.
class ChareState : public CBase_ChareState {

  public:
    //! Constructor
    ChareState();

    //! Collect chare state
    void collect() const;

  private:
    //! Chare state
    using State = tk::tuple::tagged_tuple<
                    tag::ch,   std::string     // chare name
                  , tag::id,   int             // thisIndex
                  , tag::it,   uint64_t        // iteration count
                  , tag::fn,   std::string     // member function name
                  >;

    std::vector< State > m_state;       //!< Chare state vector
};

} // tk::

#endif // ChareState_h
