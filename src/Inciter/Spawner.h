//******************************************************************************
/*!
  \file      src/Inciter/Spawner.h
  \author    J. Bakosi
  \date      Fri 11 Dec 2015 12:40:42 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Charm++ chare spawner group
  \details   Charm++ chare spawner group.
*/
//******************************************************************************
#ifndef Spawner_h
#define Spawner_h

#include <unordered_map>
#include <iostream>     // NOT NEEDED!

#include "Exception.h"
#include "PUPUtil.h"
#include "ContainerUtil.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "spawner.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

//! Spawner Charm++ chare group class
//! \details Instantiations of Spawner comprise a processor aware Charm++ chare
//!   group. When instantiated, a new object is created on each PE and not more
//!   (as opposed to individual chares or chare array object elements). The
//!   group's elements are used to spawn chare objects on every PE. See also the
//!   Charm++ interface file spawner.ci. The class is templated on classes so
//!   that the same code (parameterized by the template arguments) can be
//!   generated for interacting with different types of Charm++ proxies.
//! \see http://charm.cs.illinois.edu/manuals/html/charm++/manual.html
//! \author J. Bakosi
template< class HostProxy, class WorkerProxy, class LinSysMergerProxy >
class Spawner : public CBase_Spawner< HostProxy,
                                      WorkerProxy,
                                      LinSysMergerProxy > {

  private:
    using Group = CBase_Spawner< HostProxy, WorkerProxy, LinSysMergerProxy >;

  public:
    //! Constructor
    //! \param[in] hostproxy Host Charm++ proxy we are being called from
    //! \details The constructor creates an empty chare array of workers right
    //!   away as this object is created on each PE. This initializes the worker
    //!   proxy in m_worker, so it is impossible to call an entry (i.e.,
    //!   asynchronous) member function without m_worker having been
    //!   initialized.
    Spawner( int nchare, HostProxy& host ) :
      m_nchare( nchare ),
      m_host( host ),
      m_worker( WorkerProxy::ckNew() )
    {}

    //! \brief Create chare array elements on this PE and assign the global mesh
    //!   element IDs they will operate on
    //! \param[in] lsm Linear system merger proxy (required by the workers)
    //! \param[in] conn Reordered global mesh connectivity, i.e., node IDs, the
    //!   this PE's chares contribute to in a linear system associated to global
    //!   chare IDs it owns
    //! \param[in] chcid Map associating old node IDs (as in file) to new node
    //!   IDs (as in producing contiguous-row-id linear system contributions)
    //!   categorized by chare IDs (outer key) on this PE
    //! \details We create chare array elements by calling the insert() member
    //!   function, which allows specifying the PE on which the array element is
    //!   created and we send each chare array element the global mesh element
    //!   connectivity, i.e., node IDs, it contributes to and the old->new node
    //!   ID map.
    //! \note It is assumed here that in conn only chare ids that we own is
    //!   coming, i.e., a distribution should already have happened resulting in
    //!   the correct chare ids (and their connectivity) passed to the correct
    //!   PEs.
    void
    create( LinSysMergerProxy& lsm,
            const std::unordered_map< int, std::vector< std::size_t > >& conn,
            const std::unordered_map< int,
              std::unordered_map< std::size_t, std::size_t > >& chcid )
    {
      auto chunksize = m_nchare / CkNumPes();
      auto mynchare = chunksize;
      if (CkMyPe() == CkNumPes()-1) mynchare += m_nchare % CkNumPes();
      for (int c=0; c<mynchare; ++c) {
        // Compute chare ID
        auto cid = CkMyPe() * chunksize + c;
        // Create array element
        m_worker[ c ].insert( cid, m_host, lsm, Group::thisProxy,
                              tk::cref_find(conn,cid), tk::cref_find(chcid,cid),
                              CkMyPe() );
      }
      m_worker.doneInserting();
      Group::contribute(
        CkCallback( CkReductionTarget(Conductor,created), m_host ) );
    }

    //! Instruct all workers to continue with setup mesh data
    void setup() { m_worker.setup(); }

    //! Instruct all workers to perform initialization of the PDE
    //! \param[in] dt Size of time step
    void init( tk::real dt ) { m_worker.init( dt ); }

    //! Instruct all workers to advance PDE
    //! \param[in] stage Stage in multi-stage time stepping
    //! \param[in] dt Size of time step
    //! \param[in] it Iteration count
    //! \param[in] t Physical time
    void advance( uint8_t stage, tk::real dt, uint64_t it, tk::real t )
    { m_worker.advance( stage, dt, it, t ); }

  private:
    int m_nchare;               //!< Total number of chares across all PEs
    HostProxy m_host;           //!< Host proxy
    WorkerProxy m_worker;       //!< Worker proxy
};

} // inciter::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include "spawner.def.h"
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // Spawner_h
