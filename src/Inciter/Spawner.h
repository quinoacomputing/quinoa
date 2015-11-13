//******************************************************************************
/*!
  \file      src/Inciter/Spawner.h
  \author    J. Bakosi
  \date      Tue 10 Nov 2015 08:56:18 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Charm++ chare spawner group
  \details   Charm++ chare spawner group.
*/
//******************************************************************************
#ifndef Spawner_h
#define Spawner_h

#include <iostream>     // NOT NEEDED!

#include "Exception.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "spawner.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace inciter {

extern std::map< int, std::vector< std::size_t > > g_element;

//! Charm++ chare element map merger reducer
static CkReduction::reducerType ChareElemMapMerger;

//! Spawner Charm++ chare group class
//! \details Instantiations of Spawner comprise a processor aware Charm++ chare
//!   group. When instantiated, a new object is created on each PE and not more
//!   (as opposed to individual chares or chare array object elements). The
//!   group's elements are used to spawn chare objects on every PE. See also the
//!   Charm++ interface file spawner.ci. The class is templated on the host as
//!   well as the work proxy so that the same code (parameterized by the
//!   template arguments) can be generated for interacting with different types
//!   of Charm++ proxies.
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
    //!   asynchronous) member function, such as add() without m_worker having
    //!   been initialized.
    Spawner( int nchare, HostProxy& host ) :
      m_npe( 0 ),
      m_nchare( nchare ),
      m_host( host ),
      m_worker( WorkerProxy::ckNew() )
    {}

    //! \brief Create chare array elements on this PE and assign the global mesh
    //!   element IDs they will operate on
    //! \param[in] lsm Linear system merger proxy (required by the workers)
    //! \details We create chare array elements by calling the insert() member
    //!   function, which allows specifying the PE on which the array element is
    //!   created and we send each chare array element the global mesh element
    //!   IDs it owns. Note that the mesh element IDs we send here and now are
    //!   not necessarily the total number of mesh elements a chare array
    //!   element will work on, since other PEs can also hold element IDs
    //!   assigned to the chares we create here. This is a result of the initial
    //!   mesh graph partitioning executed in parallel on all MPI ranks,
    //!   yielding a distributed data structure storing the global mesh element
    //!   IDs associated to chares in g_element, which thus stores a different
    //!   number and a different set of element IDs associated to chares.
    void create( LinSysMergerProxy& lsm ) {
      auto chunksize = m_nchare / CkNumPes();
      auto mynchare = chunksize;
      if (CkMyPe() == CkNumPes()-1) mynchare += m_nchare % CkNumPes();
      for (int c=0; c<mynchare; ++c) {
        auto cid = CkMyPe() * chunksize + c;    // compute chare ID
        const auto it = g_element.find( cid );  // attempt to find its elements
        std::vector< std::size_t > e;
        if (it != end(g_element)) {             // if found
          e = it->second;                       // extract its elements
          g_element.erase( it );                // remove chare ID and elements
        }
        // create array element (even if e is empty, will receive it later)
        m_worker[ c ].insert( cid, m_host, lsm, Group::thisProxy, e, CkMyPe() );
      }
      m_worker.doneInserting();
      // Construct export map associating those map values (mesh element indices
      // associated to chare IDs) owned by chares we do not create, i.e.,
      // created by fellow PEs
      std::map< int, std::map< int, std::vector< std::size_t > > > exp;
      for (const auto& c : g_element) exp[ pe(c.first) ].insert( c );
      // Export chare IDs and element IDs we do not own to fellow PEs
      m_npe = exp.size();
      for (const auto& p : exp)
        Group::thisProxy[ p.first ].add( CkMyPe(), p.second );
      if (recvaliens()) done();
    }
    //! Receive mesh element indices associated to chares we created
    //! \param[in] element Mesh element indices associated to chare IDs
    //! \details Since the chare IDs the mesh elements are associated to are
    //!   global, i.e., unique across all PEs, we compute the local chare array
    //!   element index here and add the list of element IDs to the chare.
    void add( int frompe,
              const std::map< int, std::vector< std::size_t > >& element )
    {
      for (const auto& c : element) {
        Assert( pe(c.first) == CkMyPe(), "PE " + std::to_string(CkMyPe()) +
                " received a chareid-elemidx-vector pair whose chare it does"
                " not create" );
        m_worker[ c.first - CkMyPe()*m_nchare/CkNumPes() ].add( c.second );
      }
      Group::thisProxy[ frompe ].recv();
    }
    //! Acknowledge received element IDs
    void recv() {
      --m_npe;
      if (recvaliens()) done();
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
    std::size_t m_npe;          //!< Number of fellow PEs to send elem IDs to
    int m_nchare;               //!< Total number of chares across all PEs
    HostProxy m_host;           //!< Host proxy
    WorkerProxy m_worker;       //!< Worker proxy

    //! Return processing element for chare id
    //! \param[in] id Chare id
    //! \return PE that creates the chare
    int pe( int id ) {
      auto pe = id / (m_nchare / CkNumPes());
      if (pe == CkNumPes()) --pe;
      return pe;
    }

    //! Test if all fellow PEs have received my element ids contributions
    bool recvaliens() { return m_npe == 0; }

    //! Signal back to host that we have done our part in sending element IDs
    void done() {
      Group::contribute(
        CkCallback( CkReductionTarget(Conductor,created), m_host ) );
    }
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
