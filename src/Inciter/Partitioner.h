// *****************************************************************************
/*!
  \file      src/Inciter/Partitioner.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare partitioner group used to perform mesh partitioning
  \details   Charm++ chare partitioner group used to perform mesh partitioning.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation, communication as well as I/O. The
    algorithm utilizes the structured dagger (SDAG) Charm++ functionality. The
    high-level overview of the algorithm structure and how it interfaces with
    Charm++ is discussed in the Charm++ interface file
    src/Inciter/partitioner.ci.

    #### Call graph ####
    The following is a directed acyclic graph (DAG) that outlines the
    asynchronous algorithm implemented in this class The detailed discussion of
    the algorithm is given in the Charm++ interface file partitioner.ci,
    which also repeats the graph below using ASCII graphics. On the DAG orange
    fills denote global synchronization points, orange frames with white fill
    are partial synchronization points that overlap with other tasks, and dashed
    lines are potential shortcuts that allow jumping over some of the task-graph
    under some circumstances. See the detailed discussion in partitioner.ci.
    \dot
    digraph "Partitioner SDAG" {
      rankdir="LR";
      node [shape=record, fontname=Helvetica, fontsize=10];
      Own [ label="Own" tooltip="owned nodes reordered"
             URL="\ref inciter::Partitioner::reorder"];
      Req [ label="Req" tooltip="nodes requested"
             URL="\ref inciter::Partitioner::request"];
      Pre [ label="Pre" tooltip="start preparing node IDs"
            URL="\ref inciter::Partitioner::prepare" color="#e6851c"];
      Ord [ label="Ord" tooltip="Node IDs reordered"
            URL="\ref inciter::Partitioner::reordered" color="#e6851c"];
      Low [ label="Low" tooltip="lower bound received"
             URL="\ref inciter::Partitioner::lower"];
      Upp [ label="Upp" tooltip="upper bound computed"
             URL="\ref inciter::Partitioner::bounds"];
      Par [ label="Par" tooltip="partitioners participated"
             URL="\ref inciter::Partitioner::neworder"];
      Cre [ label="Cre" tooltip="create workers"
             URL="\ref inciter::Partitioner::create" color="#e6851c"];
      Own -> Pre [ style="solid" ];
      Req -> Pre [ style="solid" ];
      Pre -> Ord [ style="solid" ];
      Ord -> Low [ style="solid" ];
      Ord -> Upp [ style="solid" ];
      Ord -> Par [ style="solid" ];
      Low -> Cre [ style="solid" ];
      Upp -> Cre [ style="solid" ];
      Par -> Cre [ style="solid" ];
    }
    \enddot
    \include Inciter/partitioner.ci
*/
// *****************************************************************************
#ifndef Partitioner_h
#define Partitioner_h

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <numeric>

#include "ExodusIIMeshReader.h"
#include "ContainerUtil.h"
#include "ZoltanInterOp.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Options/PartitioningAlgorithm.h"
#include "Solver.h"
#include "DerivedData.h"
#include "UnsMesh.h"
#include "AMR/mesh_adapter.h"
#include "FaceData.h"

#include "NoWarning/partitioner.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

//! Partitioner Charm++ chare group class
//! \details Instantiations of Partitioner comprise a processor aware Charm++
//!   chare group. When instantiated, a new object is created on each PE and not
//!   more (as opposed to individual chares or chare array object elements). See
//!   also the Charm++ interface file partitioner.ci.
class Partitioner : public CBase_Partitioner {

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-parameter"
    #pragma clang diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-parameter"
    #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
  #elif defined(__INTEL_COMPILER)
    #pragma warning( push )
    #pragma warning( disable: 1478 )
  #endif
  // Include Charm++ SDAG code. See http://charm.cs.illinois.edu/manuals/html/
  // charm++/manual.html, Sec. "Structured Control Flow: Structured Dagger".
  Partitioner_SDAG_CODE
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #elif defined(__INTEL_COMPILER)
    #pragma warning( pop )
  #endif

  public:
    //! Constructor
    Partitioner( const std::vector< CkCallback >& cb,
                 const CProxy_Transporter& host,
                 const tk::CProxy_Solver& solver,
                 const CProxy_BoundaryConditions& bc,
                 const Scheme& scheme,
                 const std::map< int, std::vector< std::size_t > >& bface,
                 const std::vector< std::size_t >& triinpoel );

    //! Partition the computational mesh
    void partition( int nchare );

    //! Receive number of uniquely assigned global mesh node IDs from lower PEs
    void offset( int p, std::size_t u );

    //! Request new global node IDs for old node IDs
    void request( int p, const std::unordered_set< std::size_t >& nd );

    // Request new global node IDs for edges
    void request( int p, const tk::UnsMesh::Edges& ed );

    //! Receive new (reordered) global node IDs
    void neworder( const std::unordered_map< std::size_t, std::size_t >& nd );

    //! Receive new global node IDs associated to edge-nodes
    void neworder( const tk::UnsMesh::EdgeNodes& ed );

    //! Receive mesh node IDs associated to chares we own
    //! \param[in] n Mesh node indices associated to chare IDs
    //! \param[in] frompe PE call coming from
    void add( int frompe,
              const std::unordered_map< int, std::vector< std::size_t > >& n );

    //! Acknowledge received node IDs
    void recv();

    //! Prepare owned mesh node IDs for reordering
    void flatten();

    //! Receive lower bound of node IDs our PE operates on after reordering
    void lower( std::size_t low );

    //! \brief Compute the variance of the communication cost of merging the
    //!   linear system
    void stdCost( tk::real av );

    //! \brief Start gathering global node IDs this PE will need to receive
    //!   (instead of assign) during reordering
    void gather();

    //! \brief Query our global node IDs and edges by other PEs so they know if
    //!   they are to receive IDs for those from during reordering
    void query( int p, const std::vector< std::size_t >& nodes );

    //! Receive mask of to-be-received global mesh node IDs
    void mask( int p, const std::unordered_map< std::size_t,
                              std::vector< int > >& cn );

    //! Create worker chare array elements on this PE
    void createWorkers();

  private:
    //! Charm++ callbacks associated to compile-time tags
    tk::tuple::tagged_tuple<
        tag::refined,     CkCallback
      , tag::centroid,    CkCallback
      , tag::distributed, CkCallback
      , tag::flattened,   CkCallback
      , tag::load,        CkCallback
      , tag::avecost,     CkCallback
      , tag::stdcost,     CkCallback
      , tag::coord,       CkCallback
    > m_cb;
    //! Host proxy
    CProxy_Transporter m_host;
    //! Linear system solver proxy
    tk::CProxy_Solver m_solver;
    //! Boundary conditions proxy
    CProxy_BoundaryConditions m_bc;
    //! Discretization scheme
    Scheme m_scheme;
    //! Number of fellow PEs to send elem IDs to
    std::size_t m_npe;
    //! Queue of requested node IDs from PEs
    std::vector< std::pair< int, std::unordered_set<std::size_t> > > m_reqNodes;
    //! Queue of requested edge-node IDs from PEs
    std::vector< std::pair< int, tk::UnsMesh::Edges > > m_reqEdges;
    //! \brief Starting global mesh node ID for node reordering on this PE
    //!   during mesh node reordering
    std::size_t m_start;
    //! \brief Counter for number of offsets
    //! \details This counts the to-be-received node IDs received while
    //!   computing global mesh node ID offsets for each PE rquired for node
    //!   reordering later
    std::size_t m_noffset;
    //! \brief Counter for number of queries for global mesh node IDs
    //! \details This counts the number of queries received while
    //!   gathering the node IDs that need to be received (instead of uniquely
    //!   assigned) by each PE
    std::size_t m_nquery;
    //! \brief Counter for number of masks of to-be-received global mesh node
    //!   IDs received
    //! \details This counts the to-be-received node ID masks received while
    //!   gathering the node IDs that need to be received (instead of uniquely
    //!   assigned) by each PE
    std::size_t m_nmask;
    //! Tetrtahedron element connectivity of our chunk of the mesh
    std::vector< std::size_t > m_tetinpoel;
    //! Global element IDs we read (our chunk of the mesh)
    std::vector< long > m_gelemid;
    //! Coordinates of mesh nodes of out chunk of the mesh
    std::array< std::vector< tk::real >, 3 > m_coord;
    //! Element centroid coordinates of our chunk of the mesh
    std::array< std::vector< tk::real >, 3 > m_centroid;
    //! Total number of chares across all PEs
    int m_nchare;
    //! Lower bound of node IDs our PE operates on after reordering
    std::size_t m_lower;
    //! Upper bound of node IDs our PE operates on after reordering
    std::size_t m_upper;
    //! \brief Temporary communication map used to receive global mesh node IDs
    //! \details This map, on each PE, associates the list of global mesh point
    //!   indices to fellow PE IDs from which we will receive new node IDs (as
    //!   in producing contiguous-row-id linear system contributions) during
    //!   reordering.
    std::map< int, std::unordered_set< std::size_t > > m_ncomm;
    //! \brief Temporary communication map used to receive global mesh edges
    //! \details This map, on each PE, associates the list of global mesh edges
    //!   indices to fellow PE IDs from which we will receive new nodes IDs (as
    //!   in producing contiguous-row-id linear system contributions) associated
    //!   to edges during reordering.
    std::map< int, tk::UnsMesh::Edges > m_ecomm;
    //! \brief Communication map used for distributed mesh node reordering
    //! \details This map, on each PE, associates the list of global mesh point
    //!   indices to fellow PE IDs from which we will receive new node IDs (as
    //!   in producing contiguous-row-id linear system contributions) during
    //!   reordering. Only data that will be received from PEs with a lower
    //!   index are stored.
    std::unordered_map< int, std::unordered_set<std::size_t> > m_ncommunication;
    //! \brief Communication map used for distributed mesh edge-node reordering
    //! \details This map, on each PE, associates the list of global mesh edges
    //!   to fellow PE IDs from which we will receive new node IDs (as in
    //!   producing contiguous-row-id linear system contributions)associated to
    //!   edges during reordering. Only data that will be received from PEs with
    //!   a lower index are stored.
    std::unordered_map< int, tk::UnsMesh::Edges > m_ecommunication;
    //! \brief Unique global node IDs chares on our PE will contribute to in a
    //!   linear system
    std::set< std::size_t > m_nodeset;
    //! \brief Unique global edges whose nodes chares on our PE will contribute
    //!   to in a linear system
    tk::UnsMesh::Edges m_edgeset;
    //! \brief Map associating new node IDs (as in producing contiguous-row-id
    //!   linear system contributions) as map-values to old node IDs (as in
    //!   file) as map-keys
    std::unordered_map< std::size_t, std::size_t > m_linnodes;
    //! \brief Map associating new node IDs (as in producing contiguous-row-id
    //!   linear system contributions) as map-values to edges given by two old
    //!   node IDs (as in file) as map-keys
    tk::UnsMesh::EdgeNodes m_linedges;
    //! Global mesh element connectivity associated to chares owned
    std::unordered_map< int, std::vector< std::size_t > > m_chinpoel;
    //! \brief Maps associating old node IDs to new node IDs (as in producing
    //!   contiguous-row-id linear system contributions) categorized by chares.
    //! \details Maps associating old node IDs (as in file) as map-values to new
    //!   node IDs (as in producing contiguous-row-id linear system
    //!   contributions) as map-keys, associated to chare IDs (outer keys).
    std::unordered_map< int,
      std::unordered_map< std::size_t, std::size_t > > m_chfilenodes;
    //! \brief Maps associating new node IDs (as in producing contiguous-row-id
    //!   linear system contributions) to edges (a pair of old node IDs) in
    //!   tk::UnsMesh::EdgeNodes maps, associated to and categorized by chares.
    //! \details Maps associating new node IDs (as in producing
    //!   contiguous-row-id linear system contributions) to edges (a pair of old
    //!   node IDs, as in file) associated to chare IDs (outer key) for only
    //!   the nodes newly added as a result of initial uniform refinement.
    //! \note Used for looking up boundary conditions
    std::unordered_map< int, tk::UnsMesh::EdgeNodes > m_chedgenodes;
    //! Communication cost of linear system merging for our PE
    tk::real m_cost;
    //! \brief Map associating a set of chare IDs to old (as in file) global
    //!   mesh node IDs on the chare boundaries
    //! \details Note that a single global mesh ID can be associated to multiple
    //!   chare IDs as multiple chares can contribute to a single mesh node.
    std::unordered_map< std::size_t, std::vector< int > > m_bnodechares;
    //! \brief Map associating a set of chare IDs to edges given by two old
    //!   global mesh node IDs (old as in file)
    //! \details Note that a single edge can be associated to multiple chare IDs
    //!   as multiple chares can contribute to a single edge.
    tk::UnsMesh::EdgeChares m_edgechares;
    //! \brief Global mesh node IDs associated to chare IDs bordering the mesh
    //!   chunk held by (and associated to) chare IDs this PE owns
    //! \details msum: (M)esh chunks (S)urrounding (M)esh chunks storing mesh
    //!   nodes. Outer map-key: chare IDs this PE owns whose neighbors are
    //!   stored, inner map-key: chare IDs of those chares that hold mesh chunks
    //!   surrounding the outer-key chare's mesh, map-values: global new
    //!   (reordered, as in producing contiguous-row-id linear system
    //!   contributions) mesh node IDs along the border of chares (at which the
    //!   chares will need to communicate) during time stepping.
    std::unordered_map< int,
      std::unordered_map< int, std::unordered_set< std::size_t > > > m_msum;
    //! \brief Mesh edges given by two global mesh node IDs associated to chare
    //!   IDs bordering the mesh chunk held by (and associated to) chare IDs
    //!   this PE owns
    //! \details msum: (M)esh chunks (S)urrounding (M)esh chunks storing mesh
    //!   nodes. Outer map-key: chare IDs this PE owns whose neighbors are
    //!   stored, inner map-key: chare IDs of those chares that hold mesh chunks
    //!   surrounding the outer-key chare's mesh, map-values: mesh edges given
    //!   by two global new (reordered, as in producing contiguous-row-id linear
    //!   system contributions) mesh node IDs along the border of chares (at
    //!   which the chares will need to communicate) during time stepping.
    std::unordered_map< int,
      std::unordered_map< int, tk::UnsMesh::Edges > > m_msumed;
    //! \brief Boundary face list from side-sets.
    //!   m_bface is the list of boundary faces in the side-sets.
    std::map< int, std::vector< std::size_t > > m_bface;
    //! \brief Boundary face-node connectivity.
    std::vector< std::size_t > m_triinpoel;

    //! Compute element centroid coordinates
    void computeCentroids(
      const std::unordered_map< std::size_t, std::size_t >& lid );

    //! Construct global mesh node ids for each chare
    std::unordered_map< int, std::vector< std::size_t > >
    chareNodes( const std::vector< std::size_t >& che ) const;

    //! Distribute global mesh node IDs to their owner PEs
    void distribute( std::unordered_map< int,
                       std::vector< std::size_t > >&& n );

    //! Compute chare distribution
    std::array< int, 2 > chareDistribution() const;

    //! Reorder global mesh node IDs
    void reorder();

    //! Return processing element for chare id
    int pe( int id ) const;

    //! Associate new node IDs to old ones and return them to the requestor(s)
    void prepare();

    //! Generate compact mesh connectivity
    void generate_compact_inpoel();

    //! Uniformly refine our mesh replacing each tetrahedron with 8 new ones
    void refine( const std::vector< std::size_t >& inpoel );

    //! Compute final result of reordering
    void reordered();

    //! Compute lower and upper bounds of reordered node IDs our PE operates on
    void bounds();

    //! \brief Create chare array elements on this PE and assign the global mesh
    //!   element IDs they will operate on
    void create();

    //! Create Discretization chare array elements on this PE
    void createDiscWorkers();

    //! Compute communication cost of linear system merging for our PE
    tk::real cost( std::size_t l, std::size_t u );
};

} // inciter::

#endif // Partitioner_h
