// *****************************************************************************
/*!
  \file      src/Inciter/Sorter.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Mesh sorter for global distributed mesh reordering
  \details   Mesh sorter is Charm++ chare array and is used to do global
    distributed mesh node reordering that yields consecutive unique global node
    IDs with increasing PE IDs in asynchronous distributed-memory parallel
    fashion.
*/
// *****************************************************************************
#ifndef Sorter_h
#define Sorter_h

#include <vector>
#include <map>
#include <unordered_map>

#include "TaggedTuple.h"
#include "Tags.h"
#include "Callback.h"
#include "UnsMesh.h"
#include "UnsMesh.h"
#include "Scheme.h"

#include "NoWarning/transporter.decl.h"
#include "NoWarning/sorter.decl.h"

namespace inciter {

//! Mesh sorter for global distributed mesh node reordering
class Sorter : public CBase_Sorter {

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
  Sorter_SDAG_CODE
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #elif defined(__INTEL_COMPILER)
    #pragma warning( pop )
  #endif

  public:
    //! Constructor
    explicit Sorter( const CProxy_Transporter& transporter,
                     const tk::CProxy_Solver& solver,
                     const tk::SorterCallback& cbs,
                     const Scheme& scheme,
                     const std::vector< std::size_t >& ginpoel,
                     const tk::UnsMesh::CoordMap& coordmap,
                     const std::map< int, std::vector< std::size_t > >& bface,
                     const std::vector< std::size_t >& triinpoel,
                     const std::map< int, std::vector< std::size_t > >& bnode,
                     int nchare );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit Sorter( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Create worker chare array elements on this PE
    void createWorkers();

    //! Receive aggregated chare boundary nodes associated to chares
    void comChBndNode( CkReductionMsg* msg );

    //! Receive aggregated boundary nodes for side sets
    void comnode( CkReductionMsg* msg );

    //! \brief Receive aggregated boundary faces (and triangle connectivity)
    //!    for side sets
    void comface( CkReductionMsg* msg );

    //! \brief Receive number of uniquely assigned global mesh node IDs from
    //!   chares with lower indices
    void offset( int c, std::size_t u );

    //! Request new global node IDs for old node IDs
    void request( int c, const std::unordered_set< std::size_t >& nd );

    //! Receive lower bound of node IDs our PE operates on after reordering
    void lower( std::size_t low );

    //! Receive new (reordered) global node IDs and coordinates
    void neworder( const std::unordered_map< std::size_t,
           std::tuple< std::size_t, tk::UnsMesh::Coord > >& nodes );

    //! Create Discretization chare array elements on this PE
    void createDiscWorkers();

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_host;
      p | m_solver;
      p | m_cbs;
      p | m_scheme;
      p | m_ginpoel;
      p | m_coordmap;
      p | m_bface;
      p | m_triinpoel;
      p | m_bnode;
      p | m_nchare;
      p | m_nodeset;
      p | m_noffset;
      p | m_msum;
      p | m_reordcomm;
      p | m_start;
      p | m_newnodes;
      p | m_newcoordmap;
      p | m_reqnodes;
      p | m_lower;
      p | m_upper;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] s Sorter object reference
    friend void operator|( PUP::er& p, Sorter& s ) { s.pup(p); }
    //@}

  private:
    //! Host proxy
    CProxy_Transporter m_host;
    //! Linear system solver proxy
    tk::CProxy_Solver m_solver;
    //! Charm++ callbacks associated to compile-time tags for sorter
    tk::SorterCallback m_cbs;
    //! Discretization scheme
    Scheme m_scheme;
    //! Tetrtahedron element connectivity of our chunk of the mesh (global ids)
    std::vector< std::size_t > m_ginpoel;
    //! Coordinates associated to global node IDs of our mesh chunk
    tk::UnsMesh::CoordMap m_coordmap;
    //! List of boundary faces associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary face-node connectivity
    std::vector< std::size_t > m_triinpoel;
    //! List of boundary nodes associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Total number of sorter chares
    int m_nchare;
    //! Unique global node IDs chares on our PE will contribute to
    std::set< std::size_t > m_nodeset;
    //! \brief Counter for the number of chares from which this chare has
    //!   received node reordering offsets from
    int m_noffset;
    //! \brief Symmetric communication map associating global mesh node IDs to
    //!    chare IDs this chare shares the node IDs with
    std::map< int, std::unordered_set< std::size_t > > m_msum;
    //! \brief Communication map used for distributed mesh node reordering
    //! \details This map associates the list of global mesh point
    //!   indices to fellow chare IDs from which this chare receives new node IDs
    //!   during reordering. Only data that will be received from chares with a
    //!   lower index are stored, thus this is an asymmetric communication map.
    std::unordered_map< int, std::unordered_set< std::size_t > > m_reordcomm;
    //! \brief Starting global mesh node ID for node reordering on this chare
    //!   during mesh node reordering
    std::size_t m_start;
    //! Map associating new node IDs (value) to old node IDs (key)
    std::unordered_map< std::size_t, std::size_t > m_newnodes;
    //! Coordinates associated to global (new) node IDs during reordering
    tk::UnsMesh::CoordMap m_newcoordmap;
    //! Queue of requested node IDs from chares
    std::vector< std::pair< int, std::unordered_set<std::size_t> > > m_reqnodes;
    //! Lower bound of node IDs this chare contributes to in a linear system
    std::size_t m_lower;
    //! Upper bound of node IDs this chare contributes to in a linear system
    std::size_t m_upper;

    //! Start reordering (if user enabled it)
    void start();

    //! Start preparing for mesh node reordering in parallel
    void mask();

    //! Reorder global mesh node IDs
    void reorder();

    //! Associate new node IDs to old ones and return them to the requestor(s)
    void prepare();

    //! Compute final result of reordering
    void finish();

    //! Compute lower and upper bounds of reordered node IDs our PE operates on
    void bounds();

    //! \brief Create chare array elements on this PE and assign the global mesh
    //!   element IDs they will operate on
    void create();
};

} // inciter::

#endif // Sorter_h
