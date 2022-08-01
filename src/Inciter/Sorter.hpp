// *****************************************************************************
/*!
  \file      src/Inciter/Sorter.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "TaggedTuple.hpp"
#include "Tags.hpp"
#include "Callback.hpp"
#include "UnsMesh.hpp"
#include "UnsMesh.hpp"
#include "Scheme.hpp"
#include "CommMap.hpp"

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
    explicit Sorter( std::size_t meshid,
                     const CProxy_Transporter& transporter,
                     const tk::CProxy_MeshWriter& meshwriter,
                     const tk::SorterCallback& cbs,
                     const std::vector< Scheme >& scheme,
                     CkCallback reorderRefiner,
                     const std::vector< std::size_t >& ginpoel,
                     const tk::UnsMesh::CoordMap& coordmap,
                     const tk::UnsMesh::Chunk& el,
                     const std::map< int, std::vector< std::size_t > >& bface,
                     const std::vector< std::size_t >& triinpoel,
                     const std::map< int, std::vector< std::size_t > >& bnode,
                     const std::unordered_map< int, std::set< std::size_t > >&
                       elemblockid,
                     int nchare );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVarPrivate
    explicit Sorter( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Setup chare mesh boundary node communication map
    void setup( std::size_t npoin );
    //! \brief Incoming query for a list mesh nodes for which this chare
    //!   compiles communication maps
    void query( int fromch, const tk::AllCommMaps& bnd );
    //! Report receipt of boundary node lists
    void recvquery();
    //! Respond to boundary node list queries
    void response();
    //! Receive boundary node communication maps for our mesh chunk
    void bnd( int fromch, const tk::CommMaps& msum );
    //! Receive receipt of boundary node communication map
    void recvbnd();

    //! Start reordering (if user enabled it)
    void start();

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

    //! Create worker chare array elements on this PE
    void createWorkers();

    //! Update mesh data we hold for whoever calls this function
    void mesh( std::vector< std::size_t >& ginpoel,
               tk::UnsMesh::CoordMap& coordmap,
               std::vector< std::size_t >& triinpoel,
               std::map< int, std::vector< std::size_t > >& bnode );

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_meshid;
      p | m_host;
      p | m_meshwriter;
      p | m_cbs;
      p | m_scheme;
      p | m_reorderRefiner;
      p | m_ginpoel;
      p | m_coordmap;
      p | m_el;
      p | m_nbnd;
      p | m_bface;
      p | m_triinpoel;
      p | m_bnode;
      p | m_elemblockid;
      p | m_nchare;
      p | m_nodeset;
      p | m_noffset;
      p | m_nodech;
      p | m_chnode;
      p | m_edgech;
      p | m_chedge;
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
    //! Mesh ID
    std::size_t m_meshid;
    //! Host proxy
    CProxy_Transporter m_host;
    //! MeshWriter proxy
    tk::CProxy_MeshWriter m_meshwriter;
    //! Charm++ callbacks associated to compile-time tags for sorter
    tk::SorterCallback m_cbs;
    //! Discretization schemes (one per mesh)
    std::vector< Scheme > m_scheme;
    //! Callback to use to send reordered mesh to Refiner
    CkCallback m_reorderRefiner;
    //! Tetrtahedron element connectivity of our chunk of the mesh (global ids)
    std::vector< std::size_t > m_ginpoel;
    //! Coordinates associated to global node IDs of our mesh chunk
    tk::UnsMesh::CoordMap m_coordmap;
    //! Elements of the mesh chunk we operate on
    tk::UnsMesh::Chunk m_el;
    //! Counter for number of chares contributing to chare boundary nodes
    std::size_t m_nbnd;
    //! List of boundary faces associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bface;
    //! Boundary face-node connectivity
    std::vector< std::size_t > m_triinpoel;
    //! List of boundary nodes associated to side-set IDs
    std::map< int, std::vector< std::size_t > > m_bnode;
    //! Local tet ids associated with mesh block ids
    std::unordered_map< int, std::set< std::size_t > > m_elemblockid;
    //! Total number of sorter chares
    int m_nchare;
    //! Unique global node IDs chares on our PE will contribute to
    std::set< std::size_t > m_nodeset;
    //! \brief Counter for the number of chares from which this chare has
    //!   received node reordering offsets from
    int m_noffset;
    //! Node->chare map used to build boundary node communication maps
    std::unordered_map< std::size_t, std::vector< int > > m_nodech;
    //! Chare->node map used to build boundary node communication maps
    tk::NodeCommMap m_chnode;
    //! Edge->chare map used to build boundary edge communication maps
    std::unordered_map< tk::UnsMesh::Edge, std::vector< int >,
                        tk::UnsMesh::Hash<2>, tk::UnsMesh::Eq<2> > m_edgech;
    //! Chare->edge map used to build boundary edge communication maps
    tk::EdgeCommMap m_chedge;
    //! Communication maps associated to chare IDs
    tk::CommMaps m_msum;
    //! \brief Communication map used for distributed mesh node reordering
    //! \details This map associates the list of global mesh point
    //!   indices to fellow chare IDs from which this chare receives new node
    //!   IDs during reordering. Only data that will be received from chares
    //!   with a lower index are stored, thus this is an asymmetric
    //!   communication map.
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

    //! Start preparing for mesh node reordering in parallel
    void mask();

    //! Reorder global mesh node IDs
    void reorder();

    //! Associate new node IDs to old ones and return them to the requestor(s)
    void prepare();

    //! Compute final result of reordering
    void finish();

    //! Create Discretization chare array elements on this PE
    void createDiscWorkers();
};

} // inciter::

#endif // Sorter_h
