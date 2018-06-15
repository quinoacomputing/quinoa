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

    //! Create worker chare array elements on this PE
    void createWorkers();

    //! Query mesh nodes to identify if they are shared
    void query( int fromch, const std::vector< std::size_t >& nodes );

    //! Receive mesh node IDs found on sender chare
    void mask( int fromch, const std::vector< std::size_t >& found );

//     //! Reduction target estimating the average communication cost among all PEs
//     void aveCost( tk::real c );
// 
//     //! \brief Reduction target estimating the standard deviation of the
//     //!   communication cost among all PEs
//     void stdCost( tk::real c );

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_Sorter::pup(p);
      p | m_host;
      p | m_cbs;
      p | m_ginpoel;
      p | m_coordmap;
      p | m_bface;
      p | m_triinpoel;
      p | m_nchare;
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
    //! ...
    int m_nquery;
    //! ...
    int m_nmask;
    //! ...
    std::unordered_map< int, std::unordered_set< std::size_t > > m_msum;
//     //! Lower bound of node IDs our PE operates on after reordering
//     std::size_t m_lower;
//     //! Upper bound of node IDs our PE operates on after reordering
//     std::size_t m_upper;
    //! \brief Temporary communication map used to receive global mesh node IDs
    //! \details This map, on each PE, associates the list of global mesh point
    //!   indices to fellow PE IDs from which we will receive new node IDs (as
    //!   in producing contiguous-row-id linear system contributions) during
    //!   reordering.
    std::map< int, std::unordered_set< std::size_t > > m_ncomm;
    //! \brief Communication map used for distributed mesh node reordering
    //! \details This map, on each PE, associates the list of global mesh point
    //!   indices to fellow PE IDs from which we will receive new node IDs (as
    //!   in producing contiguous-row-id linear system contributions) during
    //!   reordering. Only data that will be received from PEs with a lower
    //!   index are stored.
    std::unordered_map< int, std::unordered_set<std::size_t> > m_ncommunication;

//     //! Receive lower bound of node IDs our PE operates on after reordering
//     void lower( std::size_t low );
// 
//     //! \brief Compute the variance of the communication cost of merging the
//     //!   linear system
//     void stdCost( tk::real av );
// 
// 
//     //! Reorder global mesh node IDs
//     void reorder();
// 
//     //! Return processing element for chare id
//     int pe( int id ) const;
// 
//     //! Associate new node IDs to old ones and return them to the requestor(s)
//     void prepare();
// 
//     //! Compute final result of reordering
//     void reordered();
// 
//     //! Compute lower and upper bounds of reordered node IDs our PE operates on
//     void bounds();

    //! \brief Create chare array elements on this PE and assign the global mesh
    //!   element IDs they will operate on
    void create();

    //! Create Discretization chare array elements on this PE
    void createDiscWorkers();
};

} // inciter::

#endif // Sorter_h
