// *****************************************************************************
/*!
  \file      src/Inciter/Refiner.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Mesh refiner for interfacing the mesh refinement library
  \details   Mesh refiner is a Charm++ chare array and is used to interface the
   mesh refinement object which does not know about parallelization and thus the
   distributed nature of the mesh it operates on, i.e., it operates on mesh
   chunks. Thus it does not do parallel communication and also does not know
   about global vs local IDs. Instead this Charm++ chare array is the one that
   does all parallel computing aspects, i.e., communcation, and using the mesh
   refiner object as a library.
*/
// *****************************************************************************
#ifndef Refiner_h
#define Refiner_h

#include <vector>

#include "AMR/mesh_adapter.h"
#include "Inciter/Options/AMRInitial.h"
#include "TaggedTuple.h"
#include "Tags.h"
#include "Callback.h"

#include "NoWarning/transporter.decl.h"
#include "NoWarning/refiner.decl.h"

namespace inciter {

//! Mesh refiner for interfacing the mesh refinement library
class Refiner : public CBase_Refiner {

  public:
    //! Constructor
    explicit Refiner( const CProxy_Transporter& transporter,
                      const tk::RefinerCallback& cbr,
                      const std::vector< std::size_t >& ginpoel,
                      const tk::UnsMesh::CoordMap& coordmap );

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_Refiner::pup(p);
      p | m_host;
      p | m_cbr;
      p | m_el;
      if (p.isUnpacking()) {
        m_inpoel = std::get< 0 >( m_el );
        m_gid = std::get< 1 >( m_el );
        m_lid = std::get< 2 >( m_el );
      }
      p | m_initref;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] r Refiner object reference
    friend void operator|( PUP::er& p, Refiner& r ) { r.pup(p); }
    //@}

  private:
    //! Host proxy
    CProxy_Transporter m_host;
    //! Charm++ callbacks associated to compile-time tags for refiner
    tk::RefinerCallback m_cbr;
    //! Elements of the mesh chunk we operate on
    //! \details The first vector is the element connectivity (local IDs), the
    //!   second vector is the global node IDs of owned elements, while the
    //!   third one is a map of global->local node IDs.
    std::tuple< std::vector< std::size_t >,
                std::vector< std::size_t >,
                std::unordered_map< std::size_t, std::size_t > > m_el;
    //! Alias to element connectivity with local node IDs in m_el
    std::vector< std::size_t >& m_inpoel = std::get<0>( m_el );
    //! Alias to global node IDs of owned elements in m_el
    std::vector< std::size_t >& m_gid = std::get<1>( m_el );
    //! \brief Alias to local node IDs associated to the global ones of owned
    //!    elements in m_el
    std::unordered_map< std::size_t, std::size_t >& m_lid = std::get<2>( m_el );
    //! Initial mesh refinement type list (in reverse order)
    std::vector< ctr::AMRInitialType > m_initref;

//     //! Mesh refiner object
//     AMR::mesh_adapter_t m_refiner;
//     //! \brief Counter during distribution of PE-boundary edges during initial
//     //!   mesh refinement
//     std::size_t m_nedge;
//     //! Counter during distribution of newly added nodes to PE-boundary edges
//     std::size_t m_nref;
//     std::size_t m_extra;
//     //! PEs we share at least a single edge with during initial mesh refinement
//     std::unordered_set< int > m_pe;
//     //! \brief Map associating the global IDs and the coordinates of a node
//     //!   added to an edge during initial mesh refinement
//     tk::UnsMesh::EdgeNodeCoord m_edgenode;
//     //! Unique set of boundary edges associated to PEs we share these edges with
//     std::unordered_map< int, tk::UnsMesh::EdgeSet > m_bndEdges;
//     //! \brief Map associating the global IDs and the coordinates of a node
//     //!   added to an edge during initial mesh refinement associated to
//     //!   a(nother) PE the edge is shared with
//     std::unordered_map< int, tk::UnsMesh::EdgeNodeCoord > m_edgenodePe;
// 

//     //! Generate boundary edges and send them to all PEs
//     void bndEdges();
// 
//     //! Receive boundary edges from all PEs (including this one)
//     void addBndEdges( int frompe, const tk::UnsMesh::EdgeSet& ed );
// 
//     //! Receive newly added mesh node IDs on our PE boundary
//     void addRefBndEdges( int frompe, const tk::UnsMesh::EdgeNodeCoord& ed );
// 
//     //! \brief Acknowledge received newly added node IDs to edges shared among
//     //!   multiple PEs
//     void recvRefBndEdges();
// 
//     //! Correct refinement to arrive at a conforming mesh across PE boundaries
//     void correctref();
// 
//     //! Decide wether to continue with another step of initial mesh refinement
//     void nextref();
// 
//     //! Partition the mesh before a (potential) refinement step
//     void partref();
// 
//     //! Optionally refine mesh
//     void refine();
// 
//     //! Finish initiel mesh refinement
//     void finishref();
// 
//     //! Do uniform mesh refinement
//     void uniformRefine();
// 
//     //! Do error-based mesh refinement
//     void errorRefine();
// 
//     //! Do mesh refinement based on user explicitly tagging edges
//     void userRefine();
// 
//     //! Do mesh refinement correcting PE-boundary edges
//     void correctRefine( const tk::UnsMesh::EdgeSet& extra );
// 
//     //! Update mesh after refinement
//     void updateMesh();
// 
//     //! Update volume mesh after mesh refinement
//     void updateVolumeMesh( const std::unordered_set< std::size_t >& old,
//                            const std::unordered_set< std::size_t >& ref );
// 
//     //! Update boundary data structures after mesh refinement
//     void updateBoundaryMesh( const std::unordered_set< std::size_t >& old,
//                              const std::unordered_set< std::size_t >& ref );
// 
//     //! Evaluate initial conditions (IC) at mesh nodes
//     tk::Fields nodeinit( std::size_t npoin,
//                          const std::pair< std::vector< std::size_t >,
//                                           std::vector< std::size_t > >& esup );
// 
};

} // inciter::

#endif // Refiner_h
