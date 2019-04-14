// *****************************************************************************
/*!
  \file      src/Inciter/Discretization.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \details   Data and functionality common to all discretization schemes
     The Discretization class contains data and functionality common to all
     discretization schemes.
*/
// *****************************************************************************
#ifndef Discretization_h
#define Discretization_h

#include "Types.h"
#include "Timer.h"
#include "Keywords.h"
#include "Fields.h"
#include "PUPUtil.h"
#include "PDFReducer.h"
#include "UnsMesh.h"

#include "NoWarning/discretization.decl.h"
#include "NoWarning/refiner.decl.h"

namespace inciter {

//! \brief Discretization Charm++ chare array holding common functinoality to
//!   all discretization schemes
class Discretization : public CBase_Discretization {

  public:
    //! Constructor
    explicit
      Discretization(
        const CProxy_DistFCT& fctproxy,
        const CProxy_Transporter& transporter,
        const tk::CProxy_MeshWriter& meshwriter,
        const std::vector< std::size_t >& ginpoel,
        const tk::UnsMesh::CoordMap& coordmap,
        const std::map< int, std::unordered_set< std::size_t > >& msum,
        int nc );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    // cppcheck-suppress uninitMemberVar
    explicit Discretization( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! Resize mesh data structures (e.g., after mesh refinement)
    void resize( const tk::UnsMesh::Chunk& chunk,
                 const tk::UnsMesh::Coords& coord,
                 const std::unordered_map< int,
                         std::vector< std::size_t > >& msum );

    //! Sum mesh volumes to nodes, start communicating them on chare-boundaries
    void vol();

    //! Set Refiner Charm++ proxy
    void setRefiner( const CProxy_Refiner& ref );

    //! Collect nodal volumes across chare boundaries
    void comvol( const std::vector< std::size_t >& gid,
                 const std::vector< tk::real >& nodevol );

    //! Sum mesh volumes and contribute own mesh volume to total volume
    void totalvol();

    //! Compute mesh cell statistics
    void stat();

    /** @name Accessors */
    ///@{
    //! Coordinates accessors as const-ref
    const tk::UnsMesh::Coords& Coord() const{ return m_coord; }
    //! Coordinates accessors as non-const ref
    tk::UnsMesh::Coords& Coord() { return m_coord; }

    //! Global ids accessors as const-ref
    const std::vector< std::size_t >& Gid() const { return m_gid; }
    //! Global ids accessors as non-const-ref
    std::vector< std::size_t >& Gid() { return m_gid; }

    //! Local ids accessors as const-ref
    const std::unordered_map< std::size_t, std::size_t >& Lid() const
    { return m_lid; }
    //! Local ids accessors as non-const-ref
    std::unordered_map< std::size_t, std::size_t >& Lid() { return m_lid; }

    //! Tetrahedron element connectivity (with local ids) accessors as const-ref
    const std::vector< std::size_t >& Inpoel() const { return m_inpoel; }
    //! \brief Tetrahedron element connectivity (with local ids) accessors as
    //!    non-const-ref
    std::vector< std::size_t >& Inpoel() { return m_inpoel; }

    //! Total mesh volume accessors const-ref
    const std::vector< tk::real >& V() const { return m_v; }
    //! Total mesh volume accessors non-const-ref
    std::vector< tk::real >& V() { return m_v; }

    //! Nodal mesh volumes accessors as const-ref
    const std::vector< tk::real >& Vol() const { return m_vol; }
    //! Nodal mesh volumes accessors as non-const-ref
    std::vector< tk::real >& Vol() { return m_vol; }

    //! Time step size accessor
    tk::real Dt() const { return m_dt; }
    //! Physical time accessor
    tk::real T() const { return m_t; }
    //! Iteration count accessor
    uint64_t It() const { return m_it; }

    //! Non-const-ref refinement iteration count accessor
    uint64_t& Itr() { return m_itr; }
    //! Non-const-ref field-output iteration count accessor
    uint64_t& Itf() { return m_itf; }

    //! Timer accessor as const-ref
    const tk::Timer& Timer() const { return m_timer; }
    //! Timer accessor as non-const-ref
    tk::Timer& Timer() { return m_timer; }

    //! Accessor to flag indicating if the mesh was refined as a value
    int refined() const { return m_refined; }
    //! Accessor to flag indicating if the mesh was refined as non-const-ref
    int& refined() { return m_refined; }

    //! Transporter proxy accessor as const-ref
    const CProxy_Transporter& Tr() const { return m_transporter; }
    //! Transporter proxy accessor as non-const-ref
    CProxy_Transporter& Tr() { return m_transporter; }

    //! Access bound Refiner class pointer
    Refiner* Ref() const {
      Assert( m_refiner[ thisIndex ].ckLocal() != nullptr,
              "Refiner ckLocal() null" );
      return m_refiner[ thisIndex ].ckLocal();
    }

    //! Access bound DistFCT class pointer
    DistFCT* FCT() const {
      Assert(m_fct[ thisIndex ].ckLocal() != nullptr, "DistFCT ckLocal() null");
      return m_fct[ thisIndex ].ckLocal();
    }

    //! Boundary node ids accessor as const-ref
    const std::unordered_map< std::size_t, std::size_t >& Bid() const
    { return m_bid; }
    //! Boundary node ids accessor as non-const-ref
    std::unordered_map< std::size_t, std::size_t >& Bid() { return m_bid; }

    //! Nodal communication map accessor as const-ref
    const std::unordered_map< int, std::vector< std::size_t > >& Msum() const
    { return m_msum; }
    //! Nodal communication map accessor as non-const-ref
    std::unordered_map< int, std::vector< std::size_t > >& Msum()
    { return m_msum; }

    //! Points surrounding points accessor as const-ref
    const std::pair< std::vector< std::size_t >, std::vector< std::size_t > >&
    Psup() const { return m_psup; }
    //! Points surrounding points accessor as non-const-ref
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >&
    Psup() { return m_psup; }
    //@}

    //! Set time step size
    void setdt( tk::real newdt );

    //! Prepare for next step
    void next();

    //! Otput one-liner status report
    void status();

    //! Output mesh and fields data (solution dump) to file(s)
    void write( const std::vector< std::size_t >& inpoel,
                const tk::UnsMesh::Coords& coord,
                const std::map< int, std::vector< std::size_t > >& bface,
                const std::map< int, std::vector< std::size_t > >& bnode,
                const std::vector< std::size_t >& triinpoel,
                const std::vector< std::string>& elemfieldnames,
                const std::vector< std::string>& nodefieldnames,
                const std::vector< std::vector< tk::real > >& elemfields,
                const std::vector< std::vector< tk::real > >& nodefields,
                CkCallback c );

    //! Return chare-node adjacency map as sets
    std::unordered_map< int, std::unordered_set< std::size_t > >
    msumset() const;

    /** @name Charm++ pack/unpack serializer member functions */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) override {
      p | m_nchare;
      p | m_it;
      p | m_itr;
      p | m_itf;
      p | m_initial;
      p | m_t;
      p | m_lastDumpTime;
      p | m_dt;
      p | m_nvol;
      p | m_fct;
      p | m_transporter;
      p | m_meshwriter;
      p | m_refiner;
      p | m_el;
      if (p.isUnpacking()) {
        m_inpoel = std::get< 0 >( m_el );
        m_gid = std::get< 1 >( m_el );
        m_lid = std::get< 2 >( m_el );
      }
      p | m_coord;
      p | m_psup;
      p | m_msum;
      p | m_v;
      p | m_vol;
      p | m_volc;
      p | m_bid;
      p | m_timer;
      p | m_refined;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Discretization object reference
    friend void operator|( PUP::er& p, Discretization& i ) { i.pup(p); }
    //@}

  private:
    //! Total number of Discretization chares
    int m_nchare;
    //! Iteration count
    uint64_t m_it;
    //! Iteration count with mesh refinement
    //! \details Used as the restart sequence number {RS} in saving output in
    //!    an ExodusII sequence
    //! \see https://www.paraview.org/Wiki/Restarted_Simulation_Readers
    uint64_t m_itr;
    //! Field output iteration count without mesh refinement
    //! \details Counts the number of field outputs to file during two
    //!   time steps with mesh efinement
    uint64_t m_itf;
    //! Flag that is nonzero during setup and zero during time stepping
    tk::real m_initial;
    //! Physical time
    tk::real m_t;
    //! Physical time at last field output
    tk::real m_lastDumpTime;
    //! Physical time step size
    tk::real m_dt;
    //! \brief Number of chares from which we received nodal volume
    //!   contributions on chare boundaries
    std::size_t m_nvol;
    //! Distributed FCT proxy
    CProxy_DistFCT m_fct;
    //! Transporter proxy
    CProxy_Transporter m_transporter;
    //! Mesh writer proxy
    tk::CProxy_MeshWriter m_meshwriter;
    //! Mesh refiner proxy
    CProxy_Refiner m_refiner;
    //! \brief Elements of the mesh chunk we operate on
    //! \details Initialized by the constructor. The first vector is the element
    //!   connectivity (local IDs), the second vector is the global node IDs of
    //!   owned elements, while the third one is a map of global->local node
    //!   IDs.
    tk::UnsMesh::Chunk m_el;
    //! Alias to element connectivity
    std::vector< std::size_t >& m_inpoel = std::get< 0 >( m_el );
    //! Alias to global node IDs of owned elements
    std::vector< std::size_t >& m_gid = std::get< 1 >( m_el );
    //! \brief Alias to local node ids associated to the global ones of owned
    //!    elements
    std::unordered_map< std::size_t, std::size_t >& m_lid = std::get< 2 >( m_el );
    //! Mesh point coordinates
    tk::UnsMesh::Coords m_coord;
    //! Points surrounding points of our chunk of the mesh
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_psup;
    //! \brief Global mesh node IDs bordering the mesh chunk held by fellow
    //!   Discretization chares associated to their chare IDs
    //! \details msum: mesh chunks surrounding mesh chunks and their neighbor
    //!   points
    std::unordered_map< int, std::vector< std::size_t > > m_msum;
    //! Nodal mesh volumes
    //! \details This is the volume of the mesh associated to nodes of owned
    //!   elements (sum of surrounding cell volumes / 4) without contributions
    //!   from other chares on chare-boundaries
    std::vector< tk::real > m_v;
    //! Volume of nodes
    //! \details This is the volume of the mesh associated to nodes of owned
    //!   elements (sum of surrounding cell volumes / 4) with contributions from
    //!   other chares on chare-boundaries
    std::vector< tk::real > m_vol;
    //! Receive buffer for volume of nodes
    //! \details This is a communication buffer used to compute the volume of
    //!   the mesh associated to nodes of owned elements (sum of surrounding
    //!   cell volumes / 4) with contributions from other chares on
    //!   chare-boundaries.
    std::vector< tk::real > m_volc;
    //! \brief Local chare-boundary mesh node IDs at which we receive
    //!   contributions associated to global mesh node IDs of mesh elements we
    //!   contribute to
    std::unordered_map< std::size_t, std::size_t > m_bid;
    //! Timer measuring a time step
    tk::Timer m_timer;
    //! 1 if mesh was refined in a time step, 0 if it was not
    int m_refined;

    //! Set mesh coordinates based on coordinates map
    tk::UnsMesh::Coords setCoord( const tk::UnsMesh::CoordMap& coordmap );
};

} // inciter::

#endif // Discretization_h
