// *****************************************************************************
/*!
  \file      src/Inciter/Discretization.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
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

namespace tk {
  class ExodusIIMeshWriter;
  class RootMeshWriter;
}

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
        const CProxy_BoundaryConditions& bc,
        const std::vector< std::size_t >& conn,
        const std::unordered_map< int,
                std::unordered_set< std::size_t > >& msum,
        const std::unordered_map< std::size_t, std::size_t >& filenodes,
        const tk::UnsMesh::EdgeNodes& edgenodes,
        int nchare );

    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wundefined-func-template"
    #endif
    //! Migrate constructor
    explicit Discretization( CkMigrateMessage* ) {}
    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

    //! Configure Charm++ reduction types
    static void registerReducers();

    //! \brief Read mesh node coordinates and optionally add new edge-nodes in
    //!   case of initial uniform refinement
    void coord();

    //! Collect nodal volumes across chare boundaries
    void comvol( const std::vector< std::size_t >& gid,
                 const std::vector< tk::real >& vol );

    //! Sum mesh volumes and contribute own mesh volume to total volume
    void totalvol();

    //! Compute mesh cell statistics
    void stat();

    /** @name Accessors
      * */
    ///@{
    const tk::UnsMesh::Coords& Coord() const{ return m_coord; }
    tk::UnsMesh::Coords& Coord() { return m_coord; }

    const std::vector< std::size_t >& Gid() const { return m_gid; }
    std::vector< std::size_t >& Gid() { return m_gid; }

    const std::unordered_map< std::size_t, std::size_t >& Lid() const
    { return m_lid; }
    std::unordered_map< std::size_t, std::size_t >& Lid() { return m_lid; }

    const std::vector< std::size_t >& Inpoel() const { return m_inpoel; }
    std::vector< std::size_t >& Inpoel() { return m_inpoel; }

    const std::vector< tk::real >& V() const { return m_v; }
    std::vector< tk::real >& V() { return m_v; }

    const std::vector< tk::real >& Vol() const { return m_vol; }
    std::vector< tk::real >& Vol() { return m_vol; }

    tk::real Dt() const { return m_dt; }
    tk::real T() const { return m_t; }
    uint64_t It() const { return m_it; }

    const tk::Timer& Timer() const { return m_timer; }
    tk::Timer& Timer() { return m_timer; }

    tk::real LastFieldWriteTime() const { return m_lastFieldWriteTime; }
    tk::real& LastFieldWriteTime() { return m_lastFieldWriteTime; }

    const CProxy_Transporter& Tr() const { return m_transporter; }
    CProxy_Transporter& Tr() { return m_transporter; }

    //! Access boundary conditions group local branch pointer
    BoundaryConditions* BC() { return m_bc.ckLocalBranch(); }

    //! Access bound DistFCT class pointer
    DistFCT* FCT() const {
      Assert(m_fct[ thisIndex ].ckLocal() != nullptr, "DistFCT ckLocal() null");
      return m_fct[ thisIndex ].ckLocal();
    }

    const std::unordered_map< std::size_t, std::size_t >& Filenodes() const
    { return m_filenodes; }
    std::unordered_map< std::size_t, std::size_t >& Filenodes()
    { return m_filenodes; }

    const std::unordered_map< std::size_t, std::size_t >& Bid() const
    { return m_bid; }
    std::unordered_map< std::size_t, std::size_t >& Bid() { return m_bid; }

    const std::unordered_map< int, std::vector< std::size_t > >& Msum() const
    { return m_msum; }
    std::unordered_map< int, std::vector< std::size_t > >& Msum()
    { return m_msum; }

    const std::string& OutFilename() const { return m_outFilename; }
    std::string& OutFilename() { return m_outFilename; }

    const std::pair< std::vector< std::size_t >, std::vector< std::size_t > >&
    Psup() const { return m_psup; }
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > >&
    Psup() { return m_psup; }
    //@}

    //! Output chare element blocks to output file
    void writeMesh();

    //! Output mesh-based fields metadata to file
    void writeNodeMeta() const;

    //! Output mesh-based element fields metadata to file
    void writeElemMeta() const;

    //! Output solution to file
    void writeNodeSolution( const tk::ExodusIIMeshWriter& ew,
                        uint64_t it,
                        const std::vector< std::vector< tk::real > >& u ) const;
    #ifdef HAS_ROOT
    void writeNodeSolution( const tk::RootMeshWriter& rmw,
                        uint64_t it,
                        const std::vector< std::vector< tk::real > >& u ) const;
    #endif
    void writeElemSolution( const tk::ExodusIIMeshWriter& ew,
                            uint64_t it,
                            const std::vector< std::vector< tk::real > >& u )
                          const;

    //! Set time step size
    void setdt( tk::real newdt );

    //! Prepare for next step
    void next();

    //! Otput one-liner status report
    void status();

    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_Discretization::pup(p);
      p | m_it;
      p | m_t;
      p | m_dt;
      p | m_lastFieldWriteTime;
      p | m_nvol;
      p | m_outFilename;
      p | m_fct;
      p | m_transporter;
      p | m_bc;
      p | m_filenodes;
      p | m_edgenodes;
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
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Discretization object reference
    friend void operator|( PUP::er& p, Discretization& i ) { i.pup(p); }
    //@}

  private:
    //! Iteration count
    uint64_t m_it;
     //! Physical time
    tk::real m_t;
    //! Physical time step size
    tk::real m_dt;
    //! Physical time at which the last field output was written
    tk::real m_lastFieldWriteTime;
    //! \brief Number of chares from which we received nodal volume
    //!   contributions on chare boundaries
    std::size_t m_nvol;
    //! Output filename
    std::string m_outFilename;
    //! Distributed FCT proxy
    CProxy_DistFCT m_fct;
    //! Transporter proxy
    CProxy_Transporter m_transporter;
    //! Boundary conditions proxy
    CProxy_BoundaryConditions m_bc;
    //! Map associating file node IDs to local node IDs
    std::unordered_map< std::size_t, std::size_t > m_filenodes;
    //! \brief Maps associating node node IDs to edges (a pair of old node IDs)
    //!   for only the nodes newly added as a result of initial uniform
    //!   refinement.
    tk::UnsMesh::EdgeNodes m_edgenodes;
    //! \brief Elements of the mesh chunk we operate on
    //! \details Initialized by the constructor. The first vector is the element
    //!   connectivity (local IDs), the second vector is the global node IDs of
    //!   owned elements, while the third one is a map of global->local node
    //!   IDs.
    std::tuple< std::vector< std::size_t >,
                std::vector< std::size_t >,
                std::unordered_map< std::size_t, std::size_t > > m_el;
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

    //! Sum mesh volumes to nodes, start communicating them on chare-boundaries
    void vol();

    //! Read coordinates of mesh nodes given
    void readCoords();

    //! \brief Add coordinates of mesh nodes newly generated to edge-mid points
    //!    during initial refinement
    void addEdgeNodeCoords();
};

} // inciter::

#endif // Discretization_h
