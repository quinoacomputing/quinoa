// *****************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Mon 18 Jul 2016 10:56:24 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Performer advances a PDE
  \details   Performer advances a PDE. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances a PDE in time.
*/
// *****************************************************************************
#ifndef Performer_h
#define Performer_h

#include <array>
#include <cstddef>
#include <iosfwd>
#include <utility>
#include <vector>
#include <cstring>
#include <cmath>
#include <unordered_map>

#include "Types.h"
#include "MeshNodes.h"
#include "DerivedData.h"
#include "VectorReducer.h"
#include "Inciter/InputDeck/InputDeck.h"

#include "NoWarning/conductor.decl.h"
#include "NoWarning/performer.decl.h"

namespace tk { class ExodusIIMeshWriter; }

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern CkReduction::reducerType VerifyBCMerger;

//! Performer Charm++ chare used to advance a PDE in time
class Performer : public CBase_Performer {

  private:
    using ConductorProxy = CProxy_Conductor;
    using LinSysMergerProxy = tk::CProxy_LinSysMerger< CProxy_Conductor,
                                                       CProxy_Performer >;
  public:
    //! Constructor
    explicit
      Performer( const CProxy_Conductor& conductor,
                 const LinSysMergerProxy& lsm,
                 const std::vector< std::size_t >& conn,
                 const std::unordered_map< std::size_t, std::size_t >& cid );

    #if defined(__GNUC__)
      #pragma GCC diagnostic push
      #pragma GCC diagnostic ignored "-Weffc++"
    #endif

    //! Migrate constructor
    explicit Performer( CkMigrateMessage* ) {}

    #if defined(__GNUC__)
      #pragma GCC diagnostic pop
    #endif

    //! \brief Configure Charm++ reduction types for concatenating BC nodelists
    //! \details Since this is a [nodeinit] routine, see linsysmerger.ci, the
    //!   Charm++ runtime system executes the routine exactly once on every
    //!   logical node early on in the Charm++ init sequence. Must be static as
    //!   it is called without an object. See also: Section "Initializations at
    //!   Program Startup" at in the Charm++ manual
    //!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
    static void registerVerifyBCMerger()
    { VerifyBCMerger = CkReduction::addReducer( tk::mergeVector ); }

    //! Initialize mesh IDs, element connectivity, coordinates
    void setup();

    //! Request owned node IDs on which a Dirichlet BC is set by the user
    void requestBCs();

    //! Look up and return old node ID for new one
    void oldID( int frompe, const std::vector< std::size_t >& newids );

    //! Look up boundary condition values at node IDs for all PDEs
    void bcval( int frompe, const std::vector< std::size_t >& nodes );

    //! Initialize communication and mesh data
    void init( tk::real dt );

    //! Update solution vector
    void updateSolution( const std::vector< std::size_t >& gid,
                         const std::vector< tk::real >& sol );

    //! Advance equations to next stage in multi-stage time stepping
    void advance( uint8_t stage, tk::real dt, uint64_t it, tk::real t );


    /** @name Pack/Unpack: Serialize Performer object for Charm++ */
    ///@{
    //! \brief Pack/Unpack serialize member function
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    void pup( PUP::er &p ) {
      CBase_Performer::pup(p);
      p | m_it;
      p | m_itf;
      p | m_t;
      p | m_stage;
      p | m_nsol;
      p | m_outFilename;
      p | m_conductor;
      p | m_linsysmerger;
      p | m_cid;
      p | m_el;
      if (p.isUnpacking()) { m_inpoel = m_el.first; m_gid = m_el.second; }
      p | m_lid;
      p | m_coord;
      p | m_psup;
      p | m_u; p | m_uf; p | m_un;
      p | m_lhsd; p | m_lhso;
    }
    //! \brief Pack/Unpack serialize operator|
    //! \param[in,out] p Charm++'s PUP::er serializer object reference
    //! \param[in,out] i Performer object reference
    //! \author J. Bakosi
    friend void operator|( PUP::er& p, Performer& i ) { i.pup(p); }

  private:
    using ncomp_t = kw::ncomp::info::expect::type;

    uint64_t m_it;                      //!< Iteration count
    uint64_t m_itf;                     //!< Field output iteration count
    tk::real m_t;                       //!< Physical time
    uint8_t m_stage;                    //!< Stage in multi-stage time stepping
    std::size_t m_nsol;                 //!< Counter for solution nodes updated
    std::string m_outFilename;          //!< Output filename
    ConductorProxy m_conductor;         //!< Conductor proxy
    LinSysMergerProxy m_linsysmerger;   //!< Linear system merger proxy
    //! \brief Map associating old node IDs (as in file) to new node IDs (as in
    //!   producing contiguous-row-id linear system contributions)
    std::unordered_map< std::size_t, std::size_t > m_cid;
    //! \brief Elements of the mesh chunk we operate on
    //! \details Initialized by the constructor. The first vector is the element
    //!   connectivity (local IDs), while the second vector is the global node
    //!   IDs of owned elements.
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_el;
    //! Alias to element connectivity in m_el
    decltype(m_el.first)& m_inpoel = m_el.first;
    //! Alias to global node IDs of owned elements in m_el
    decltype(m_el.second)& m_gid = m_el.second;
    //!< Local node ids associated to the global ones of owned elements
    std::unordered_map< std::size_t, std::size_t > m_lid;
    //! Mesh point coordinates
    std::array< std::vector< tk::real >, 3 > m_coord;
    //! Points surrounding points of our chunk of the mesh
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_psup;
    //! Unknown/solution vector: global mesh point row ids and values
    tk::MeshNodes m_u, m_uf, m_un;
    //! Sparse matrix sotring the diagonals and off-diagonals of nonzeros
    tk::MeshNodes m_lhsd, m_lhso;

    //! Send off global row IDs to linear system merger, setup global->local IDs
    void setupIds();

    //! Extract node IDs from element side sets and match to BCs
    std::vector< std::size_t > queryBCs();

    //! Query old node IDs for a list of new node IDs
    std::vector< std::size_t > old( const std::vector< std::size_t >& newids );

    //! Send node list to our LinSysMerger branch which is then used to set BCs
    void sendBCs( const std::vector< std::size_t >& bc );

    //! Read coordinates of mesh nodes given
    void readCoords();

    //! Compute left-hand side matrix of PDE
    void lhs();

    //! Compute righ-hand side vector of PDE
    void rhs( tk::real mult,
              tk::real dt,
              const tk::MeshNodes& sol,
              tk::MeshNodes& rhs );

    //! Output chare mesh to file
    void writeMesh();

    //! Output chare mesh chare id field to file
    void writeChareId( const tk::ExodusIIMeshWriter& ew, uint64_t it ) const;

    //! Output solution to file
    void writeSolution( const tk::ExodusIIMeshWriter& ew,
                        uint64_t it,
                        const std::vector< std::vector< tk::real > >& u ) const;

    //! Output mesh-based fields metadata to file
    void writeMeta() const;

    //! Output mesh-based fields to file
    void writeFields( tk::real time );
};

} // inciter::

#endif // Performer_h
