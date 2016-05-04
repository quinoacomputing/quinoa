//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Wed 04 May 2016 10:44:49 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Performer advances a PDE
  \details   Performer advances a PDE. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances a PDE in time.
*/
//******************************************************************************
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
#include "Inciter/InputDeck/InputDeck.h"

#include "NoWarning/conductor.decl.h"
#include "NoWarning/performer.decl.h"

namespace tk { class ExodusIIMeshWriter; }

namespace inciter {

extern ctr::InputDeck g_inputdeck;

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

    //! Migrate constructor
    explicit Performer( CkMigrateMessage* ) :
      // WARNING: This is a "blind" copy of the standard constructor initializer
      // list - it must be changed for migration to be correct.
      m_it( 0 ),
      m_itf( 0 ),
      m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
      m_stage( 0 ),
      m_nsol( 0 ),
      m_outFilename( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + "." +
                     std::to_string( thisIndex ) ),
      m_conductor(),
      m_linsysmerger(),
      m_cid(),
      m_el(),     // fills m_inpoel and m_gid
      m_lid(),
      m_coord(),
      m_psup( tk::genPsup( m_inpoel, 4, tk::genEsup(m_inpoel,4) ) ),
      m_u( 0, g_inputdeck.get< tag::component >().nprop() ),
      m_uf( 0, g_inputdeck.get< tag::component >().nprop() ),
      m_un( 0, g_inputdeck.get< tag::component >().nprop() ),
      m_lhsd( 0, g_inputdeck.get< tag::component >().nprop() ),
      m_lhso( 0, g_inputdeck.get< tag::component >().nprop() )
    {}

    //! Initialize mesh IDs, element connectivity, coordinates
    void setup();

    //! Initialize communication and mesh data
    void init( tk::real dt );

    //! Update solution vector
    void updateSolution( const std::vector< std::size_t >& gid,
                         const std::vector< tk::real >& sol );

    //! Advance equations to next stage in multi-stage time stepping
    void advance( uint8_t stage, tk::real dt, uint64_t it, tk::real t );

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
    //! Alias to global node IDs of owned elements in in m_el
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
