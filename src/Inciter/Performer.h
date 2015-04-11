//******************************************************************************
/*!
  \file      src/Inciter/Performer.h
  \author    J. Bakosi
  \date      Fri 10 Apr 2015 03:58:33 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Performer advances the Euler equations
  \details   Performer advances the Euler equations. There are a potentially
    large number of Performer Charm++ chares created by Conductor. Each
    performer gets a chunk of the full load (part of the mesh) and does the
    same: initializes and advances the Euler equations in time.
*/
//******************************************************************************
#ifndef Performer_h
#define Performer_h

#include <iostream>     // NOT REALLY NEEDED

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <inciter.decl.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#include <DerivedData.h>
#include <Inciter/InputDeck/InputDeck.h>
#include <ExodusIIMeshReader.h>
#include <ExodusIIMeshWriter.h>
#include <performer.decl.h>

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern
  std::pair< std::vector< std::size_t >, std::vector< std::size_t > > g_esup;
extern std::vector< int > g_tetinpoel;
extern std::vector< std::size_t > g_colors;
extern
  std::vector< std::map< std::size_t, std::vector< std::size_t > > > g_comm;

//! Performer Charm++ chare used to advance the Euler equations in time
template< class Proxy >
class Performer : public CBase_Performer< Proxy > {

  private:
    using Array = CBase_Performer< Proxy >;

  public:
    //! Constructor
    //! \param[in] proxy Host (Conductor) proxy to call back to
    explicit Performer( Proxy& proxy ) :
      m_proxy( proxy ),
      m_export( g_comm.size() > Array::thisIndex ?
                g_comm[ static_cast<std::size_t>(Array::thisIndex) ] :
                std::map< std::size_t, std::vector< std::size_t > >() )
    {
      // Initialize communication maps
      initImports();
      // Initialize global->local and local->global node ids, global element ids
      std::map< std::size_t, std::size_t > lnode;
      std::vector< std::size_t > gnode;
      std::set< std::size_t > gelem;
      std::tie( lnode, gnode, gelem ) = initIds();
      // Initialize element connectivity
      initInpoel( lnode, gelem );
      // Generate derived data structures
      initDerivedData();
      // Read coordinates of owned and received mesh nodes
      initCoords( gnode );
      // Output chare mesh and nodal chare id field to file
      writeChareId( gnode.size() );
      // Signal host that initialization is complete
      m_proxy.init();
    }

    //! Migrate constructor
    Performer( CkMigrateMessage* ) {}

  private:
    Proxy m_proxy;                      //!< Host proxy

    //! Communication maps
    std::map< std::size_t, std::vector< std::size_t > > m_export;
    std::map< std::size_t, std::vector< std::size_t > > m_import;

    //! Mesh node coordinates
    std::vector< tk::real > m_x;
    std::vector< tk::real > m_y;
    std::vector< tk::real > m_z;

    //! Tetrahedron connectivity with local point ids
    std::vector< int > m_tetinpoel;

    //! Data structures derived from element connectivity
    std::pair< std::vector< std::size_t >, std::vector< std::size_t > > m_edsup;

    //! Initialize import map
    void initImports() {
      std::size_t h = 0;
      for (const auto& m : g_comm) {
        for (const auto& x : m)
          if (Array::thisIndex == x.first)
            for (auto p : x.second)
              m_import[ h ].push_back( p );
        ++h;
      }
    }

    //! \brief Initialize global->local and local->global node ids and global
    //!   element ids
    //! \return A tuple containing at 0 vector os local node ids, at 1 vector of
    //!   global node ids, at 2 set of unique global element ids
    //! \author J. Bakosi
    std::tuple< std::map< std::size_t, std::size_t >,
                std::vector< std::size_t >,
                std::set< std::size_t > >
    initIds() {
      // Generate global->local and local->global node ids for owned nodes
      std::vector< std::size_t > gnode;
      std::map< std::size_t, std::size_t > lnode;
      std::size_t gid = 0, lid = 0;
      for (auto n : g_colors) {
        if (n == Array::thisIndex) {
          gnode.push_back( gid );
          lnode[ gid ] = lid;
          ++lid;
        }
        ++gid;
      }
      // Create set of unique owned global element ids
      std::set< std::size_t > gelem;
      for (auto p : gnode)      // at this stage: only consider owned nodes
        for (auto i=g_esup.second[p]+1; i<=g_esup.second[p+1]; ++i)
          gelem.insert( g_esup.first[i] );
      // Add received node ids to global->local and local->global node ids
      for (auto& m : m_import)
        for (auto p : m.second) {
          gnode.push_back( p );
          lnode[ p ] = lid;
          ++lid;
        }
      Assert( m_import.empty() ? lnode.size() == gnode.size() : true,
              "The size of global->local and local->global node id map must "
              "be equal if there is only a single chare" );
      return std::make_tuple( lnode, gnode, gelem );
    }

    //! Initialize element connectivity
    //! \param[in] lnode Global->local node id map
    //! \param[in] gelem Set of unique owned global element ids
    //! \author J. Bakosi
    void initInpoel( const std::map< std::size_t, std::size_t >& lnode,
                     const std::set< std::size_t >& gelem ) {
      // Lambda for finding local node id for global node id
      auto localId = [ &lnode ]( int gid ) -> int {
        const auto it = lnode.find( static_cast< std::size_t >( gid ) );
        if (it != lnode.end()) return static_cast< int >( it->second );
        else Throw( "Can't find global node id " + std::to_string(gid) );
      };
      // Extract tetrahedron element connectivity of those elements that
      // contain at least a single owned node id and store using local node ids
      for (auto e : gelem) {
        m_tetinpoel.push_back( localId( g_tetinpoel[e*4] ) );
        m_tetinpoel.push_back( localId( g_tetinpoel[e*4+1] ) );
        m_tetinpoel.push_back( localId( g_tetinpoel[e*4+2] ) );
        m_tetinpoel.push_back( localId( g_tetinpoel[e*4+3] ) );
      }
    }

    //! Initialize data structures derived from mesh connectivity
    void initDerivedData() {
      // Generate edges surrounding points
      m_edsup = tk::genEdsup( m_tetinpoel, 4, tk::genEsup(m_tetinpoel,4) );
    }

    //! Read coordinates of owned and received mesh nodes
    void initCoords( const std::vector< std::size_t >& gnode ) {
      tk::UnsMesh inmesh;
      tk::ExodusIIMeshReader
        er( g_inputdeck.get< tag::cmd, tag::io, tag::input >(), inmesh );
      for (auto p : gnode) er.readNode( p, m_x, m_y, m_z );
    }

    //! Output chare mesh and nodal chare id field to file
    //! \param[in] npoin Number points in chare mesh
    void writeChareId( std::size_t npoin ) {
      // Create mesh object initializing element connectivity and point coords
      tk::UnsMesh outmesh( m_tetinpoel, m_x, m_y, m_z );
      // Create ExodusII writer
      tk::ExodusIIMeshWriter
        ew( g_inputdeck.get< tag::cmd, tag::io, tag::output >() +
              "_" + std::to_string(Array::thisIndex),
            outmesh );
      // Write chare mesh
      ew.write();
      // Write chare id field to mesh
      std::vector< tk::real >
        chareid( npoin, static_cast< tk::real >( Array::thisIndex ) );
      ew.writeVarNames( { "Chare Id" } );
      ew.writeTimeStamp( 1, 1.0 );
      ew.writeNodeScalar( 1, 1, chareid );
    }
};

} // inciter::

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#define CK_TEMPLATES_ONLY
#include <performer.def.h>
#undef CK_TEMPLATES_ONLY

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

#endif // Performer_h
