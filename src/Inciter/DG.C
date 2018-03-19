// *****************************************************************************
/*!
  \file      src/Inciter/DG.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     DG advances a system of PDEs with the discontinuous Galerkin scheme
  \details   DG advances a system of partial differential equations (PDEs) using
    discontinuous Galerkin (DG) finite element (FE) spatial discretization (on
    tetrahedron elements) combined with Runge-Kutta (RK) time stepping.
  \see The documentation in DG.h.
*/
// *****************************************************************************

#include <algorithm>

#include "DG.h"
#include "Discretization.h"
#include "DGPDE.h"
#include "Solver.h"
#include "DiagReducer.h"
#include "DerivedData.h"
#include "Diagnostics.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "ExodusIIMeshWriter.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< DGPDE > g_dgpde;

static CkReduction::reducerType DiagMerger;

} // inciter::

using inciter::DG;

DG::DG( const CProxy_Discretization& disc,
        const tk::CProxy_Solver& solver,
        const FaceData& fd ) :
  m_solver( solver ),
  m_nfac( 0 ),
  m_nadj( 0 ),
  m_itf( 0 ),
  m_disc( disc ),
  m_fd( fd ),
  m_u( m_disc[thisIndex].ckLocal()->Inpoel().size()/4,
       g_inputdeck.get< tag::component >().nprop() ),
  m_un( m_u.nunk(), m_u.nprop() ),
  m_vol( 0.0 ),
  m_geoFace( tk::genGeoFaceTri( fd.Ntfac(), fd.Inpofa(),
                                m_disc[thisIndex].ckLocal()->Coord()) ),
  m_geoElem( tk::genGeoElemTet( m_disc[thisIndex].ckLocal()->Inpoel(),
                                m_disc[thisIndex].ckLocal()->Coord() ) ),
  m_lhs( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_facecnt( fd.Inpofa().size()/3 ),
  m_msumset( msumset() ),
  m_esuelTet( tk::genEsuelTet( m_disc[thisIndex].ckLocal()->Inpoel(),
                tk::genEsup( m_disc[thisIndex].ckLocal()->Inpoel(), 4 ) ) ),
  m_ipface(),
  m_potBndFace(),
  m_bndFace(),
  m_ghostData(),
  m_ghostReq(),
  m_ghost()
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
  // Activate SDAG waits for face adjacency map (ghost data) calculation
  wait4ghost();

  auto d = Disc();

  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();
  const auto& inpofa = fd.Inpofa();

  m_esuf = fd.Esuf();

  // Invert inpofa to enable searching for faces based on (global) node triplets
  Assert( inpofa.size() % 3 == 0, "Inpofa must contain triplets" );
  for (std::size_t f=0; f<inpofa.size()/3; ++f)
    m_ipface.insert( {{{ gid[ inpofa[f*3+0] ],
                         gid[ inpofa[f*3+1] ],
                         gid[ inpofa[f*3+2] ] }}} );

  // At this point m_ipface has node-id-triplets (faces) on the internal
  // chare-domain and on the physical boundary but not on chare boundaries,
  // hence the name internal + physical boundary faces.

  // Lambda to find the boundary chare for a face (given by 3 global node IDs)
  auto facechare = [&]( const tk::UnsMesh::Face& t ) -> int {
    for (const auto& n : m_msumset) {  // for all neighbor chares
      if ( n.second.find(t[0]) != end(n.second) &&
           n.second.find(t[1]) != end(n.second) &&
           n.second.find(t[2]) != end(n.second) )
       return n.first;
    }
    Throw( "Face not found on chare-node communication map" );
  };

  // Build a set of faces (each face given by 3 global node IDs) associated to
  // chares we potentially share boundary faces with
  for (std::size_t e=0; e<m_esuelTet.size()/4; ++e) {   // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f)     // for all tet faces
      if (m_esuelTet[mark+f] == -1) {   // if face has no outside-neighbor tet
        tk::UnsMesh::Face t {{ gid[ inpoel[ mark + tk::lpofa[f][0] ] ],
                               gid[ inpoel[ mark + tk::lpofa[f][1] ] ],
                               gid[ inpoel[ mark + tk::lpofa[f][2] ] ] }};
        // if does not exist in ipface, store as a potential chare-boundary face
        // associated to neighbor chare
        if (m_ipface.find(t) == end(m_ipface))
          m_potBndFace[ facechare(t) ].insert( t );
      }
  }

  // Note that while the (potential) boundary-face adjacency map (m_potBndFace)
  // is derived from the node adjacency map (m_msumset), their sizes is not
  // necessarily the same. This is because while a node can be shared at a
  // single corner or along an edge, that does not necessarily share a face as
  // well (in other words, shared nodes or edges can exist that are not part of
  // a shared face). So the chares we communicate with across faces are not
  // necessarily the same as the chares we would communicate nodes with.
  //
  // Since the sizes of the node and face adjacency maps are not the same, while
  // sending the faces on chare boundaries would be okay, however, the receiver
  // would not necessarily know how many chares it must receive from. To solve
  // this problem we send to chares which we share at least a single node with,
  // i.e., rely on the node-adjacency map, but we either send a set of faces
  // which do share faces on the chare boundary or an empty set if there is not
  // a single face that shares a face with the destination chare (only single
  // nodes or edges). The assumption here is, of course, that the size of the
  // face adjacency map is always smaller than or equal to that of the node
  // adjacency map, which is always true. Since the receive side already knows
  // how many fellow chares it must receive shared node ids from, we use that to
  // detect completion of the number of receives in comfac(). This simplifies
  // the communication pattern and code for a small price of sending a few empty
  // sets (for those chare boundaries that only share individual nodes but not
  // faces).

  // Send sets of faces adjacent to chare boundaries to fellow workers (if any)
  if (m_msumset.empty())        // in serial, skip setting up ghosts altogether
    adj();
  else
    for (const auto& c : m_msumset) {   // for all chares we share nodes with
      decltype(m_potBndFace)::mapped_type fs;      // will store face set
      auto f = m_potBndFace.find( c.first );       // find face set for chare
      if (f != end(m_potBndFace)) fs = std::move( f->second ); // if found, send
      thisProxy[ c.first ].comfac( thisIndex, fs );
    }
}

std::unordered_map< int, std::unordered_set< std::size_t > >
DG::msumset() const
// *****************************************************************************
// Convert chare-node adjacency map to hold sets instead of vectors
//! \return Chare-node adjacency map that holds sets instead of vectors
// *****************************************************************************
{
  auto d = Disc();

  std::unordered_map< int, std::unordered_set< std::size_t > > m;
  for (const auto& n : d->Msum())
    m[ n.first ].insert( n.second.cbegin(), n.second.cend() );

  return m;
}

void
DG::comfac( int fromch, const tk::UnsMesh::FaceSet& infaces )
// *****************************************************************************
// Receive unique set of faces we potentially share with/from another chare
// *****************************************************************************
{
  // Attempt to find sender chare among chares we potentially share faces with.
  // Note that it is feasible that a sender chare called us but we do not have a
  // set of faces associated to that chare. This can happen if we only share a
  // single node or an edge but note a face with that chare.
  auto b = m_potBndFace.find( fromch );
  if (b != end(m_potBndFace)) {
    auto& bndface = m_bndFace[ fromch ];  // will associate to sender chare
    // try to find incoming faces among our faces we potentially share with
    // fromch, if found, generate and assign new local ID to face (associated to
    // sender chare)
    for (const auto& t : infaces) {
      if (b->second.find(t) != end(b->second))
        bndface[ t ] = m_facecnt++;
    }
    // if at this point we have not found any face among our faces we
    // potentially share with fromch, there is no need to keep an empty set of
    // faces associated to fromch as we only share nodes or edges with it, but
    // not faces
    if (bndface.empty()) m_bndFace.erase( fromch );
  }

  if (++m_nfac == m_msumset.size()) {
    // At this point m_bndFace is complete on this PE. This means that
    // starting from the sets of faces we potentially share with fellow chares
    // (m_potBndFace), we now only have those faces we actually share faces with
    // (through which we need to communicate later). Also, m_bndFace not only
    // has the unique faces associated to fellow chares, but also a newly
    // assigned local face ID. We continue by starting setting up ghost data.
    setupGhost();
    // Besides setting up our own ghost data, we also issue requests (for ghost
    // data) to those chares which we share faces with.
    for (const auto& cf : m_bndFace)
      thisProxy[ cf.first ].reqGhost( thisIndex );
  }
}

void
DG::setupGhost()
// *****************************************************************************
// Setup own ghost data on this chare
// *****************************************************************************
{
  auto d = Disc();
  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();

  // Collect tet ids, their face connectivity (given by 3 global node IDs, each
  // triplet for potentially mulitple faces on the chare boundary), and their
  // elem geometry data (see GhostData) associated to fellow chares adjacent to
  // chare boundaries. Once received by fellow chares, these tets will become
  // known as ghost elements and their data as ghost data.
  for (std::size_t e=0; e<m_esuelTet.size()/4; ++e) {  // for all our tets
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {  // for all cell faces
      if (m_esuelTet[mark+f] == -1) {  // if face has no outside-neighbor tet
        auto A = gid[ inpoel[ mark + tk::lpofa[f][0] ] ];
        auto B = gid[ inpoel[ mark + tk::lpofa[f][1] ] ];
        auto C = gid[ inpoel[ mark + tk::lpofa[f][2] ] ];
        auto c = findchare( {{A,B,C}} );
        // It is possible that we do not find the chare for this face. We are
        // looping through all of our tets and interrogating all faces that do
        // not have neighboring tets but we only care about chare-boundary faces
        // here as only those need ghost data. (m_esuelTet may also contain
        // physical boundary faces.)
        if (c > -1) {
          // Will store ghost data associated to neighbor chare
          auto& ghost = m_ghostData[ c ];
          // Store tet id adjacent to chare boundary as key for ghost data
          auto& tuple = ghost[ e ];
          // If tetid e has not yet been encountered, store geometry (only once)
          auto& nodes = std::get< 0 >( tuple );
          if (nodes.empty()) std::get< 1 >( tuple ) = m_geoElem[ e ];
          // (Always) store face node IDs on chare boundary, even if tetid e has
          // already been stored. Thus we store potentially multiple faces along
          // the same chare-boundary. This happens, e.g., when the boundary
          // between chares is zig-zaggy enough to have 2 or even 3 faces of the
          // same tet.
          nodes.push_back( A );
          nodes.push_back( B );
          nodes.push_back( C );
        }
      }
    }
  }

  // If our own ghost data is empty, we do not expect to requests, otherwise
  // tell the runtime system that we have finished setting up our ghost data.
  if (m_ghostData.empty()) adj(); else ownghost_complete();
}

void
DG::reqGhost( int fromch )
// *****************************************************************************
// Receive requests for ghost data
// *****************************************************************************
{
  // Buffer up requestor chare IDs
  m_ghostReq.push_back( fromch );

  // If every chare we communicate with has requested ghost data from us, we may
  // fulfill the requests, but only if we have already setup our ghost data.
  if (m_ghostReq.size() == m_bndFace.size()) reqghost_complete();
}

void
DG::sendGhost()
// *****************************************************************************
// Send all of our ghost data to fellow chares
// *****************************************************************************
{
  for (const auto& c : m_ghostData)
    thisProxy[ c.first ].comGhost( thisIndex, c.second );
}

int
DG::findchare( const tk::UnsMesh::Face& t )
// *****************************************************************************
// Find chare for face (given by 3 global node IDs
//! \return chare ID if found, -1 if not
// *****************************************************************************
{
  for (const auto& cf : m_bndFace)
    if (cf.second.find(t) != end(cf.second))
      return cf.first;
  return -1;
}

void
DG::comGhost( int fromch, const GhostData& ghost )
// *****************************************************************************
// Receive ghost data on chare boundaries from fellow chare
// *****************************************************************************
{
  // nodelist with fromch, currently only used for an assert
  const auto& nl = tk::cref_find( m_msumset, fromch );
  IGNORE(nl);

  // Store ghost data coming from chare
  auto ghostcnt = m_u.nunk();    // will start new local cell ids from
  for (const auto& g : ghost) {  // loop over incoming ghost data
    auto e = g.first;  // (ghost) tet id on the other side of chare boundary
    const auto& nodes = std::get< 0 >( g.second );  // node IDs of face(s)
    const auto& geo = std::get< 1 >( g.second );    // ghost elem geometry data
    Assert( nodes.size() % 3 == 0, "Face node IDs must be triplets" );
    Assert( geo.size() == m_geoElem.nprop(), "Ghost geometry size mismatch" );
    for (std::size_t n=0; n<nodes.size()/3; ++n) {  // face(s) of ghost e
      auto A = nodes[ n*3+0 ];  // global node IDs of face on chare boundary
      auto B = nodes[ n*3+1 ];
      auto C = nodes[ n*3+2 ];
      // must find face(A,B,C) in nodelist of chare-boundary adjacent to fromch
      Assert( nl.find(A)!=end(nl) && nl.find(B)!=end(nl) && nl.find(C)!=end(nl),
           "Ghost face not found in chare-node adjacency map on receiving end" );
      // must find face(A,B,C) in boundary-face adjacency map
      auto c = findchare({{A,B,C}});
      IGNORE(c);
      Assert( c > -1, "Ghost face not found in boundary-face adjacency map on "
                      "receiving end" );
      // if ghost tet id not yet encountered on boundary with fromch
      if ( m_ghost.find(e) == end(m_ghost) ) {
        m_ghost[e] = ghostcnt++;  // assign new local tet id to remote ghost id
        m_geoElem.push_back( geo );  // store ghost elem geometry
      }
    }
  }

  if (++m_nadj == m_bndFace.size()) adj();
}

void
DG::adj()
// *****************************************************************************
// Continue after face adjacency communication map is complete on this chare
// *****************************************************************************
{
//   std::cout << "\nGhosts on " << thisIndex << " (remote:local): ";
//   for (const auto& g : m_ghost) std::cout << g.first << ":" << g.second << ' ';
//   std::cout << '\n';

//std::cout << thisIndex << "b: " << m_geoElem.nunk() << ", " << m_ghost.size() << '\n';

  // Enlarge lhs, rhs, and solution to accommodate ghost cells on chare boundaries
  m_u.enlarge( m_ghost.size() );
  m_un.enlarge( m_ghost.size() );
  m_lhs.enlarge( m_ghost.size() );
  m_rhs.enlarge( m_ghost.size() );

  std::vector< int > esufch;
  esufch.resize(2*m_facecnt - m_esuf.size());

  if (!m_ghostData.empty())
  {
    for (const auto& n : m_bndFace)
    {
      auto igd = m_ghostData.find( n.first );
      if ( igd != end(m_ghostData) )
      {
        const auto& ngd = igd->second;
        for (const auto& i : ngd)
        {
          // tet-id adjacent to chare-face
          auto e = i.first;
          // node-ids of chare-face(s)
          const auto& nodes = std::get< 0 >( i.second );
          Assert( nodes.size() % 3 == 0, "Face node IDs must be triplets" );
          for (std::size_t in=0; in<nodes.size(); ++in)
          {
            auto A = nodes[ in*3+0 ];
            auto B = nodes[ in*3+1 ];
            auto C = nodes[ in*3+2 ];
            const std::array< std::size_t, 3 >& t = {{A, B, C}};

            // find if this node-triplet exists on the current m_bndFace
            auto it = n.second.find(t);
            if ( it != end(n.second) )
            {
              // a matching face in m_ghostData and m_bndFace is found
              // now esufch can be updated
              auto f = it->second - m_esuf.size()/2;

              esufch[2*f + 0] = static_cast< int >(e);
              esufch[2*f + 1] = static_cast< int >(m_ghost[e]);
            }
          }
        }
      }
    }
  }

  m_esuf.insert( std::end(m_esuf), std::begin(esufch), std::end(esufch) );

//std::cout << thisIndex << "a: " << m_geoElem.nunk() << ", " << m_lhs.nunk() << '\n';

  // Signal the runtime system that all workers have received their adjacency
  m_solver.ckLocalBranch()->created();
}

void
DG::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details Since this is a [nodeinit] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  DiagMerger = CkReduction::addReducer( tk::mergeDiag );
}

void
DG::setup( tk::real v )
// *****************************************************************************
// Setup rows, query boundary conditions, output mesh, etc.
//! \param[in] v Total mesh volume
// *****************************************************************************
{
  auto d = Disc();

  // Store total mesh volume
  m_vol = v;
  // Output chare mesh to file
  d->writeMesh();
  // Output fields metadata to output file
  d->writeElemMeta();

  // Compute left-hand side of discrete PDEs
  lhs();

  // zero initial solution vector
  // ...

  // Set initial conditions for all PDEs
  for (const auto& eq : g_dgpde) eq.initialize( m_geoElem, m_u,  d->T() );
  m_un = m_u;

  // Output initial conditions to file (regardless of whether it was requested)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() ) writeFields( d->T() );

  // Start time stepping
  dt();
}

void
DG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
  auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  // use constant dt if configured
  if (std::abs(const_dt - def_const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    // find the minimum dt across all PDEs integrated
    // ...
    mindt = 0.1;        // stub for now to overwrite numeric_limits::max

    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  }

  // Contribute to minimum dt across all chares the advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(DG,commGhostData), thisProxy) );
}

void
DG::commGhostData( tk::real newdt )
// *****************************************************************************
// Communication step to send chare-boundary ghost data to neighboring chares
// *****************************************************************************
{
  if (!m_ghostData.empty())
  {
    for(const auto& n : m_ghostData)
    {
      std::vector< std::size_t > geid;
      std::vector< std::vector< tk::real > > gd;
      for(const auto& i : n.second)
      {
        geid.push_back( i.first );
        gd.push_back( m_u[i.first] );
      }
      thisProxy[ n.first ].comrhs( geid, gd );
    }
  }

  advance( newdt );
}

void
DG::comrhs(const std::vector< std::size_t >& geid,
           const std::vector< std::vector< tk::real > >& V)
// *****************************************************************************
//  Receive chare-boundary ghost data from neighboring chares
//! \param[in] geid Global element IDs of the ghost element for which we receive
//!   data
//! \param[in] V Ghost element data
//! \details This function receives contributions to m_u, m_geoElem.
// *****************************************************************************
{
  Assert( V.size() == geid.size(), "Size mismatch in DG::comrhs()" );

  for (std::size_t i=0; i<geid.size(); ++i)
  {
    auto leid = tk::cref_find( m_ghost, geid[i] );
    Assert( leid < m_u.nunk(), "Indexing out of bounds in DG::comrhs()" );
    for (std::size_t c=0; c<m_u.nprop(); ++c)
    {
      m_u(leid, c, 0) = V[i][c];
    }
  }
}

void
DG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  d->setdt( newdt );

  // Compute rhs for next time step
  rhs();

  // Advance solution/time-stepping
  solve( newdt );

  // Prepare for next time step
  next();
}

void
DG::writeFields( tk::real time )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] time Physical time
// *****************************************************************************
{
  auto d = Disc();

  // Save time stamp at which the last field write happened
  d->LastFieldWriteTime() = time;

  // Increase field output iteration count
  ++m_itf;

  // Collect element field output
  std::vector< std::vector< tk::real > > elemfields;
  auto u = m_u;   // make a copy as eq::output() may overwrite its arg
  for (const auto& eq : g_dgpde)
  {
    auto output = eq.fieldOutput( time, m_vol, m_geoElem, u );
    for (auto& o : output) o.resize( o.size()-m_ghost.size() );
    elemfields.insert( end(elemfields), begin(output), end(output) );
  }

  // Create ExodusII writer
  tk::ExodusIIMeshWriter ew( d->OutFilename(), tk::ExoWriter::OPEN );
  // Write time stamp
  ew.writeTimeStamp( m_itf, time );
  // Write node fields to file
  d->writeElemSolution( ew, m_itf, elemfields );
}

bool
DG::diagnostics()
// *****************************************************************************
// Compute diagnostics, e.g., residuals
//! \return True if diagnostics have been computed
// *****************************************************************************
{
  auto d = Disc();

  const auto ncomp = g_inputdeck.get< tag::component >().nprop();

  std::vector< std::vector< tk::real > >
    diag( NUMDIAG, std::vector< tk::real >( ncomp, 0.0 ) );

  // Compute diagnostics
  // ...

  // Contribute to diagnostics across all PEs
  auto stream = tk::serialize( diag );
  contribute( stream.first, stream.second.get(), DiagMerger,
    CkCallback(CkIndex_Transporter::diagnostics(nullptr), d->Tr()) );

  return true;
}

void
DG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  // Optionally output field and particle data
  if ( !((d->It()+1) % g_inputdeck.get< tag::interval, tag::field >()) &&
       !g_inputdeck.get< tag::cmd, tag::benchmark >() )
  {
    writeFields( d->T() + d->Dt() );
  }

  // Output final field data to file (regardless of whether it was requested)
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  if ( (std::fabs(d->T() + d->Dt() - term) < eps || (d->It()+1) >= nstep) &&
       (!g_inputdeck.get< tag::cmd, tag::benchmark >()) )
  {
    writeFields( d->T()+d->Dt() );
  }
}

void
DG::lhs()
// *****************************************************************************
// Compute left-hand side of discrete transport equations
// *****************************************************************************
{
  // Compute left-hand side matrix for all equations
  for (const auto& eq : g_dgpde)
    eq.lhs( m_geoElem, m_lhs );
}

void
DG::rhs()
// *****************************************************************************
// Compute right-hand side of discrete transport equations
// *****************************************************************************
{
  auto d = Disc();

  for (const auto& eq : g_dgpde)
    eq.rhs( d->T(), m_geoFace, m_fd, m_u, m_rhs );
}

void
DG::solve( tk::real deltat )
// *****************************************************************************
// Explicit time-stepping using forward Euler to discretize time-derivative
// *****************************************************************************
{
  m_u = m_un + deltat * m_rhs/m_lhs;

  m_un = m_u;
}

void
DG::next()
// *****************************************************************************
// Prepare for next step
// *****************************************************************************
{
  auto d = Disc();

  // Output field data to file
  out();
  // Compute diagnostics, e.g., residuals
  //auto diag = diagnostics();
  // Increase number of iterations and physical time
  d->next();
  // Output one-liner status report
  d->status();

  // Evaluate whether to continue with next step
  /*if (!diag)*/ eval();
}

void
DG::eval()
// *****************************************************************************
// Evaluate whether to continue with next step
// *****************************************************************************
{
  auto d = Disc();

  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();

  // If neither max iterations nor max time reached, continue, otherwise finish
  if (std::fabs(d->T()-term) > eps && d->It() < nstep)
    dt();
  else
    contribute( CkCallback( CkReductionTarget(Transporter,finish), d->Tr() ) );
}

#include "NoWarning/dg.def.h"
