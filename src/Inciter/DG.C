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
  m_msumset( msumset() ),
  m_ghost(),
  m_chBndFace()
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
  // Activate SDAG waits for face adjacency map (ghost data) calculation
  wait4adj();

  auto d = Disc();

  const auto& gid = d->Gid();
  const auto& inpoel = d->Inpoel();
  const auto& inpofa = fd.Inpofa();
  auto esup = tk::genEsup( inpoel, 4 );
  auto esuel = tk::genEsuelTet( inpoel, esup );

  // Invert inpofa to enable searching for faces based on (global) node triplets
  Assert( inpofa.size() % 3 == 0, "Inpofa must contain triplets" );
  Faces faces;
  for (std::size_t f=0; f<inpofa.size()/3; ++f) {
    auto A = gid[ inpofa[f*3+0] ];
    auto B = gid[ inpofa[f*3+1] ];
    auto C = gid[ inpofa[f*3+2] ];
    faces.insert( {{{A,B,C}}} );
  }

  // At this point faces should have a set of node-id-triplets (faces) of
  // internal faces and physical boundary faces but not chare boundary faces.

  // Build map associating face id to (global) node ID triplets on chare boundary
  auto facecnt = faces.size();  // will start new face IDs from
  for (std::size_t e=0; e<esuel.size()/4; ++e) {
    auto mark = e*4;
    for (std::size_t f=0; f<4; ++f) {
      if (esuel[mark+f] == -1) {
        auto A = gid[ inpoel[ mark + tk::lpofa[f][0] ] ];
        auto B = gid[ inpoel[ mark + tk::lpofa[f][1] ] ];
        auto C = gid[ inpoel[ mark + tk::lpofa[f][2] ] ];

        // if does not exist in inpofa, assign new face ID on chare boundary
        NodeTriplet t{{ A, B, C }};
        if (faces.find( t ) == end(faces)) {
          m_chBndFace[ t ] = facecnt++;

//  std::cout << thisIndex << ": {" << A << ',' << B << ',' << C << '}'  << '\n';

//           // Attempt to find t on one of our chare boundaries based on msum_set
//           bool found = false;
//           for (const auto& n : m_msumset) {  // for all neighbor chares
//             auto i = n.second.find( A );
//             auto j = n.second.find( B );
//             auto k = n.second.find( C );
//             // if all face nodes are on chare boundary
//             if ( i != end(n.second) && j != end(n.second) && k != end(n.second) )
//              found = true;
//           }
//           Assert( found, "Not in msum_set: " + std::to_string(t[0]) + ',' +
//                     std::to_string(t[1]) + ',' + std::to_string(t[2]) + '\n' );
        }

      }
    }
  }

//  std::cout << thisIndex << ": " << facecnt-faces.size() << '\n';

  // At this point m_chBndFace should have new (local) face IDs assigned to
  // node-triplets (faces) only along our chare-boundary.

//   std::cout << thisIndex << "bndface: ";
//   for (const auto& f : m_chBndFace)
//     std::cout << f.first[0] << ',' << f.first[1] << ',' << f.first[2] << ' ';
//   std::cout << '\n';


//std::cout << thisIndex << "c: " << m_geoElem.nunk() << '\n';

  // Collect tet ids, their face connectivity (3 global node IDs for potentially
  // mulitple faces on the chare boundary), and their elem geometry data (see
  // GhostData) associated to fellow chares adjacent to chare boundaries. Once
  // received by fellow chares, these tets will become known as ghost elements.
  std::unordered_map< int, GhostData > msum_el;
  for (const auto& n : m_msumset) {  // for all neighbor chares
    for (std::size_t e=0; e<esuel.size()/4; ++e) {  // for all cells in our chunk
      auto mark = e*4;
      for (std::size_t f=0; f<4; ++f) {  // for all cell faces
        if (esuel[mark+f] == -1) {  // if face has no tet on the other side
          // get global node IDs of face
          auto A = gid[ inpoel[ mark + tk::lpofa[f][0] ] ];
          auto B = gid[ inpoel[ mark + tk::lpofa[f][1] ] ];
          auto C = gid[ inpoel[ mark + tk::lpofa[f][2] ] ];
          auto i = n.second.find( A );
          auto j = n.second.find( B );
          auto k = n.second.find( C );
          // if all face nodes are on chare boundary
          if ( i != end(n.second) && j != end(n.second) && k != end(n.second) ) {
            // will store ghost associated to neighbor chare
            auto& ghost = msum_el[ n.first ];
            // store tet id adjacent to chare boundary as key for ghost data
            auto& tuple = ghost[ e ];
            // if e has not yet been encountered, store geometry (only once)
            auto& nodes = std::get< 0 >( tuple );
            if (nodes.empty()) std::get< 1 >( tuple ) = m_geoElem[ e ];
            // (always) store face node IDs on chare boundary, even if we have e
            nodes.push_back( A );
            nodes.push_back( B );
            nodes.push_back( C );
            // The search below should succeed, because we are looking for node
            // triplets (faces) along our chare boundary and m_chBndFace should
            // contain all faces along our chare boundary (adjacent to all
            // chares we communicate with). However, this is not the case.
            auto bf = m_chBndFace.find( {{A,B,C}} );
            if (bf == end(m_chBndFace)) std::cout << A << ',' << B << ',' << C << " not found on sending chare " << thisIndex << '\n';
          }
        }
      }
    }
  }

  ownadj_complete();

  // Note that while the face adjacency map is derived from the node adjacency
  // map, the size of the face adjacency communication map (msum_el computed
  // above) does not necessarily equal to the node adjacency map (d->Msum()),
  // because while a node can be shared at a single corner or along an edge,
  // that does not necessarily share a face as well (in other words, shared
  // nodes or edges can exist that are not part of a shared face). So the chares
  // we communicate with across faces are not necessarily the same as the chares
  // we would communicate nodes with.
  //
  // Since the sizes of the node and face adjacency maps are not the same,
  // simply sending the tet ids adjacent to chare boundaries would be okay, but
  // the receiving size would not necessarily know how many chares it must
  // receive tet ids from. To solve this problem we send to chares that which we
  // share at least a single node, which is the size of the node adjacency map,
  // d->Msum(), but we either send a list of ghosts which share faces on the
  // chare boundary or an empty ghost list if there is not a single tet that
  // shares a face with the destination chare (only single nodes or edges). The
  // assumption here is, of course, that the size of the face adjacency map is
  // always smaller than or equal to that of the node adjacency map. Since the
  // receive side already knows how many fellow chares it must receive shared
  // node ids from, we can use that to detect completion of the number of
  // receives. This simplifies the communication pattern and code for a small
  // price of sending a few approximately empty messages (for those chare
  // boundaries that only share individual nodes but not faces).

  // Send ghost data adjacent to chare boundaries to fellow workers (if any)
  if (d->Msum().empty())
    comadj_complete();
  else
    for (const auto& c : d->Msum()) {
      decltype(msum_el)::mapped_type ghost;
      auto e = msum_el.find( c.first );
      if (e != end(msum_el)) ghost = std::move( e->second );
      thisProxy[ c.first ].comadj( thisIndex, ghost );
    }
}

std::unordered_map< int, std::unordered_set< std::size_t > >
DG::msumset() const
// *****************************************************************************
// Convert vectors to sets inside node adjacency map, Discretization::m_msum
// *****************************************************************************
{
  auto d = Disc();

  std::unordered_map< int, std::unordered_set< std::size_t > > m;
  for (const auto& n : d->Msum())
    m[ n.first ].insert( n.second.cbegin(), n.second.cend() );

  std::cout << thisIndex << "msum_set: ";
  for (const auto& n : m_msumset)
    for (auto i : n.second)
      std::cout << i << ' ';
  std::cout << '\n';

  return m;
}

void
DG::comadj( int fromch, const GhostData& ghost )
// *****************************************************************************
// Receive ghost data on chare boundaries from fellow chare
// *****************************************************************************
{
  auto d = Disc();

  // Store ghosts sharing a face with our mesh chunk categorized by fellow chares
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
      const auto& nl = tk::cref_find( m_msumset, fromch );// nodelist with fromch
      Assert( nl.find(A)!=end(nl) && nl.find(B)!=end(nl) && nl.find(C)!=end(nl),
              "Ghost face not found on receiving end" );
      // add new tet as element surrounding face(A,B,C)
      // The search below should succeed, because we are looking for node
      // triplets (faces) along our chare boundary and m_chBndFace should
      // contain all faces along our chare boundary (adjacent to all chares we
      // communicate with). However, this is not the case.
      auto bf = m_chBndFace.find( {{A,B,C}} );
      if (bf == end(m_chBndFace)) std::cout << A << ',' << B << ',' << C << " not found on receiving chare " << thisIndex << '\n';
      //auto f = tk::cref_find( m_chBndFace, {{A,B,C}} );
      //m_fd.Esuf()
      // if ghost tet id not yet encountered on boundary with fromch
      if ( m_ghost.find(e) == end(m_ghost) ) {
        m_ghost[e] = ghostcnt++;  // assign new local tet id to remote ghost id
        m_geoElem.push_back( geo );  // store ghost elem geometry
      }
    }
  }

  if (++m_nadj == d->Msum().size()) comadj_complete();
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
              CkCallback(CkReductionTarget(DG,advance), thisProxy) );
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
