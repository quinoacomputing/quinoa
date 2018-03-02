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

#include "DG.h"
#include "Discretization.h"
#include "DGPDE.h"
#include "Solver.h"
#include "DiagReducer.h"
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
  m_lhs( m_disc[thisIndex].ckLocal()->Inpoel().size()/4, 0.0 ),
  m_rhs( m_u.nunk(), m_u.nprop() )
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
  // Activate SDAG waits for face adjacency map calculation
  wait4adj();

  auto d = Disc();

  // Convert vectors to sets from inside node adjacency map, d->Msum()
  std::unordered_map< int, std::unordered_set< std::size_t > > msum_set;
  for (const auto& n : d->Msum())
    msum_set[ n.first ].insert( n.second.cbegin(), n.second.cend() );

  // Collect tet ids associated to fellow chares adjacent to chare boundaries
  std::unordered_map< int, std::vector< std::size_t > > msum_el;
  auto& belem = fd.Belem();
  for (const auto& n : msum_set) {
    for (std::size_t e=0; e<belem.size(); ++e) {
      int counter = 0;
      for (std::size_t en=0; en<4; ++en) {
        auto i = n.second.find( d->Inpoel()[ e*4+en ] );
        if (i != end(n.second)) ++counter;
      }
      // if tet has at least 3 nodes on the chare boundary, it shares a face
      if (counter == 3) msum_el[ n.first ].push_back( e );
    }
  }

  // Note that while the face adjacency map is derived from the node adjacency
  // map, the size of the face adjacency communication map (msum_el computed
  // above) does not necessarily equal to the node adjacency map (d->Msum()),
  // because while a node can be shared at a single corner or along an edge, but
  // that does not necessarily share a face as well. So the chares we
  // communicate with across faces are not necessarily the same as the chares we
  // would communicate nodes with.
  //
  // Since the sizes of the node and face adjacency maps are not the same,
  // simply sending the tet ids adjacent to chare boundaries would be okay, but
  // the receiving size would not necessarily know how many chares it must
  // receive tet ids from. To solve this problem we send to chares that which we
  // share at least a single node, which is the size of the node adjacency map,
  // d->Msum(), but we either send a tet id list which share faces on the chare
  // boundary or an empty elem list if there is not a single tet that shares a
  // face with the destination chare (only single nodes or edges). The
  // assumption here is, of course, that the size of the face adjacency map is
  // always smaller than or equal to that of the node adjacency map. Since the
  // receive side already knows how many fellow chares it must receive shared
  // node ids from, we can use that to detect completion of the number of
  // receives. This simplifies the communication pattern and code for a small
  // price of sending a few approximately empty messages (for those chare
  // boundaries that only share individual nodes but not faces).

  ownadj_complete();

  // Send tet ids adjacent to chare boundaries to fellow workers (if any)
  if (d->Msum().empty())
    comadj_complete();
  else
    for (const auto& c : d->Msum()) {
      decltype(msum_el)::mapped_type elems;
      auto e = msum_el.find( c.first );
      if (e != end(msum_el)) elems = std::move( e->second );
      thisProxy[ c.first ].comadj( thisIndex, elems );
    }

  // Compute face geometry
  m_geoFace = tk::genGeoFaceTri(fd.Ntfac(), fd.Inpofa(), d->Coord());

  // Compute element geometry
  m_geoElem = tk::genGeoElemTet(d->Inpoel(), d->Coord());
}

void
DG::comadj( int fromch, const std::vector< std::size_t >& elems )
// *****************************************************************************
// Receive tet ids on chare boundaries from fellow chare
// *****************************************************************************
{
  auto d = Disc();

  // Store tets sharing a face with our mesh chunk categorized by fellow chares
  if (!elems.empty()) {
    auto& elemlist = m_msum_el[ fromch ];
    elemlist.insert( end(elemlist), elems.cbegin(), elems.cend() );
  }

  if (++m_nadj == d->Msum().size()) comadj_complete();
}

void
DG::adj()
// *****************************************************************************
// Continue after face adjacency communication map is complete
// *****************************************************************************
{
  std::cout << "\nAdj:";
  for (const auto& c : m_msum_el) {
    std::cout << thisIndex << ": " << c.first << ": ";
    for (auto e : c.second) std::cout << e << ' ';
  }
  std::cout << '\n';

  // Signal the runtime system that all workers have received their adjacency
  //contribute( CkCallback( CkReductionTarget(Transporter,comfinal), d->Tr() ) );
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
  m_ax = 0.1;
  m_ay = 0.1;
  m_az = 0.0;

  std::size_t nelem = m_u.nunk();

  for (std::size_t e=0; e<nelem; ++e)
  {
    auto xcc = m_geoElem(e,1,0);
    auto ycc = m_geoElem(e,2,0);
    auto zcc = m_geoElem(e,3,0);
    IGNORE(zcc);

    tk::real u = 1.0 * exp( -((xcc-0.25)*(xcc-0.25) 
                            + (ycc-0.25)*(ycc-0.25))/(2.0 * 0.005) );
    m_u(e,0,0) = u; // scalar
    m_un(e,0,0) = u; // scalar
  }

  // Set initial conditions for all PDEs
  for (const auto& eq : g_dgpde) eq.initialize( d->Coord(), m_u, d->T() );

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
  elemfields.push_back( m_u.extract(0,0) );

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
  std::size_t nelem = m_u.nunk();

  for (std::size_t e=0; e<nelem; ++e)
  {
    m_lhs[e] = m_geoElem(e,0,0);
  }
}

void
DG::rhs()
// *****************************************************************************
// Compute right-hand side of discrete transport equations
// *****************************************************************************
{
  auto& esuf = m_fd.Esuf();
  auto& bface = m_fd.Bface();

  // set rhs to zero
  m_rhs.fill(0.0);

  // compute internal surface flux integrals
  for (auto f=m_fd.Nbfac(); f<m_fd.Ntfac(); ++f)
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    auto farea = m_geoFace(f,0,0);

    std::array< tk::real, 3 > fn {{ m_geoFace(f,1,0),
                                    m_geoFace(f,2,0),
                                    m_geoFace(f,3,0) }};

    // need to use tk::Fields::extract() here somehow
    std::vector< tk::real > ul { m_u(el,0,0) };
    std::vector< tk::real > ur { m_u(er,0,0) };

    //--- upwind fluxes
    auto flux = upwindFlux(ul, ur, fn);
    
    m_rhs(el,0,0) -= farea * flux[0];
    m_rhs(er,0,0) += farea * flux[0];
  }

  // compute boundary surface flux integrals

  // symmetry boundary condition
  auto bc = bface.find(1);

  if (bc != bface.end())
  {
    for (const auto& f : bc->second)
    {
      std::size_t el = static_cast< std::size_t >(esuf[2*f]);

      Assert( esuf[2*f+1] == -1,
              "outside boundary element not -1" );

      auto farea = m_geoFace(f,0,0);

      std::array< tk::real, 3 > fn {{ m_geoFace(f,1,0),
                                      m_geoFace(f,2,0),
                                      m_geoFace(f,3,0) }};

      // need to use tk::Fields::extract() here?
      std::vector< tk::real > ul { m_u(el,0,0) };
      std::vector< tk::real > ur { ul[0] };

      //--- upwind fluxes
      auto flux = upwindFlux(ul, ur, fn);
      
      m_rhs(el,0,0) -= farea * flux[0];
    }
  }

  // inlet boundary condition
  bc = bface.find(2);

  if (bc != bface.end())
  {
    for (const auto& f : bc->second)
    {
      std::size_t el = static_cast< std::size_t >(esuf[2*f]);

      Assert( esuf[2*f+1] == -1,
              "outside boundary element not -1" );

      auto farea = m_geoFace(f,0,0);

      std::array< tk::real, 3 > fn {{ m_geoFace(f,1,0),
                                      m_geoFace(f,2,0),
                                      m_geoFace(f,3,0) }};

      // need to use tk::Fields::extract() here?
      std::vector< tk::real > ul { m_u(el,0,0) };
      std::vector< tk::real > ur { 0.0 };

      //--- upwind fluxes
      auto flux = upwindFlux(ul, ur, fn);
      
      m_rhs(el,0,0) -= farea * flux[0];
    }
  }

  // outlet boundary condition
  bc = bface.find(3);

  if (bc != bface.end())
  {
    for (const auto& f : bc->second)
    {
      std::size_t el = static_cast< std::size_t >(esuf[2*f]);

      Assert( esuf[2*f+1] == -1,
              "outside boundary element not -1" );

      auto farea = m_geoFace(f,0,0);

      std::array< tk::real, 3 > fn {{ m_geoFace(f,1,0),
                                      m_geoFace(f,2,0),
                                      m_geoFace(f,3,0) }};

      // need to use tk::Fields::extract() here?
      std::vector< tk::real > ul { m_u(el,0,0) };
      std::vector< tk::real > ur { ul[0] };

      //--- upwind fluxes
      auto flux = upwindFlux(ul, ur, fn);
      
      m_rhs(el,0,0) -= farea * flux[0];
    }
  }
}

std::vector< tk::real >
DG::upwindFlux( std::vector< tk::real > ul,
                std::vector< tk::real > ur,
                std::array< tk::real, 3 > fn )
// *****************************************************************************
// Riemann solver using upwind method
//! \param[in] ul Left unknown/state vector
//! \param[in] ur Right unknown/state vector
//! \param[in] fn Face unit normal vector
//! \return Riemann solution using upwind method
// *****************************************************************************
{
    std::vector< tk::real > flux(ul.size(),0);

    // wave speed
    tk::real swave = m_ax*fn[0] + m_ay*fn[1] + m_az*fn[2];

    // upwinding
    tk::real splus  = 0.5 * (swave + fabs(swave));
    tk::real sminus = 0.5 * (swave - fabs(swave));

    flux[0] = splus * ul[0] + sminus * ur[0];

    return flux;
}

void
DG::solve( tk::real deltat )
// *****************************************************************************
// Explicit time-stepping using forward Euler to discretize time-derivative
// *****************************************************************************
{
  std::size_t nelem = m_u.nunk();

  for (std::size_t e=0; e<nelem; ++e)
  {
    m_u(e,0,0) = m_un(e,0,0) + deltat/m_lhs[e] * m_rhs(e,0,0);
  }

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
