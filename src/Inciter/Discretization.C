// *****************************************************************************
/*!
  \file      src/Inciter/Discretization.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \details   Data and functionality common to all discretization schemes
  \see       Discretization.h and Discretization.C for more info.
*/
// *****************************************************************************

#include "Tags.h"
#include "Reorder.h"
#include "Vector.h"
#include "DerivedData.h"
#include "Discretization.h"
#include "ExodusIIMeshReader.h"
#include "ExodusIIMeshWriter.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/Options/Scheme.h"
#include "CGPDE.h"
#include "DGPDE.h"
#include "Print.h"

#ifdef HAS_ROOT
  #include "RootMeshWriter.h"
#endif

namespace inciter {

static CkReduction::reducerType PDFMerger;
extern std::vector< CGPDE > g_cgpde;
extern std::vector< DGPDE > g_dgpde;
extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::Discretization;

Discretization::Discretization(
  const CProxy_DistFCT& fctproxy,
  const CProxy_Transporter& transporter,
  const std::vector< std::size_t >& conn,
  const tk::UnsMesh::CoordMap& coordmap,
  const std::map< int, std::unordered_set< std::size_t > >& msum,
  int nchare ) :
  m_it( 0 ),
  m_t( g_inputdeck.get< tag::discr, tag::t0 >() ),
  m_dt( g_inputdeck.get< tag::discr, tag::dt >() ),
  m_lastFieldWriteTime( -1.0 ),
  m_nvol( 0 ),
  m_outFilename( g_inputdeck.get< tag::cmd, tag::io, tag::output >() + '.' +
                 std::to_string( thisIndex )
                 #ifdef HAS_ROOT
                 + (g_inputdeck.get< tag::selected, tag::filetype >() ==
                     tk::ctr::FieldFileType::ROOT ? ".root" : "")
                 #endif
                ),
  m_fct( fctproxy ),
  m_transporter( transporter ),
  m_el( tk::global2local( conn ) ),     // fills m_inpoel, m_gid, m_lid
  m_coord( setCoord( coordmap ) ),
  m_psup( tk::genPsup( m_inpoel, 4, tk::genEsup(m_inpoel,4) ) ),
  m_v( m_gid.size(), 0.0 ),
  m_vol( m_gid.size(), 0.0 ),
  m_volc(),
  m_bid(),
  m_timer()
// *****************************************************************************
//  Constructor
//! \param[in] transporter Host (Transporter) proxy
//! \param[in] fctproxy Distributed FCT proxy
//! \param[in] conn Vector of mesh element connectivity owned (global IDs)
//! \param[in] coordmap Coordinates of mesh nodes and their global IDs
//! \param[in] msum Global mesh node IDs associated to chare IDs bordering the
//!   mesh chunk we operate on
//! \param[in] nchare Total number of Discretization chares
//! \details "Contiguous-row-id" here means that the numbering of the mesh nodes
//!   (which corresponds to rows in the linear system) are (approximately)
//!   contiguous (as much as this can be done with an unstructured mesh) as the
//!   problem is distirbuted across PEs, held by Solver objects. This ordering
//!   is in start contrast with "as-in-file" ordering, which is the ordering of
//!   the mesh nodes as it is stored in the file from which the mesh is read in.
//!   The as-in-file ordering is highly non-contiguous across the distributed
//!   problem.
// *****************************************************************************
{
  Assert( m_psup.second.size()-1 == m_gid.size(),
          "Number of mesh points and number of global IDs unequal" );

  // Enable migration at AtSync()
  usesAtSync = true;

  // Convert neighbor nodes to vectors from sets
  for (const auto& n : msum) {
    auto& v = m_msum[ n.first ];
    v.insert( end(v), begin(n.second), end(n.second) );
  }

  // Count the number of mesh nodes at which we receive data from other chares
  // and compute map associating boundary-chare node ID associated to global ID
  std::vector< std::size_t > c;
  for (const auto& n : m_msum) for (auto i : n.second) c.push_back( i );
  tk::unique( c );
  m_bid = tk::assignLid( c );

  // Allocate receive buffer for nodal volumes
  m_volc.resize( m_bid.size(), 0.0 );

  // Insert DistFCT chare array element if FCT is needed. Note that even if FCT
  // is configured false in the input deck, at this point, we still need the FCT
  // object as FCT is still being performed, only its results are ignored. See
  // also, e.g., MatCG::next().
  const auto sch = g_inputdeck.get< tag::discr, tag::scheme >();
  const auto nprop = g_inputdeck.get< tag::component >().nprop();
  if ((sch == ctr::SchemeType::MatCG || sch == ctr::SchemeType::DiagCG))
    m_fct[ thisIndex ].insert( m_transporter, nchare, m_gid.size(), nprop,
                               m_msum, m_bid, m_lid, m_inpoel, CkMyPe() );

  contribute( CkCallback(CkReductionTarget(Transporter,disccreated),
              m_transporter) );
}

void
Discretization::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//!  \details Since this is a [nodeinit] routine, see cg.ci, the
//!   Charm++ runtime system executes the routine exactly once on every
//!   logical node early on in the Charm++ init sequence. Must be static as
//!   it is called without an object. See also: Section "Initializations at
//!   Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  PDFMerger = CkReduction::addReducer( tk::mergeUniPDFs );
}

tk::UnsMesh::Coords
Discretization::setCoord( const tk::UnsMesh::CoordMap& coordmap )
// *****************************************************************************
// ...
// *****************************************************************************
{
  Assert( coordmap.size() == m_gid.size(), "Size mismatch" );
  Assert( coordmap.size() == m_lid.size(), "Size mismatch" );

  tk::UnsMesh::Coords coord;
  coord[0].resize( coordmap.size() );
  coord[1].resize( coordmap.size() );
  coord[2].resize( coordmap.size() );

  for (const auto& p : coordmap) {
    auto i = tk::cref_find( m_lid, p.first );
    coord[0][i] = p.second[0];
    coord[1][i] = p.second[1];
    coord[2][i] = p.second[2];
  }

  return coord;
}

void
Discretization::vol()
// *****************************************************************************
// Sum mesh volumes to nodes, start communicating them on chare-boundaries
// *****************************************************************************
{
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  // Compute nodal volumes on our chunk of the mesh
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }};
    // compute element Jacobi determinant * 5/120 = element volume / 4
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto J = tk::triple( ba, ca, da ) * 5.0 / 120.0;
    ErrChk( J > 0, "Element Jacobian non-positive: PE:" +
                   std::to_string(CkMyPe()) + ", node IDs: " +
                   std::to_string(m_gid[N[0]]) + ',' +
                   std::to_string(m_gid[N[1]]) + ',' +
                   std::to_string(m_gid[N[2]]) + ',' +
                   std::to_string(m_gid[N[3]]) + ", coords: (" +
                   std::to_string(x[N[0]]) + ", " +
                   std::to_string(y[N[0]]) + ", " +
                   std::to_string(z[N[0]]) + "), (" +
                   std::to_string(x[N[1]]) + ", " +
                   std::to_string(y[N[1]]) + ", " +
                   std::to_string(z[N[1]]) + "), (" +
                   std::to_string(x[N[2]]) + ", " +
                   std::to_string(y[N[2]]) + ", " +
                   std::to_string(z[N[2]]) + "), (" +
                   std::to_string(x[N[3]]) + ", " +
                   std::to_string(y[N[3]]) + ", " +
                   std::to_string(z[N[3]]) + ')' );
    // scatter add V/4 to nodes
    for (std::size_t j=0; j<4; ++j) m_vol[N[j]] += J;
  }

  // Store nodal volumes without contributions from other chares on
  // chare-boundaries
  m_v = m_vol;

  // Send our nodal volume contributions to neighbor chares
  if (m_msum.empty())
    contribute( CkCallback(CkReductionTarget(Transporter,vol), m_transporter) );
  else
    for (const auto& n : m_msum) {
      std::vector< tk::real > v;
      for (auto i : n.second) v.push_back( m_vol[ tk::cref_find(m_lid,i) ] );
      thisProxy[ n.first ].comvol( n.second, v );
    }
}

void
Discretization::comvol( const std::vector< std::size_t >& gid,
                        const std::vector< tk::real >& nodevol )
// *****************************************************************************
//  Receive nodal volumes on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive volume contributions
//! \param[in] nodevol Partial sums of nodal volume contributions to
//!    chare-boundary nodes
//! \details This function receives contributions to m_vol, which stores the
//!   nodal volumes. While m_vol stores own contributions, m_volc collects the
//!   neighbor chare contributions during communication. This way work on m_vol
//!   and m_volc is overlapped. The two are combined in totalvol().
// *****************************************************************************
{
  Assert( nodevol.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( m_bid, gid[i] );
    Assert( bid < m_volc.size(), "Indexing out of bounds" );
    m_volc[ bid ] += nodevol[i];
  }

  if (++m_nvol == m_msum.size()) {
    m_nvol = 0;
    contribute( CkCallback(CkReductionTarget(Transporter,vol), m_transporter) );
  }
}

void
Discretization::totalvol()
// *****************************************************************************
// Sum mesh volumes and contribute own mesh volume to total volume
// *****************************************************************************
{
  // Combine own and communicated contributions of nodal volumes
  for (const auto& b : m_bid) {
    auto lid = tk::cref_find( m_lid, b.first );
    m_vol[ lid ] += m_volc[ b.second ];
  }

  // Sum mesh volume to host
  tk::real tvol = 0.0;
  for (auto v : m_v) tvol += v;
  contribute( sizeof(tk::real), &tvol, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,totalvol), m_transporter) );
}

void
Discretization::stat()
// *****************************************************************************
// Compute mesh cell statistics
// *****************************************************************************
{
  const auto& x = m_coord[0];
  const auto& y = m_coord[1];
  const auto& z = m_coord[2];

  auto MIN = -std::numeric_limits< tk::real >::max();
  auto MAX = std::numeric_limits< tk::real >::max();
  std::vector< tk::real > min{ MAX, MAX };
  std::vector< tk::real > max{ MIN, MIN };
  std::vector< tk::real > sum{ 0.0, 0.0, 0.0, 0.0 };
  tk::UniPDF edgePDF( 1e-4 );
  tk::UniPDF volPDF( 1e-4 );

  // Compute edge length statistics
  // Note that while the min and max edge lengths are independent of the number
  // of CPUs (by the time they are aggregated across all chares), the sum of
  // the edge lengths and the edge length PDF are not. This is because the
  // edges on the chare-boundary are counted multiple times and we
  // conscientiously do not make an effort to precisely compute this, because
  // that would require communication and more complex logic. Since these
  // statistics are intended as simple average diagnostics, we ignore these
  // small differences. For reproducible average edge lengths and edge length
  // PDFs, run the mesh in serial.
  for (std::size_t p=0; p<m_gid.size(); ++p)
    for (auto i=m_psup.second[p]+1; i<=m_psup.second[p+1]; ++i) {
       const auto dx = x[ m_psup.first[i] ] - x[ p ];
       const auto dy = y[ m_psup.first[i] ] - y[ p ];
       const auto dz = z[ m_psup.first[i] ] - z[ p ];
       const auto length = std::sqrt( dx*dx + dy*dy + dz*dz );
       if (length < min[0]) min[0] = length;
       if (length > max[0]) max[0] = length;
       sum[0] += 1.0;
       sum[1] += length;
       edgePDF.add( length );
    }

  // Compute mesh cell volume statistics
  for (std::size_t e=0; e<m_inpoel.size()/4; ++e) {
    const std::array< std::size_t, 4 > N{{ m_inpoel[e*4+0], m_inpoel[e*4+1],
                                           m_inpoel[e*4+2], m_inpoel[e*4+3] }};
    const std::array< tk::real, 3 >
      ba{{ x[N[1]]-x[N[0]], y[N[1]]-y[N[0]], z[N[1]]-z[N[0]] }},
      ca{{ x[N[2]]-x[N[0]], y[N[2]]-y[N[0]], z[N[2]]-z[N[0]] }},
      da{{ x[N[3]]-x[N[0]], y[N[3]]-y[N[0]], z[N[3]]-z[N[0]] }};
    const auto L = std::cbrt( tk::triple( ba, ca, da ) / 6.0 );
    if (L < min[1]) min[1] = L;
    if (L > max[1]) max[1] = L;
    sum[2] += 1.0;
    sum[3] += L;
    volPDF.add( L );
  }

  // Contribute to mesh statistics across all Discretization chares
  contribute( min, CkReduction::min_double,
    CkCallback(CkReductionTarget(Transporter,minstat), m_transporter) );
  contribute( max, CkReduction::max_double,
    CkCallback(CkReductionTarget(Transporter,maxstat), m_transporter) );
  contribute( sum, CkReduction::sum_double,
    CkCallback(CkReductionTarget(Transporter,sumstat), m_transporter) );

  // Serialize PDFs to raw stream
  auto stream = tk::serialize( { edgePDF, volPDF } );
  // Create Charm++ callback function for reduction of PDFs with
  // Transporter::pdfstat() as the final target where the results will appear.
  CkCallback cb( CkIndex_Transporter::pdfstat(nullptr), m_transporter );
  // Contribute serialized PDF of partial sums to host via Charm++ reduction
  contribute( stream.first, stream.second.get(), PDFMerger, cb );
}

void
Discretization::writeMesh()
// *****************************************************************************
// Output chare element blocks to file
// *****************************************************************************
{
  if (!g_inputdeck.get< tag::cmd, tag::benchmark >()) {

    #ifdef HAS_ROOT
    auto filetype = g_inputdeck.get< tag::selected, tag::filetype >();

    if (filetype == tk::ctr::FieldFileType::ROOT) {

      tk::RootMeshWriter rmw( m_outFilename, 0 );
      rmw.writeMesh( tk::UnsMesh( m_inpoel, m_coord ) );

    } else
    #endif
    {
      // Create ExodusII writer
      tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::CREATE );
      // Write chare mesh initializing element connectivity and point coords
      ew.writeMesh( tk::UnsMesh( m_inpoel, m_coord ) );
    }    
  }
}

void
Discretization::writeNodeSolution(
  const tk::ExodusIIMeshWriter& ew,
  uint64_t it,
  const std::vector< std::vector< tk::real > >& u ) const
// *****************************************************************************
// Output solution to file
//! \param[in] ew ExodusII mesh-based writer object
//! \param[in] it Iteration count
//! \param[in] u Vector of fields to write to file
// *****************************************************************************
{
  int varid = 0;
  for (const auto& f : u) ew.writeNodeScalar( it, ++varid, f );
}

#ifdef HAS_ROOT
void
Discretization::writeNodeSolution(
  const tk::RootMeshWriter& rmw,
  uint64_t it,
  const std::vector< std::vector< tk::real > >& u ) const
// *****************************************************************************
// Output solution to file
//! \param[in] rmw Root mesh-based writer object
//! \param[in] it Iteration count
//! \param[in] u Vector of fields to write to file
// *****************************************************************************
{
  int varid = 0;
  for (const auto& f : u) rmw.writeNodeScalar( it, ++varid, f );
}
#endif

void
Discretization::writeElemSolution(
  const tk::ExodusIIMeshWriter& ew,
  uint64_t it,
  const std::vector< std::vector< tk::real > >& u ) const
// *****************************************************************************
// Output solution to file
//! \param[in] ew ExodusII mesh-based writer object
//! \param[in] it Iteration count
//! \param[in] u Vector of element fields to write to file
// *****************************************************************************
{
  int varid = 0;
  for (const auto& f : u) ew.writeElemScalar( it, ++varid, f );
}

void
Discretization::writeNodeMeta() const
// *****************************************************************************
// Output mesh-based fields metadata to file
// *****************************************************************************
{
  if (!g_inputdeck.get< tag::cmd, tag::benchmark >()) {

    #ifdef HAS_ROOT
    auto filetype = g_inputdeck.get< tag::selected, tag::filetype >();

    if (filetype == tk::ctr::FieldFileType::ROOT) {
 
      tk::RootMeshWriter rmw( m_outFilename, 1 );

      // Collect nodal field output names from all PDEs
      std::vector< std::string > names;
      for (const auto& eq : g_cgpde) {
        auto n = eq.fieldNames();
        names.insert( end(names), begin(n), end(n) );
      }

      // Write node field names
      rmw.writeNodeVarNames( names );

    } else
    #endif
    {

      // Create ExodusII writer
      tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::OPEN );

      // Collect nodal field output names from all PDEs
      std::vector< std::string > names;
      for (const auto& eq : g_cgpde) {
        auto n = eq.fieldNames();
        names.insert( end(names), begin(n), end(n) );
      }

      // Write node field names
      ew.writeNodeVarNames( names );
    }

  }
}

void
Discretization::writeElemMeta() const
// *****************************************************************************
// Output mesh-based element fields metadata to file
// *****************************************************************************
{
  if (!g_inputdeck.get< tag::cmd, tag::benchmark >())
  {
    // Create ExodusII writer
    tk::ExodusIIMeshWriter ew( m_outFilename, tk::ExoWriter::OPEN );

    // Collect elemental field output names from all PDEs
    std::vector< std::string > names;
    for (const auto& eq : g_dgpde) {
      auto n = eq.fieldNames();
      names.insert( end(names), begin(n), end(n) );
    }

    // Write element field names
    ew.writeElemVarNames( names );
  }
}

void
Discretization::setdt( tk::real newdt )
// *****************************************************************************
// Set time step size
//! \param[in] newdt Size of the new time step
// *****************************************************************************
{
  m_dt = newdt;

  // Truncate the size of last time step
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  if (m_t+m_dt > term) m_dt = term - m_t;;
}

void
Discretization::next()
// *****************************************************************************
// Prepare for next step
// *****************************************************************************
{
  ++m_it;
  m_t += m_dt;
}

void
Discretization::status()
// *****************************************************************************
// Output one-liner status report
// *****************************************************************************
{
  // Query after how many time steps user wants TTY dump
  const auto tty = g_inputdeck.get< tag::interval, tag::tty >();

  if (thisIndex==0 && !(m_it%tty)) {

    const auto term = g_inputdeck.get< tag::discr, tag::term >();
    const auto t0 = g_inputdeck.get< tag::discr, tag::t0 >();
    const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
    const auto field = g_inputdeck.get< tag::interval,tag::field >();
    const auto diag = g_inputdeck.get< tag::interval, tag::diag >();
    const auto verbose = g_inputdeck.get< tag::cmd, tag::verbose >();

    // estimate time elapsed and time for accomplishment
    tk::Timer::Watch ete, eta;
    m_timer.eta( term-t0, m_t-t0, nstep, m_it, ete, eta );
 
    tk::Print print( verbose ? std::cout : std::clog );
 
    // Output one-liner
    print << std::setfill(' ') << std::setw(8) << m_it << "  "
          << std::scientific << std::setprecision(6)
          << std::setw(12) << m_t << "  "
          << m_dt << "  "
          << std::setfill('0')
          << std::setw(3) << ete.hrs.count() << ":"
          << std::setw(2) << ete.min.count() << ":"
          << std::setw(2) << ete.sec.count() << "  "
          << std::setw(3) << eta.hrs.count() << ":"
          << std::setw(2) << eta.min.count() << ":"
          << std::setw(2) << eta.sec.count() << "  ";
  
    // Augment one-liner with output indicators
    if (!(m_it % field)) print << 'F';
    if (!(m_it % diag)) print << 'D';
  
    print << std::endl;
  }

  // Migreate here if needed
  AtSync();
}

#include "NoWarning/discretization.def.h"
