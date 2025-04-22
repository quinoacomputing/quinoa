// *****************************************************************************
/*!
  \file      src/LinearSolver/BiCG.cpp
  created 2025, Christopher Long
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare array for distributed conjugate gradients.
  \details   Charm++ chare array for asynchronous distributed
    conjugate gradients linear solver.
  \see Y. Saad, Iterative Methods for Sparse Linear Systems: Second Edition,
    ISBN 9780898718003, 2003, Algorithm 6.18, conjugate gradients to solve the
    linear system A * x = b, reproduced here:

    Compute r0:=b-A*x0, p0:=r0
    Pick rhat0 = r0
    Let rho0 = rhat0^T \dot r0
    For j=0,1,..., until convergence, do
      alpha := rho_j / (rhat0, Ap_j)
      h := x_j + alpha p_j
      s := r_j - alpha A p_j
      Exit if s is sufficient
      t = As
      w = (t,s)/(t,t)   //replaces beta

      x_{j+1} := h + ws 
      r_{j+1} := s - wt 
      Exit if r_{j+1} is sufficient
      rho_{j+1} := (rhat0, r_{j+1})  
      beta := (rho_{j+1}/rho_{j})*(alpha/w)  
      p_{j+1} := r_{j+1} + beta(p_j - wAp_j)
    end
*/
// *****************************************************************************

#include <numeric>
#include <iostream>

#include "Exception.hpp"
#include "BiCG.hpp"
#include "Vector.hpp"

using tk::BiCG;

BiCG::BiCG(
  const CSR& A,
  const std::vector< tk::real >& x,
  const std::vector< tk::real >& b,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const NodeCommMap& nodecommmap ) :
  m_A( A ),
  m_x( x ),
  m_b( b ),
  m_gid( gid ),
  m_lid( lid ),
  m_nodeCommMap( nodecommmap ),
  m_r( m_A.rsize(), 0.0 ),
  m_r0( m_A.rsize(), 0.0 ),
  m_rc(),
  m_nr( 0 ),
  m_bc(),
  m_bcc(),
  m_bcmask( m_A.rsize(), 1.0 ),
  m_nb( 0 ),
  m_p( m_A.rsize(), 0.0 ),
  m_q( m_A.rsize(), 0.0 ),
  m_t( m_A.rsize(), 0.0 ),
  m_qc(),
  m_tc(), //can I reuse this?
  m_nt( 0 ),
  m_nq( 0 ),
  m_initres(),
  m_solved(),
  m_normb( 0.0 ),
  m_it( 0 ),
  m_maxit( 0 ),
  m_rho( 0.0 ),
  m_rho0( 0.0 ),
  m_alpha( 0.0 ),
  m_omega( 0.0 ),
  m_converged( false ),
  m_xc(),
  m_x2c(),
  m_nx( 0 ),
  m_nx2( 0 )
// *****************************************************************************
//  Constructor
//! \param[in] A Left hand side matrix of the linear system to solve in Ax=b
//! \param[in] x Solution (initial guess) of the linear system to solve in Ax=b
//! \param[in] b Right hand side of the linear system to solve in Ax=b
//! \param[in] gid Global node ids
//! \param[in] lid Local node ids associated to global ones
//! \param[in] nodecommmap Global mesh node IDs shared with other chares
//!   associated to their chare IDs
// *****************************************************************************
{
  // Fill in gid and lid for serial solve
  if (gid.empty() || lid.empty() || nodecommmap.empty()) {
    m_gid.resize( m_A.rsize()/m_A.Ncomp() );
    std::iota( begin(m_gid), end(m_gid), 0 );
    for (auto g : m_gid) m_lid[g] = g;
  }

  Assert( m_A.rsize() == m_gid.size()*A.Ncomp(), "Size mismatch" );
  Assert( m_x.size() == m_gid.size()*A.Ncomp(), "Size mismatch" );
  Assert( m_b.size() == m_gid.size()*A.Ncomp(), "Size mismatch" );
}

void
BiCG::setup( CkCallback c )
// *****************************************************************************
//  Setup solver
//! \param[in] c Call to continue with after initialization is complete
//! \details This function initiates computing the residual (r=b-A*x), its dot
//!   product, and the rhs norm.
// *****************************************************************************
{
  m_initres = c;

  // initiate computing A * x (for the initial residual)
  thisProxy[ thisIndex ].wait4res();
  residual();

  // initiate computing norm of right hand side
  dot( m_b, m_b,
       CkCallback( CkReductionTarget(BiCG,normb), thisProxy ) );
}

void
BiCG::dot( const std::vector< tk::real >& a,
                         const std::vector< tk::real >& b,
                         CkCallback c )
// *****************************************************************************
//  Initiate computation of dot product of two vectors
//! \param[in] a 1st vector of dot product
//! \param[in] b 2nd vector of dot product
//! \param[in] c Callback to target with the final result
// *****************************************************************************
{
  Assert( a.size() == b.size(), "Size mismatch" );

  tk::real D = 0.0;
  auto ncomp = m_A.Ncomp();
  for (std::size_t i=0; i<a.size()/ncomp; ++i) {
    auto incomp = i*ncomp;
    if (not slave(m_nodeCommMap,m_gid[i],thisIndex))
      for (std::size_t d=0; d<ncomp; ++d)
        D += a[incomp+d] * b[incomp+d];
  }

  contribute( sizeof(tk::real), &D, CkReduction::sum_double, c );
}

void
BiCG::normb( tk::real n )
// *****************************************************************************
// Compute the norm of the right hand side
//! \param[in] n Norm of right hand side (aggregated across all chares)
// *****************************************************************************
{
  m_normb = std::sqrt(n);
  normb_complete();
}

void
BiCG::residual()
// *****************************************************************************
//  Initiate A * x for computing the initial residual, r = b - A * x
// *****************************************************************************
{
  // Compute own contribution to r = A * x
  m_A.mult( m_x, m_r, m_bcmask );

  // Send partial product on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comres_complete();
  } else {
    auto ncomp = m_A.Ncomp();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > rc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nr( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nr[d] = m_r[ i*ncomp+d ];
        rc[j++] = std::move(nr);
      }
      thisProxy[c].comres( std::vector<std::size_t>(begin(n),end(n)), rc );
    }
  }

  ownres_complete();
}

void
BiCG::comres( const std::vector< std::size_t >& gid,
                            const std::vector< std::vector< tk::real > >& rc )
// *****************************************************************************
//  Receive contributions to A * x on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] rc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( rc.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i)
    m_rc[ gid[i] ] += rc[i];

  if (++m_nr == m_nodeCommMap.size()) {
    m_nr = 0;
    comres_complete();
  }
}

void
BiCG::initres()
// *****************************************************************************
// Finish computing the initial residual, r = b - A * x
// *****************************************************************************
{
  // Combine own and communicated contributions to r = A * x
  auto ncomp = m_A.Ncomp();
  for (const auto& [gid,r] : m_rc) {
    auto i = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<ncomp; ++c) m_r[i*ncomp+c] += r[c];
  }
  tk::destroy( m_rc );

  // Finish computing the initial residual, r = b - A * x
  for (auto& r : m_r) r *= -1.0;
  m_r += m_b;

  // Initialize p
  m_p = m_r;
  //Initialize r0
  m_r0 = m_r;
  // initiate computing the dot product of the initial residual, rho = (r,r)
  dot( m_r0, m_r,
       CkCallback( CkReductionTarget(BiCG,rho), thisProxy ) );
}

void
BiCG::rho( tk::real r )
// *****************************************************************************
// Compute rho = (r,r)
//! \param[in] r Dot product, rho = (r,r) (aggregated across all chares)
// *****************************************************************************
{
  // store dot product of residual
  m_rho = r;

  // send back rhs norm to caller
  m_initres.send( CkDataMsg::buildNew( sizeof(tk::real), &m_normb ) );
}

void
BiCG::init(
  const std::vector< tk::real >& x,
  const std::vector< tk::real >& b,
  const std::unordered_map< std::size_t,
          std::vector< std::pair< bool, tk::real > > >& bc,
  std::size_t ignorebc,
  CkCallback cb )
// *****************************************************************************
//  Initialize linear solve: set initial guess and boundary conditions
//! \param[in] x Initial guess
//! \param[in] b Right hand side vector
//! \param[in] bc Local node ids and associated Dirichlet BCs
//! \param[in] ignorebc True if applyin BCs should be skipped
//! \param[in] cb Call to continue with when initialized and ready for a solve
//! \details This function allows setting the initial guess and boundary
//!   conditions, followed by computing the initial residual and the rhs norm.
// *****************************************************************************
{
  // Optionally set initial guess
  if (not x.empty()) m_x = x;

  // Optionally update rhs
  if (not b.empty()) m_b = b;

  if (ignorebc) {

    setup( cb );

  } else {

    // Store incoming BCs
    m_bc = bc;

    // Get ready to communicate boundary conditions. This is necessary because
    // there can be nodes a chare contributes to but does not apply BCs on. This
    // happens if a node is in the node communication map but not on the list of
    // incoming BCs on this chare. To have all chares share the same view on all
    // BC nodes, we send the global node ids together with the Dirichlet BCs at
    // which BCs are set to those fellow chares that also contribute to those BC
    // nodes. Only after this communication step we apply the BCs on the matrix,
    // which then will correctly setup the BC rows that exist on multiple chares
    // (which now will be the same as the results of making the BCs consistent
    // across all chares that contribute.
    thisProxy[ thisIndex ].wait4bc();

    // Send boundary conditions to those who contribute to those rows
    if (m_nodeCommMap.empty()) {
      combc_complete();
    } else {
      for (const auto& [c,n] : m_nodeCommMap) {
        std::unordered_map< std::size_t,
          std::vector< std::pair< bool, tk::real > > > expbc;
        for (auto g : n) {
          auto i = tk::cref_find( m_lid, g );
          auto j = bc.find(i);
          if (j != end(bc)) expbc[g] = j->second;
        }
        thisProxy[c].combc( expbc );
      }
    }

    ownbc_complete( cb );

  }
}

void
BiCG::combc(
  const std::unordered_map< std::size_t,
     std::vector< std::pair< bool, tk::real > > >& bc )
// *****************************************************************************
//  Receive contributions to boundary conditions on chare-boundaries
//! \param[in] bc Contributions to boundary conditions
// *****************************************************************************
{
  for (const auto& [g,dirbc] : bc) m_bcc[ tk::cref_find(m_lid,g) ] = dirbc;

  if (++m_nb == m_nodeCommMap.size()) {
    m_nb = 0;
    combc_complete();
  }
}

void
BiCG::apply( CkCallback cb )
// *****************************************************************************
//  Apply boundary conditions
//! \param[in] cb Call to continue with after applying the BCs is complete
// *****************************************************************************
{
  // Merge own and received contributions to boundary conditions
  for (const auto& [i,dirbc] : m_bcc) m_bc[i] = dirbc;
  tk::destroy( m_bcc );

  auto ncomp = m_A.Ncomp();

  // Setup Dirichlet BC map as contiguous mask
  for (const auto& [i,bc] : m_bc)
    for (std::size_t j=0; j<ncomp; ++j)
      m_bcmask[i*ncomp+j] = 0.0;

  // Apply Dirichlet BCs on matrix and rhs
  for (const auto& [i,dirbc] : m_bc) {
    for (std::size_t j=0; j<ncomp; ++j) {
      if (dirbc[j].first) {
        m_A.dirichlet( i, m_gid, m_nodeCommMap, j );
        m_b[i*ncomp+j] = dirbc[j].second;
      }
    }
  }

  // Recompute initial residual (r=b-A*x), its dot product, and the rhs norm
  setup( cb );
}

void
BiCG::solve( std::size_t maxit, tk::real tol, CkCallback c )
// *****************************************************************************
//  Solve linear system
//! \param[in] maxit Max iteration count
//! \param[in] tol Stop tolerance
//! \param[in] c Call to continue with after solve is complete
// *****************************************************************************
{
  m_maxit = maxit;
  m_tol = tol;
  m_solved = c;
  m_it = 0;

  next();
}

void
BiCG::next()
// *****************************************************************************
//  Start next linear solver iteration
// *****************************************************************************
{
  if (m_it == 0) {
     m_alpha = 0.0;
     m_omega=0.0;
     
  }else{
     m_alpha =( m_rho/m_rho0 ) * ( m_alpha/m_omega ) ; //alpha functions as Beta??
  }
  m_rho0 = m_rho;

  // compute p = r + alpha * p
  for (std::size_t i=0; i<m_p.size(); ++i) m_p[i] = m_r[i] + m_alpha * ( m_p[i] - m_omega * m_q[i] );

  // initiate computing q = A * p
  thisProxy[ thisIndex ].wait4q();
  qAp();
}


void
BiCG::qAp()
// *****************************************************************************
//  Initiate computing q = A * p
// *****************************************************************************
{
  // Compute own contribution to q = A * p
  m_A.mult( m_p, m_q, m_bcmask );

  // Send partial product on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comq_complete();
  } else {
    auto ncomp = m_A.Ncomp();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > qc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nq( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nq[d] = m_q[ i*ncomp+d ];
        qc[j++] = std::move(nq);
      }
      thisProxy[c].comq( std::vector<std::size_t>(begin(n),end(n)), qc );
    }
  }

  ownq_complete();
}

void
BiCG::comq( const std::vector< std::size_t >& gid,
                          const std::vector< std::vector< tk::real > >& qc )
// *****************************************************************************
//  Receive contributions to q = A * p on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] qc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( qc.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i)
    m_qc[ gid[i] ] += qc[i];

  if (++m_nq == m_nodeCommMap.size()) {
    m_nq = 0;
    comq_complete();
  }
}

void
BiCG::q()
// *****************************************************************************
// Finish computing q = A * p
// *****************************************************************************
{
  // Combine own and communicated contributions to q = A * p
  auto ncomp = m_A.Ncomp();
  for (const auto& [gid,q] : m_qc) {
    auto i = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<ncomp; ++c)
      m_q[i*ncomp+c] += q[c];
  }
  tk::destroy( m_qc );

  //BiCGStab uses (rhat_0,q), initiate here
  dot( m_r0, m_q,
       CkCallback( CkReductionTarget(BiCG,pq), thisProxy ) );
}

void
BiCG::pq( tk::real d )
// *****************************************************************************
// Compute the dot product (p,q)
//! \param[in] d Dot product of (p,q) (aggregated across all chares)
// *****************************************************************************
{
  // If (p,q)=0, then p and q are orthogonal and the system either has a trivial
  // solution, x=x0, or the BCs are incomplete or wrong, in either case the
  // solve cannot continue.
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  if (std::abs(d) < eps) {
    m_it = m_maxit;
    m_alpha = 0.0;
  } else {
    m_alpha = m_rho / d;
  }

  // compute s = r - alpha * q
  for (std::size_t i=0; i<m_r.size(); ++i) m_r[i] -= m_alpha * m_q[i];

  // initiate computing norm of residual: (r,r)
  dot( m_r, m_r,
       CkCallback( CkReductionTarget(BiCG,normres), thisProxy ) );
}

void
BiCG::normres( tk::real r )
// *****************************************************************************
// Compute norm of residual: (r,r)
//! \param[in] r Dot product, (r,r) (aggregated across all chares)
// *****************************************************************************
{
  m_rho = r;

  // Advance solution: h = x + alpha * p
  for (std::size_t i=0; i<m_x.size(); ++i) m_x[i] += m_alpha * m_p[i];

  // Communicate solution
  thisProxy[ thisIndex ].wait4x();

  // Send solution on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comx_complete();
  } else {
    auto ncomp = m_A.Ncomp();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > xc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nx( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nx[d] = m_x[ i*ncomp+d ];
        xc[j++] = std::move(nx);
      }
      thisProxy[c].comx( std::vector<std::size_t>(begin(n),end(n)), xc );
    }
  }

  ownx_complete();
}

void
BiCG::comx( const std::vector< std::size_t >& gid,
                          const std::vector< std::vector< tk::real > >& xc )
// *****************************************************************************
//  Receive contributions to final solution on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] xc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( xc.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) m_xc[ gid[i] ] += xc[i];

  if (++m_nx == m_nodeCommMap.size()) {
    m_nx = 0;
    comx_complete();
  }
}

void
BiCG::x()
// *****************************************************************************
// Assemble solution on chare boundaries
// *****************************************************************************
{
  // Assemble solution on chare boundaries by averaging
  auto ncomp = m_A.Ncomp();
  for (const auto& [g,x] : m_xc) {
    auto i = tk::cref_find(m_lid,g);
    for (std::size_t d=0; d<ncomp; ++d) m_x[i*ncomp+d] += x[d];
    auto c = tk::count(m_nodeCommMap,g);
    for (std::size_t d=0; d<ncomp; ++d) m_x[i*ncomp+d] /= c;
  }
  tk::destroy( m_xc );

  //++m_it;//Don't iterate counter yet!
  auto normb = m_normb > 1.0e-14 ? m_normb : 1.0;
  auto normr = std::sqrt( m_rho );

  if ( m_it < m_maxit && normr > m_tol*normb ) { //If we are not solved, continue, else exit
                                                 //Here is where we significantly diverge from regular CG as we dont' go to "next()" yet
    //next();
  // initiate computing t = A * s  ;

   // std::cout << "Tol not met:  "  << normr << " in first pass, going to second stage \n";
    thisProxy[ thisIndex ].wait4t();
    tAs();

  } else {
   // std::cout << "Tol Met:  "  << normr << " in first pass \n";
    m_converged = m_it == m_maxit && normr > m_tol*normb ? false : true;
    m_solved.send( CkDataMsg::buildNew( sizeof(tk::real), &normr ) );

  }
}

void
BiCG::tAs()
// *****************************************************************************
//  Initiate computing t = A * s
// *****************************************************************************
{
  // Compute own contribution to t = A * s
  m_A.mult( m_r, m_t, m_bcmask );

  // Send partial product on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comt_complete();
  } else {
    auto ncomp = m_A.Ncomp();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > tc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nt( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nt[d] = m_t[ i*ncomp+d ];
        tc[j++] = std::move(nt);
      }
      thisProxy[c].comt( std::vector<std::size_t>(begin(n),end(n)), tc );
    }
  }

  ownt_complete();
}

void
BiCG::comt( const std::vector< std::size_t >& gid,
                          const std::vector< std::vector< tk::real > >& tc )
// *****************************************************************************
//  Receive contributions to t = A * s on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] tc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( tc.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  for (std::size_t i=0; i<gid.size(); ++i)
    m_tc[ gid[i] ] += tc[i];

  if (++m_nt == m_nodeCommMap.size()) {
    m_nt = 0;
    comt_complete();
  }
}


void
BiCG::t()
// *****************************************************************************
// Finish computing t = A * s
// *****************************************************************************
{
  // Combine own and communicated contributions to t = A * s
  auto ncomp = m_A.Ncomp();
  for (const auto& [gid,t] : m_tc) {
    auto i = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<ncomp; ++c)
      m_t[i*ncomp+c] += t[c];
  }
  tk::destroy( m_tc );

  //omega = (t,s)/(t,t).  Compute numerator
  dot( m_t, m_r,
       CkCallback( CkReductionTarget(BiCG,ts), thisProxy ) );
}
void
BiCG::ts( tk::real d )
// *****************************************************************************
// Compute the dot product (t,s)
//! \param[in] d Dot product of (t,s) (aggregated across all chares)
// *****************************************************************************
{
  // If (t,s)=0, then p and q are orthogonal and the system either has a trivial
  // solution, x=x0, or the BCs are incomplete or wrong, in either case the
  // solve cannot continue.
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  if (std::abs(d) < eps) {
    m_it = m_maxit;
    m_alpha = 0.0;
  }
  else {
    m_omega = d; //numerator set
  }
  dot( m_t, m_t, //compute denominator
       CkCallback( CkReductionTarget(BiCG,tt), thisProxy ) );

//  // compute s = r - alpha * q
//  for (std::size_t i=0; i<m_r.size(); ++i) m_r[i] -= m_alpha * m_q[i];

  // initiate computing norm of residual: (r,r)
  //dot( m_r, m_r,
  //     CkCallback( CkReductionTarget(BiCG,normres), thisProxy ) );
}
void
BiCG::tt( tk::real d )
// *****************************************************************************
// Compute the dot product (t,t)
//! \param[in] d Dot product of (t,t) (aggregated across all chares)
// *****************************************************************************
{
  // If (t,t)=0, then p and q are orthogonal and the system either has a trivial
  // solution, x=x0, or the BCs are incomplete or wrong, in either case the
  // solve cannot continue.
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  if (std::abs(d) < eps) {
    m_it = m_maxit;
    m_alpha = 0.0;
   // std::cout << "Uhoh..  bad alpha/omega? \n";
  }
  else {
    m_omega = m_omega/d; //d is denominator
  }
  //std::cout << "Omega is: " << m_omega << "\n";
  //compute x = h + omega * s
  for (std::size_t i=0; i<m_x.size(); ++i) m_x[i] += m_omega * m_r[i];
  // Now update r:  compute r = s - omega * t
  for (std::size_t i=0; i<m_r.size(); ++i) m_r[i] -= m_omega * m_t[i];
/*  for (std::size_t i = 0; i < m_r.size(); ++i) {
	  std::cout << i << "  " << m_r[i] << " " << m_r0[i] << "\n" ;
  }
*/	  // initiate computing norm of residual: (r0,r) for assigning to new rho
  dot( m_r0, m_r,
       CkCallback( CkReductionTarget(BiCG,normresomega), thisProxy ) );
}
void
BiCG::normresomega( tk::real r )
// *****************************************************************************
// Compute norm of residual: (r,r) for second half of bicgstab
//! \param[in] r Dot product, (r,r) (aggregated across all chares)
// *****************************************************************************
{
  m_rho = r;
  
  // Communicate solution
  thisProxy[ thisIndex ].wait4x2();

  // Send solution on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty()) {
    comx2_complete();  //Does this need to be unique?
  } else {
    auto ncomp = m_A.Ncomp();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > x2c( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nx2( ncomp );
        auto i = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<ncomp; ++d) nx2[d] = m_x[ i*ncomp+d ];
        x2c[j++] = std::move(nx2);
      }
      thisProxy[c].comx2( std::vector<std::size_t>(begin(n),end(n)), x2c );
    }
  }

  ownx2_complete();
}
void
BiCG::x2()
// *****************************************************************************
// Assemble solution on chare boundaries
// *****************************************************************************
{
  // Assemble solution on chare boundaries by averaging
  auto ncomp = m_A.Ncomp();
  for (const auto& [g,x] : m_x2c) {
    auto i = tk::cref_find(m_lid,g);
    for (std::size_t d=0; d<ncomp; ++d) m_x[i*ncomp+d] += x[d];
    auto c = tk::count(m_nodeCommMap,g);
    for (std::size_t d=0; d<ncomp; ++d) m_x[i*ncomp+d] /= c;
  }
  tk::destroy( m_x2c );

  ++m_it;
  auto normb = m_normb > 1.0e-14 ? m_normb : 1.0;
  auto normr = std::sqrt( m_rho );

  if ( m_it < m_maxit && normr > m_tol*normb ) { //If we are not solved, continue, else exit

    //std::cout << "Tol not met in 2nd pass, iterating:  " << m_it << "  " << normr << " \n";
    next();

  } else {

    //std::cout << "Tol met in 2nd pass! Iteration count:  " << m_it << "  " << normr << " \n";
    m_converged = m_it == m_maxit && normr > m_tol*normb ? false : true;
    m_solved.send( CkDataMsg::buildNew( sizeof(tk::real), &normr ) );

  }
}

void
BiCG::comx2( const std::vector< std::size_t >& gid,
                          const std::vector< std::vector< tk::real > >& x2c )
// *****************************************************************************
//  Receive contributions to final solution on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] xc Partial contributions at chare-boundary nodes
// *****************************************************************************
{
  Assert( x2c.size() == gid.size(), "Size mismatch" );

  for (std::size_t i=0; i<gid.size(); ++i) m_x2c[ gid[i] ] += x2c[i];

  if (++m_nx2 == m_nodeCommMap.size()) {
    m_nx2 = 0;
    comx2_complete();
  }
}
#include "NoWarning/bicg.def.h"
