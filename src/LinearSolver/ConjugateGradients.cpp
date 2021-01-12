// *****************************************************************************
/*!
  \file      src/LinearSolver/ConjugateGradients.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare array for distributed conjugate gradients.
  \details   Charm++ chare array for asynchronous distributed
    conjugate gradients linear solver.
*/
// *****************************************************************************

#include <numeric>
#include <iostream>     // NOT NEEDED

#include "Exception.hpp"
#include "ConjugateGradients.hpp"
#include "Vector.hpp"

using tk::ConjugateGradients;

ConjugateGradients::ConjugateGradients(
  std::size_t dof,
  const std::vector< std::size_t >& gid,
  const std::unordered_map< std::size_t, std::size_t >& lid,
  const NodeCommMap& nodecommmap,
  const CSR& A,
  const std::vector< tk::real >& x,
  const std::vector< tk::real >& b ) :
  m_gid( gid ),
  m_lid( lid ),
  m_nodeCommMap( nodecommmap ),
  m_A( A ),
  m_x( x ),
  m_b( b ),
  m_r( gid.size()*dof, 0.0 ),
  m_rc(),
  m_nr( 0 ),
  m_p( gid.size()*dof, 0.0 ),
  m_q( gid.size()*dof, 0.0 )
// *****************************************************************************
//  Constructor
//! \param[in] dof Number of scalars per unknown (degrees of freedom, DOF)
//! \param[in] gid Global node ids
//! \param[in] lid Local node ids associated to global ones
//! \param[in] nodecommmap Global mesh node IDs shared with other chares
//!   associated to their chare IDs
//! \param[in] b Right hand side of the linear system to solve in Ax=b
// *****************************************************************************
{
  Assert( m_A.rsize() == gid.size()*dof, "Size mismatch" );
  Assert( m_x.size() == gid.size()*dof, "Size mismatch" );
  Assert( m_b.size() == gid.size()*dof, "Size mismatch" );

  // compute norm of right hand side
  CkCallback normb( CkReductionTarget(ConjugateGradients,normb),thisProxy );
  dot( m_b, m_b, normb );

  // compute A * x (for the initial residual)
  thisProxy[ thisIndex ].wait4res();
  residual();

  m_A.write_as_matlab( std::cout );
}

void
ConjugateGradients::dot( const std::vector< tk::real >& a,
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

  tk::real d = 0.0;
  for (std::size_t i=0; i<a.size(); ++i)
    if (!slave(m_nodeCommMap,m_gid[i],thisIndex))
      d += a[i]*b[i];

  contribute( sizeof(tk::real), &d, CkReduction::sum_double, c );
}

void
ConjugateGradients::normb( tk::real n )
// *****************************************************************************
// Compute the norm of the right hand side
//! \param[in] n Norm of right hand side (aggregated across all chares)
// *****************************************************************************
{
  std::cout << thisIndex << " normb: " << std::sqrt(n) << '\n';
}

void
ConjugateGradients::residual()
// *****************************************************************************
//  Initiate A * x for computing the initial residual, r = b - A * x
// *****************************************************************************
{
  // Compute own contribution to r = A * x
  m_A.mult( m_x, m_r );

  // Send partial product on chare-boundary nodes to fellow chares
  if (m_nodeCommMap.empty())
    comres_complete();
  else {
    auto dof = m_A.DOF();
    for (const auto& [c,n] : m_nodeCommMap) {
      std::vector< std::vector< tk::real > > rc( n.size() );
      std::size_t j = 0;
      for (auto g : n) {
        std::vector< tk::real > nr( dof );
        auto lid = tk::cref_find( m_lid, g );
        for (std::size_t d=0; d<dof; ++d) nr[d] = m_r[ lid*dof+d ];
        rc[j++] = std::move(nr);
      }
      thisProxy[c].comres( std::vector<std::size_t>(begin(n),end(n)), rc );
    }
  }

  ownres_complete();
}

void
ConjugateGradients::comres( const std::vector< std::size_t >& gid,
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
ConjugateGradients::initres()
// *****************************************************************************
// Compute the initial residual, r = b - A * x
// *****************************************************************************
{
  // Combine own and communicated contributions to r = A * x
  auto dof = m_A.DOF();
  for (const auto& [gid,r] : m_rc) {
    auto lid = tk::cref_find( m_lid, gid );
    for (std::size_t c=0; c<dof; ++c) m_r[lid*dof+c] += r[c];
  }

  // Finish computing initial residual, r = b - A * x
  for (auto& r : m_r) r *= -1.0;
  m_r += m_b;

  std::size_t j=0;
  for (auto r : m_r) std::cout << thisIndex << ", " << m_gid[j++] << ": " << r << '\n';
}

#include "NoWarning/conjugategradients.def.h"
