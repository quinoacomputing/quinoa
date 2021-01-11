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

using tk::ConjugateGradients;

ConjugateGradients::ConjugateGradients(
  std::size_t size,
  std::size_t dof,
  const std::pair< std::vector< std::size_t >,
                   std::vector< std::size_t > >& psup,
  const std::vector< tk::real >& b ) :
  m_A( dof, psup ),
  m_x( size*dof, 0.0 ),
  m_r( size*dof, 0.0 ),
  m_p( size*dof, 0.0 ),
  m_q( size*dof, 0.0 )
// *****************************************************************************
//  Constructor
//! \param[in] size Number of unknowns on this chare
//! \param[in] dof Number of scalars per unknown (degrees of freedom, DOF)
//! \param[in] psup Points surrounding points
//! \param[in] b Right hand side in Ax=b
// *****************************************************************************
{
  Assert( b.size() == size*dof, "Size mismatch" );
  Assert( m_A.rsize() == size*dof, "Size mismatch" );

  // compute norm of right hand side
  dot(b, b, CkCallback(CkReductionTarget(ConjugateGradients,normb),thisProxy));

  // // compute initial residual and norm of right hand side
  // tk::real normb = 0.0;
  // for (std::size_t i=0; i<m_A.rsize(); ++i) {
  //   normb += b[i] * b[i];
  //   // r = b - A * x
  //   A->r[i] = b[i];
  //   for (std::size_t j=A->ia[i]-1; j<A->ia[i+1]-1; ++j)
  //     A->r[i] -= A->a[j] * x[A->ja[j]-1];
  // }
  // normb = sqrt( normb );
}

void
ConjugateGradients::dot( const std::vector< tk::real >& a,
                         const std::vector< tk::real >& b,
                         CkCallback c )
// *****************************************************************************
//  Initiate computation of dot product of two vectors
//! \param[in] a 1st vector of dot product
//! \param[in] b 2nd vector of dot product
// *****************************************************************************
{
  Assert( a.size() == b.size(), "Size mismatch" );
  auto d = std::inner_product( begin(a), end(a), begin(b), 0.0 );
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

#include "NoWarning/conjugategradients.def.h"
