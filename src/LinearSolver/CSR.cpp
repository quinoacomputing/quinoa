// *****************************************************************************
/*!
  \file      src/LinearSolver/CSR.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Compressed sparse row (CSR) storage for a sparse matrix
  \details   Compressed sparse row (CSR) storage for a sparse matrix.
*/
// *****************************************************************************

#include "CSR.hpp"
#include "Exception.hpp"

using tk::CSR;

CSR::CSR( std::size_t DOF,
          std::size_t size,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& psup ) :
  dof( DOF ),
  rnz( size ),
  ia( size*DOF+1 )
// *****************************************************************************
//  Constructor: Create a CSR matrix for a size x size sparse symmetric matrix
//  with DOF degrees of freedom, storing only the upper triangular part.
//! \param[in] DOF Number of scalar components (degrees of freedom)
//! \param[in] size Number of scalar components (degrees of freedom)
//! \param[in] psup Points surrounding points of mesh graph, see tk::genPsup
// *****************************************************************************
{
  // Calculate number of nonzeros in each block row (rnz), total number of
  // nonzeros (nnz), and fill in row indices (ia)
  std::size_t nnz, i, j;
  for (ia[0]=1, nnz=i=0; i<size; ++i) {
    // add up and store nonzeros of row i (only upper triangular part)
    for (rnz[i]=1, j=psup.second[i]+1; j<=psup.second[i+1]; ++j)
      ++rnz[i];

    // add up total number of nonzeros
    nnz += rnz[i] * DOF;

    // fill up row index
    for (std::size_t k=0; k<DOF; ++k)
      ia[i*DOF+k+1] = ia[i*DOF+k] + rnz[i];
  }

  // Allocate storage for matrix values and column indices
  a.resize( nnz );
  ia.resize( nnz );
}

const tk::real&
CSR::operator()( std::size_t row, std::size_t col, std::size_t pos ) const
// *****************************************************************************
// Return non-const reference to sparse matrix entry at a position specified
// using relative addressing
//! \param[in] row Block row
//! \param[in] col Block column
//! \param[in] pos Position in block
//! \return Const reference to matrix entry at position specified
// *****************************************************************************
{
  auto rdof = row * dof;

  for (std::size_t n=0, j=ia[rdof+pos]-1; j<=ia[rdof+pos+1]-2; ++j, ++n)
    if (col*dof+pos+1 == ja[j])
      return a[ia[rdof+pos]-1+n];

  Throw("Sparse matrix index not found");
}
