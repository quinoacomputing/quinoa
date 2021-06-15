// *****************************************************************************
/*!
  \file      src/LinearSolver/CSR.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Compressed sparse row (CSR) storage for a sparse matrix
  \details   Compressed sparse row (CSR) storage for a sparse matrix.
*/
// *****************************************************************************

#include "Exception.hpp"
#include "CSR.hpp"

using tk::CSR;

CSR::CSR( std::size_t nc,
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& psup )
try :
  ncomp( nc ),
  rnz( psup.second.size()-1 ),
  ia( rnz.size()*ncomp+1 )
// *****************************************************************************
//  Constructor: Create a CSR symmetric matrix with ncomp scalar components per
//  non-zero matrix entry, storing only the upper triangular part
//! \param[in] ncomp Number of scalar components (degrees of freedom)
//! \param[in] psup Points surrounding points of mesh graph, see tk::genPsup
// *****************************************************************************
{
  Assert( ncomp > 0, "Sparse matrix ncomp must be positive" );
  Assert( rnz.size() > 0, "Sparse matrix size must be positive" );

  const auto& psup1 = psup.first;
  const auto& psup2 = psup.second;

  // Calculate number of nonzeros in each block row (rnz), total number of
  // nonzeros (nnz), and fill in row indices (ia)
  std::size_t nnz, i;
  for (ia[0]=1, nnz=i=0; i<psup2.size()-1; ++i) {
    // add up and store nonzeros of row i (only upper triangular part)
    std::size_t j;
    for (rnz[i]=1, j=psup2[i]+1; j<=psup2[i+1]; ++j)
      ++rnz[i];

    // add up total number of nonzeros
    nnz += rnz[i] * ncomp;

    // fill up row index
    for (std::size_t k=0; k<ncomp; ++k)
      ia[i*ncomp+k+1] = ia[i*ncomp+k] + rnz[i];
  }

  // Allocate storage for matrix values and column indices
  a.resize( nnz, 0.0 );
  ja.resize( nnz );

  // fill column indices
  for (i=0; i<rnz.size(); ++i)
    for (std::size_t k=0; k<ncomp; ++k) {
      auto itmp = i*ncomp+k;
      ja[ia[itmp]-1] = itmp+1;  // put in column index of diagonal
      for (std::size_t n=1, j=psup2[i]+1; j<=psup2[i+1]; ++j) {
        // put in column index of an off-diagonal
	ja[ia[itmp]-1+(n++)] = psup1[j]*ncomp+k+1;
      }
    }

  // (bubble-)sort column indices
  for (i=0; i<rnz.size(); ++i)
    for (std::size_t k=0; k<ncomp; ++k)
      for (std::size_t j=psup2[i]+1; j<=psup2[i+1]; ++j)
         for (std::size_t l=1; l<rnz[i]; ++l)   // sort column indices of row i
            for (std::size_t e=0; e<rnz[i]-l; ++e)
              if (ja[ia[i*ncomp+k]-1+e] > ja[ia[i*ncomp+k]+e])
	        std::swap( ja[ia[i*ncomp+k]-1+e], ja[ia[i*ncomp+k]+e] );

} // Catch std::exception
  catch (std::exception& se) {
    // (re-)throw tk::Excpetion
    Throw( std::string("RUNTIME ERROR in CSR constructor: ") + se.what() );
  }

const tk::real&
CSR::operator()( std::size_t row, std::size_t col, std::size_t pos ) const
// *****************************************************************************
//  Return const reference to sparse matrix entry at a position specified
//  using relative addressing
//! \param[in] row Block row
//! \param[in] col Block column
//! \param[in] pos Position in block
//! \return Const reference to matrix entry at position specified
// *****************************************************************************
{
  auto rncomp = row * ncomp;

  for (std::size_t n=0, j=ia[rncomp+pos]-1; j<ia[rncomp+pos+1]-1; ++j, ++n)
    if (col*ncomp+pos+1 == ja[j])
      return a[ia[rncomp+pos]-1+n];

  Throw("Sparse matrix index not found");
}

void
CSR::dirichlet( std::size_t i,
                const std::vector< std::size_t >& gid,
                const NodeCommMap& nodecommap,
                std::size_t pos )
// *****************************************************************************
//  Set Dirichlet boundary condition at a node
//! \param[in] i Local id at which to set Dirichlet BC
//! \param[in] gid Local->global node id map
//! \param[in] nodecommap Node communication map with global node ids
//! \param[in] pos Position in block
//! \details In parallel there can be multiple contributions to a single node
//!   on the mesh, and correspondingly, a single matrix row can be partially
//!   represented on multiple partitions. Setting a Dirichlet BC entails
//!   zeroing out the row of the matrix and putting 1/N into the diagonal entry,
//!   where N is the number of partitions that contribute to the mesh node
//!   (matrix row). As a result, when the matrix participates in a matrix-vector
//!   product, where the partial contributions across all partitions are
//!   aggregated, the diagonal will contain 1 after the sum across partitions.
//! \note Both gid and nodecommap are optional - unused in serial. If nodecommap
//!   is empty, serial is assumed and gid is unused.
// *****************************************************************************
{
  // Lambda to count the number of contributions to a node at which to set BC
  auto count = [&]( std::size_t g ){
    return 1.0 + std::count_if( nodecommap.cbegin(), nodecommap.cend(),
                   [&](const auto& s) {
                     return s.second.find(g) != s.second.cend(); } ); };

  auto incomp = i * ncomp;
  auto diag = nodecommap.empty() ? 1.0 : 1.0/count(gid[i]);

  for (std::size_t j=ia[incomp+pos]-1; j<ia[incomp+pos+1]-1; ++j)
    if (incomp+pos+1==ja[j]) a[j] = diag; else a[j] = 0.0;
}

void
CSR::mult( const std::vector< real >& x, std::vector< real >& r ) const
// *****************************************************************************
//  Multiply CSR matrix with vector from the right: r = A * x
//! \param[in] x Vector to multiply matrix with from the right
//! \param[in] r Result vector of product r = A * x
//! \note This is only complete in serial. In parallel, this computes the own
//!   contributions to the product, so it must be followed by communication
//!   combining the rows stored on multiple partitions.
// *****************************************************************************
{
  std::fill( begin(r), end(r), 0.0 );

  for (std::size_t i=0; i<rnz.size()*ncomp; ++i)
    for (std::size_t j=ia[i]-1; j<ia[i+1]-1; ++j)
      r[i] += a[j] * x[ja[j]-1];
}

std::ostream&
CSR::write_stored( std::ostream& os ) const
// *****************************************************************************
//  Write out CSR as stored
//! \param[in,out] os Output stream to write to
//! \param[in] csr CSR matrix to write
//! \return Updated output stream
// *****************************************************************************
{
  os << "size (npoin) = " << rnz.size() << '\n';
  os << "ncomp = " << ncomp << '\n';
  os << "rsize (size*ncomp) = " << rnz.size() * ncomp << '\n';
  os << "nnz = " << a.size() << '\n';

  std::size_t i;

  os << "rnz[npoin=" << rnz.size() << "] = { ";
  for (i=0; i<rnz.size()-1; ++i) os << rnz[i] << ", ";
  os << rnz[i] << " }\n";

  os << "ia[rsize+1=" << rnz.size()*ncomp+1 << "] = { ";
  for (i=0; i<ia.size()-1; ++i) os << ia[i] << ", ";
  os << ia[i] << " }\n";

  os << "ja[nnz=" << ja.size() << "] = { ";
  for (i=0; i<ja.size()-1; ++i) os << ja[i] << ", ";
  os << ja[i] << " }\n";

  os << "a[nnz=" << a.size() << "] = { ";
  for (i=0; i<a.size()-1; ++i) os << a[i] << ", ";
  os << a[i] << " }\n";

  return os;
}

std::ostream&
CSR::write_structure( std::ostream& os ) const
// *****************************************************************************
//  Write out CSR nonzero structure
//! \param[in,out] os Output stream to write to
//! \param[in] csr CSR matrix to write
//! \return Updated output stream
// *****************************************************************************
{
  for (std::size_t i=0; i<rnz.size()*ncomp; ++i) {
    // leading zeros
    for (std::size_t j=1; j<ja[ia[i]-1]; ++j) os << ". ";
    for (std::size_t n=ia[i]-1; n<ia[i+1]-1; ++n) {
      // zeros between nonzeros
      if (n>ia[i]-1) for (std::size_t j=ja[n-1]; j<ja[n]-1; ++j) os << ". ";
      // nonzero
      os << "o ";
    }
    // trailing zeros
    for (std::size_t j=ja[ia[i+1]-2]; j<rnz.size()*ncomp; ++j) os << ". ";
    os << '\n';
  }

  return os;
}

std::ostream&
CSR::write_matrix( std::ostream& os ) const
// *****************************************************************************
//  Write out CSR as a real matrix
//! \param[in,out] os Output stream to write to
//! \param[in] csr CSR matrix to write
//! \return Updated output stream
// *****************************************************************************
{
  for (std::size_t i=0; i<rnz.size()*ncomp; ++i) {
    for (std::size_t j=1; j<ja[ia[i]-1]; ++j) os << "0\t";
    for (std::size_t n=ia[i]-1; n<ia[i+1]-1; ++n) {
      if (n>ia[i]-1) for (std::size_t j=ja[n-1]; j<ja[n]-1; ++j ) os << "0\t";
      os << a[n] << '\t';
    }
    for (std::size_t j=ja[ia[i+1]-2]; j<rnz.size()*ncomp; ++j) os << "0\t";
    os << '\n';
  }

  return os;
}

std::ostream&
CSR::write_matlab( std::ostream& os ) const
// *****************************************************************************
//  Write out CSR in Matlab/Octave format
//! \param[in,out] os Output stream to write to
//! \param[in] csr CSR matrix to write
//! \return Updated output stream
// *****************************************************************************
{
  os << "A = [ ";
  for (std::size_t i=0; i<rnz.size()*ncomp; ++i) {
    for (std::size_t j=1; j<ja[ia[i]-1]; ++j) os << "0 ";
    for ( std::size_t n=ia[i]-1; n<ia[i+1]-1; ++n) {
      if (n > ia[i]-1)
        for (std::size_t j=ja[n-1]; j<ja[n]-1; ++j) os << "0 ";
      os << a[n] << ' ';
    }
    for (std::size_t j=ja[ia[i+1]-2]; j<rnz.size()*ncomp; ++j) os << "0 ";
    os << ";\n";
  }
  os << "]\n";

  return os;
}
