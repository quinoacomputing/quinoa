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
          const std::pair< std::vector< std::size_t >,
                           std::vector< std::size_t > >& psup )
try :
  dof( DOF ),
  rnz( psup.second.size()-1 ),
  ia( rnz.size()*DOF+1 )
// *****************************************************************************
//  Constructor: Create a CSR symmetric matrix with DOF degrees of freedom,
//  storing only the upper triangular part.
//! \param[in] DOF Number of scalar components (degrees of freedom)
//! \param[in] psup Points surrounding points of mesh graph, see tk::genPsup
// *****************************************************************************
{
  Assert( dof > 0, "Sparse matrix DOF must be positive" );
  Assert( rnz.size() > 0, "Sparse matrix size must be positive" );

  auto& psup1 = psup.first;
  auto& psup2 = psup.second;

  // Calculate number of nonzeros in each block row (rnz), total number of
  // nonzeros (nnz), and fill in row indices (ia)
  std::size_t nnz, i;
  for (ia[0]=1, nnz=i=0; i<psup2.size()-1; ++i) {
    // add up and store nonzeros of row i (only upper triangular part)
    std::size_t j;
    for (rnz[i]=1, j=psup2[i]+1; j<=psup2[i+1]; ++j)
      ++rnz[i];

    // add up total number of nonzeros
    nnz += rnz[i] * DOF;

    // fill up row index
    for (std::size_t k=0; k<DOF; ++k)
      ia[i*DOF+k+1] = ia[i*DOF+k] + rnz[i];
  }

  // Allocate storage for matrix values and column indices
  a.resize( nnz, 0.0 );
  ja.resize( nnz );

  // fill column indices
  for (i=0; i<rnz.size(); ++i)
    for (std::size_t k=0; k<dof; ++k) {
      auto itmp = i*dof+k;
      ja[ia[itmp]-1] = itmp+1;  // put in column index of main diagonal
      for (std::size_t n=1, j=psup2[i]+1; j<=psup2[i+1]; ++j) {
        // put in column index of an off-diagonal
	ja[ia[itmp]-1+(n++)] = psup1[j]*DOF+k+1;
      }
    }

  // (bubble-)sort column indices
  for (i=0; i<rnz.size(); ++i)
    for (std::size_t k=0; k<dof; ++k)
      for (std::size_t j=psup2[i]+1; j<=psup2[i+1]; ++j)
         for (std::size_t l=1; l<rnz[i]; ++l)   // sort column indices of row i
            for (std::size_t e=0; e<rnz[i]-l; ++e)
              if (ja[ia[i*dof+k]-1+e] > ja[ia[i*dof+k]+e])
	        std::swap( ja[ia[i*dof+k]-1+e], ja[ia[i*dof+k]+e] );

} // Catch std::exception
  catch (std::exception& se) {
    // (re-)throw tk::Excpetion
    Throw( std::string("RUNTIME ERROR in CSR constructor: ") + se.what() );
  }

const tk::real&
CSR::operator()( std::size_t row, std::size_t col, std::size_t pos ) const
// *****************************************************************************
//  Return non-const reference to sparse matrix entry at a position specified
//  using relative addressing
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

std::ostream&
CSR::write_as_stored( std::ostream& os ) const
// *****************************************************************************
//  Write out CSR as stored
//! \param[in,out] os Output stream to write to
//! \param[in] csr CSR matrix to write
//! \return Updated output stream
// *****************************************************************************
{
  os << "size (npoin) = " << rnz.size() << '\n';
  os << "dof = " << dof << '\n';;
  os << "rsize (size*dof) = " << rnz.size() * dof << '\n';
  os << "nnz = " << a.size() << '\n';

  std::size_t i;

  os << "rnz[npoin=" << rnz.size() << "] = { ";
  for (i=0; i<rnz.size()-1; ++i) os << rnz[i] << ", ";
  os << rnz[i] << " }\n";

  os << "ia[rsize+1=" << rnz.size()*dof+1 << "] = { ";
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
CSR::write_as_structure( std::ostream& os ) const
// *****************************************************************************
//  Write out CSR nonzero structure
//! \param[in,out] os Output stream to write to
//! \param[in] csr CSR matrix to write
//! \return Updated output stream
// *****************************************************************************
{
  for (std::size_t i=0; i<rnz.size()*dof; ++i) {
    // leading zeros
    for (std::size_t j=1; j<ja[ia[i]-1]; ++j) os << ". ";
    for (std::size_t n=ia[i]-1; n<ia[i+1]-1; ++n) {
      // zeros between nonzeros
      if (n>ia[i]-1) for (std::size_t j=ja[n-1]; j<ja[n]-1; ++j) os << ". ";
      // nonzero
      os << "o ";
    }
     // trailing zeros
    for (std::size_t j=ja[ia[i+1]-2]; j<rnz.size()*dof; ++j) os << ". ";
    os << '\n';
  }

  return os;
}

std::ostream&
CSR::write_as_matrix( std::ostream& os ) const
// *****************************************************************************
//  Write out CSR as a real matrix
//! \param[in,out] os Output stream to write to
//! \param[in] csr CSR matrix to write
//! \return Updated output stream
// *****************************************************************************
{
  for (std::size_t i=0; i<rnz.size()*dof; ++i) {
    for (std::size_t j=1; j<ja[ia[i]-1]; ++j) os << "0\t";
    for (std::size_t n=ia[i]-1; n<ia[i+1]-1; ++n) {
      if (n>ia[i]-1) for (std::size_t j=ja[n-1]; j<ja[n]-1; ++j ) os << "0\t";
      os << a[n] << '\t';
    }
    for (std::size_t j=ja[ia[i+1]-2]; j<rnz.size()*dof; ++j) os << "0\t";
    os << '\n';
  }

  return os;
}

std::ostream&
CSR::write_as_matlab( std::ostream& os ) const
// *****************************************************************************
//  Write out CSR in Matlab/Octave format
//! \param[in,out] os Output stream to write to
//! \param[in] csr CSR matrix to write
//! \return Updated output stream
// *****************************************************************************
{
  os << "A = [ ";
  for (std::size_t i=0; i<rnz.size(); ++i) {
    for (std::size_t j=1; j<ja[ia[i]-1]; ++j) os << "0 ";
    for ( std::size_t n=ia[i]-1; n<ia[i+1]-1; ++n) {
      if (n > ia[i]-1)
        for (std::size_t j=ja[n-1]; j<ja[n]-1; ++j) os << "0 ";
      os << a[n] << ' ';
    }
    for (std::size_t j=ja[ia[i+1]-2]; j<ia.size()-1; ++j) os << "0 ";
    os << ";\n";
  }
  os << "]\n";

  return os;
}