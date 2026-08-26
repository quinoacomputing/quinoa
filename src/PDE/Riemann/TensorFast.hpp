// *****************************************************************************
/*!
  \file      src/PDE/Riemann/TensorFast.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     BLAS-free 3x3 tensor rotations for the const-order (_constP) path
  \details   tk::rotateTensor and tk::unrotateTensor in Base/Vector.hpp compute
             R*mat*R^T and R^T*mat*R by issuing two cblas_dgemm calls with
             m=n=k=3. The BLAS dispatch, packing and internal locking cost two
             to three orders of magnitude more than the 27 multiply-adds it
             performs. Profiling a DGP2 water-air case put ~25% of total runtime
             in OpenBLAS dgemm reached from these two functions.

             These are drop-in replacements with identical semantics, for use
             by the _constP path only. The originals are left untouched so the
             p-adaptive path is bit-for-bit unaffected.
*/
// *****************************************************************************
#ifndef TensorFast_h
#define TensorFast_h

#include "Vector.hpp"

namespace tk {

//! Multiply two 3x3 matrices: C = A*B
inline void mat3mul( const tk::real A[9], const tk::real B[9], tk::real C[9] )
{
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) {
      tk::real s = 0.0;
      for (std::size_t k=0; k<3; ++k) s += A[i*3+k] * B[k*3+j];
      C[i*3+j] = s;
    }
}

//! Rotate a second order tensor: returns R*mat*R^T
//! \details Semantically identical to tk::rotateTensor, without BLAS.
inline std::array< std::array< tk::real, 3 >, 3 >
rotateTensorFast( const std::array< std::array< tk::real, 3 >, 3 >& mat,
                  const std::array< tk::real, 3 >& r )
{
  if (is_matrix_empty(mat)) return {{}};

  std::array< std::array< tk::real, 3 >, 3 > rotMatrix = getRotationMatrix(r);

  tk::real M[9], R[9], Rt[9], T[9], O[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) {
      M[i*3+j] = mat[i][j];
      R[i*3+j] = rotMatrix[i][j];
      Rt[j*3+i] = rotMatrix[i][j];
    }

  mat3mul( M, Rt, T );   // T = mat*R^T
  mat3mul( R, T, O );    // O = R*T

  return {{ {O[0], O[1], O[2]}, {O[3], O[4], O[5]}, {O[6], O[7], O[8]} }};
}

//! Un-rotate a second order tensor: returns R^T*mat*R
//! \details Semantically identical to tk::unrotateTensor, without BLAS.
inline std::array< std::array< tk::real, 3 >, 3 >
unrotateTensorFast( const std::array< std::array< tk::real, 3 >, 3 >& mat,
                    const std::array< tk::real, 3 >& r )
{
  if (is_matrix_empty(mat)) return {{}};

  std::array< std::array< tk::real, 3 >, 3 > rotMatrix = getRotationMatrix(r);

  tk::real M[9], R[9], Rt[9], T[9], O[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) {
      M[i*3+j] = mat[i][j];
      R[i*3+j] = rotMatrix[i][j];
      Rt[j*3+i] = rotMatrix[i][j];
    }

  mat3mul( M, R, T );    // T = mat*R
  mat3mul( Rt, T, O );   // O = R^T*T

  return {{ {O[0], O[1], O[2]}, {O[3], O[4], O[5]}, {O[6], O[7], O[8]} }};
}

} // tk::

#endif // TensorFast_h
