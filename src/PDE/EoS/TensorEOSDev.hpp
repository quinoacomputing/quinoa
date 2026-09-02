// *****************************************************************************
/*!
  \file      src/PDE/EoS/TensorEOSDev.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     BLAS/LAPACK-free 3x3 strain-tensor helpers for the EoS device path
  \details   Device-callable mirrors of the strain tensor helpers in
             Base/Vector.hpp that the elastic energy terms need:

               tk::getRightCauchyGreen        -> rightCauchyGreenEOS
               tk::getIsochorRightCauchyGreen -> isochorRightCauchyGreenEOS
               tk::getDevHencky               -> devHenckyEOS

             The host originals are left untouched, so the p-adaptive path is
             bit-for-bit unaffected.

             Two constraints shape this file:

             1. It cannot use Kokkos. EOS.hpp is included by host-only targets
                compiled with mpicxx, where Kokkos' __CUDACC__ guard rejects
                Kokkos_Core.hpp. Hence EOS_FN rather than
                KOKKOS_INLINE_FUNCTION, and raw C arrays rather than
                Kokkos::Array. This is why these are not simply reused from
                Riemann/RiemannDevice.hpp, which is Kokkos::Array-based.

             2. It cannot use std::array. std::array::operator[] is constexpr
                in libstdc++ and nvcc rejects calls to it from device code
                without --expt-relaxed-constexpr, which this build does not
                pass. The existing device soundspeed bodies sidestep this by
                accepting a std::array they never index; the elastic terms
                must index it, so the interface here is tk::real[3][3].
                Same reason fmax/fabs are used instead of std::max/std::abs.

  \warning   NOT bit-identical to the host originals. tk::getRightCauchyGreen
             inverts g.g^T with LAPACKE_dgetrf/dgetri (LU with partial
             pivoting); inv3x3EOS below uses a closed-form adjugate/
             determinant inverse. Different algorithm, therefore different
             rounding. Measured over 2e5 random inverse deformation gradients
             per regime, perturbations of +/-1e-6 to +/-1e-1 about identity
             (cond(g.g^T) <= 2.7):

               absolute deviation from the host result
                 inv(g.g^T)  <= 8.9e-16      devHencky  <= 2.2e-16
               error against an 80-bit reference, in units of double eps
                 LAPACK    avg 0.44 ulp, max 1.40 ulp
                 adjugate  avg 0.65 ulp, max 2.86 ulp

             So the adjugate carries roughly 1.5x LAPACK's rounding error,
             both at the few-ulp level. For ill-conditioned inputs
             (cond ~1e9) both lose all accuracy and the adjugate is ~15%
             worse; the host path is already unusable there.
*/
// *****************************************************************************
#ifndef TensorEOSDev_h
#define TensorEOSDev_h

#include <cmath>

#include "Types.hpp"
#include "EoS/EOSDeviceFn.hpp"

namespace tk {

//! Test whether a 3x3 tensor is (numerically) empty
//! \details Mirror of tk::is_matrix_empty in Base/Vector.hpp
EOS_FN bool isMatrixEmptyEOS( const tk::real mat[3][3] )
{
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      if (fabs(mat[i][j]) > 1.0e-16) return false;
  return true;
}

//! Determinant of a 3x3 tensor
//! \details Mirror of tk::determinant in Base/Vector.hpp, same term ordering
EOS_FN tk::real determinantEOS( const tk::real a[3][3] )
{
  return ( a[0][0] * (a[1][1]*a[2][2]-a[1][2]*a[2][1])
         - a[0][1] * (a[1][0]*a[2][2]-a[1][2]*a[2][0])
         + a[0][2] * (a[1][0]*a[2][1]-a[1][1]*a[2][0]) );
}

//! Compute A = g.g^T for a 3x3 tensor
//! \details Replaces the cblas_dgemm(NoTrans,Trans) in
//!   tk::getRightCauchyGreen and tk::getDevHencky. The accumulation order
//!   matches the naive triple loop, as in tk::mat3mul in Riemann/TensorFast.hpp
EOS_FN void selfTransposeProductEOS( const tk::real g[3][3], tk::real A[3][3] )
{
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) {
      tk::real s = 0.0;
      for (std::size_t k=0; k<3; ++k) s += g[i][k] * g[j][k];
      A[i][j] = s;
    }
}

//! Closed-form inverse of a 3x3 tensor via adjugate / determinant
//! \details Replaces LAPACKE_dgetrf + LAPACKE_dgetri. See the file-level
//!   warning: this is NOT bit-identical to the LU inverse it replaces.
EOS_FN void inv3x3EOS( const tk::real A[3][3], tk::real Ainv[3][3] )
{
  const tk::real a=A[0][0], b=A[0][1], c=A[0][2],
                 d=A[1][0], e=A[1][1], f=A[1][2],
                 g=A[2][0], h=A[2][1], i=A[2][2];

  // cofactors (transposed, i.e. the adjugate)
  const tk::real c00 =  (e*i - f*h);
  const tk::real c01 = -(b*i - c*h);
  const tk::real c02 =  (b*f - c*e);
  const tk::real c10 = -(d*i - f*g);
  const tk::real c11 =  (a*i - c*g);
  const tk::real c12 = -(a*f - c*d);
  const tk::real c20 =  (d*h - e*g);
  const tk::real c21 = -(a*h - b*g);
  const tk::real c22 =  (a*e - b*d);

  const tk::real det = a*c00 + b*c10 + c*c20;
  const tk::real invdet = 1.0/det;

  Ainv[0][0] = c00*invdet; Ainv[0][1] = c01*invdet; Ainv[0][2] = c02*invdet;
  Ainv[1][0] = c10*invdet; Ainv[1][1] = c11*invdet; Ainv[1][2] = c12*invdet;
  Ainv[2][0] = c20*invdet; Ainv[2][1] = c21*invdet; Ainv[2][2] = c22*invdet;
}

//! Right Cauchy-Green strain tensor C = inv(g.g^T)
//! \details Mirror of tk::getRightCauchyGreen in Base/Vector.hpp
EOS_FN void rightCauchyGreenEOS( const tk::real g[3][3], tk::real C[3][3] )
{
  // Return empty if input matrix is empty
  if (isMatrixEmptyEOS(g)) {
    for (std::size_t i=0; i<3; ++i)
      for (std::size_t j=0; j<3; ++j) C[i][j] = 0.0;
    return;
  }

  tk::real A[3][3];
  selfTransposeProductEOS( g, A );   // A = g.g^T
  inv3x3EOS( A, C );                 // C = inv(A)
}

//! Volume-preserving part of the right Cauchy-Green strain tensor
//! \details Mirror of tk::getIsochorRightCauchyGreen in Base/Vector.hpp.
//!   Note the host original does NOT short-circuit on an empty input itself;
//!   it relies on getRightCauchyGreen returning zeros, which then makes detC
//!   zero and the division produce nan. That behaviour is preserved here
//!   rather than "fixed", so the two agree on degenerate input.
EOS_FN void isochorRightCauchyGreenEOS( const tk::real g[3][3],
                                        tk::real Ct[3][3] )
{
  rightCauchyGreenEOS( g, Ct );
  auto detC = std::pow(determinantEOS(Ct), 1.0/3.0);
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) Ct[i][j] /= detC;
}

//! Deviatoric Hencky strain tensor
//! \details Mirror of tk::getDevHencky in Base/Vector.hpp. Uses the
//!   approximation of section 3.4 of Barton, P. T. (2019), An
//!   interface-capturing Godunov method for the simulation of compressible
//!   solid-fluid problems. Journal of Computational Physics, 390, 25-50.
//!   Operation ordering follows the host original exactly.
EOS_FN void devHenckyEOS( const tk::real g[3][3], tk::real devH[3][3] )
{
  // Return empty if input matrix is empty
  if (isMatrixEmptyEOS(g)) {
    for (std::size_t i=0; i<3; ++i)
      for (std::size_t j=0; j<3; ++j) devH[i][j] = 0.0;
    return;
  }

  // Get right volume-preserving part of Cauchy-Green strain tensor
  tk::real C[3][3];
  isochorRightCauchyGreenEOS( g, C );

  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) devH[i][j] = 0.0;

  // get g.g^T (i.e. inv(C))
  tk::real CInv[3][3];
  selfTransposeProductEOS( g, CInv );

  // volume-preserving part
  auto detCInv = std::pow(determinantEOS(CInv), 1.0/3.0);
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) CInv[i][j] /= detCInv;

  // Compute (C-CInv)/4
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      devH[i][j] = 0.25*(C[i][j]-CInv[i][j]);

  // Compute deviatoric part
  tk::real trH = devH[0][0] + devH[1][1] + devH[2][2];

  for (std::size_t i=0; i<3; ++i) devH[i][i] -= trH/3.0;
}

} // tk::

#endif // TensorEOSDev_h
