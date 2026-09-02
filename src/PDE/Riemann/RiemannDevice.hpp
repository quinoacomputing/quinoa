// *****************************************************************************
/*!
  \file      src/PDE/Riemann/RiemannDevice.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Device-callable helpers for the const-order (_constP) surface path
  \details   Kokkos::Array-based mirrors of the tensor and state helpers that
             HLLCMultiMatConstP::flux needs. The host originals in
             Base/Vector.hpp and MultiMat/MiscMultiMatFns.cpp are left
             untouched, so the p-adaptive path is bit-for-bit unaffected.

             Every function here is a straight transcription of its host
             counterpart with the same operation ordering, so results are
             bitwise identical. Deviations from the host source are limited to:
               - std::array          -> Kokkos::Array
               - std::abs            -> fabs   (std::abs is constexpr in
                                                libstdc++; nvcc rejects it)
               - g_inputdeck lookups -> explicit parameters (the input deck is
                                        not reachable from device code)
               - stressCmp[i][j]     -> stressCmpDev(i,j), closed form of the
                                        same table
*/
// *****************************************************************************
#ifndef RiemannDevice_h
#define RiemannDevice_h

#include <Kokkos_Core.hpp>

#include "Types.hpp"

namespace tk {

//! 3-vector in device-callable storage
using Vec3Dev = Kokkos::Array< tk::real, 3 >;
//! 3x3 tensor in device-callable storage
using Mat3Dev = Kokkos::Array< Kokkos::Array< tk::real, 3 >, 3 >;

//! Zero-initialised 3x3 tensor
KOKKOS_INLINE_FUNCTION Mat3Dev zeroMat3Dev()
{
  Mat3Dev m;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) m[i][j] = 0.0;
  return m;
}

//! Test whether a 3x3 tensor is (numerically) empty
//! \details Mirror of tk::is_matrix_empty in Base/Vector.hpp
KOKKOS_INLINE_FUNCTION bool is_matrix_emptyDev( const Mat3Dev& mat )
{
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      if (fabs(mat[i][j]) > 1.0e-16)
        return false;
  return true;
}

//! Matrix-vector product m*v
//! \details Mirror of tk::matvec in Base/Vector.hpp
KOKKOS_INLINE_FUNCTION Vec3Dev matvecDev( const Mat3Dev& m, const Vec3Dev& v )
{
  Vec3Dev mv{ 0, 0, 0 };
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      mv[i] += m[i][j]*v[j];
  }
  return mv;
}

//! Transpose a 3x3 tensor
//! \details Mirror of tk::transpose3by3 in Base/Vector.hpp
KOKKOS_INLINE_FUNCTION Mat3Dev transpose3by3Dev( const Mat3Dev& mat )
{
  Mat3Dev transMat;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      transMat[i][j] = mat[j][i];
  return transMat;
}

//! Rotation matrix given a normal direction
//! \details Mirror of tk::getRotationMatrix in Base/Vector.hpp. Expression
//!   ordering preserved exactly; only std::abs -> fabs.
KOKKOS_INLINE_FUNCTION Mat3Dev getRotationMatrixDev( const Vec3Dev& r )
{
  Mat3Dev rotMat;
  tk::real rx = r[0];
  tk::real ry = r[1];
  tk::real rz = r[2];
  if (fabs(ry+rz) <= fabs(ry-rz))
  {
    tk::real norm = std::sqrt(2*(1-rx*ry-rx*rz-ry*rz));
    rotMat[0][0] = rx;
    rotMat[0][1] = ry;
    rotMat[0][2] = rz;
    rotMat[1][0] = (ry-rz)/norm;
    rotMat[1][1] = (rz-rx)/norm;
    rotMat[1][2] = (rx-ry)/norm;
    rotMat[2][0] = (rx*(ry+rz)-ry*ry-rz*rz)/norm;
    rotMat[2][1] = (ry*(rx+rz)-rx*rx-rz*rz)/norm;
    rotMat[2][2] = (rz*(rx+ry)-rx*rx-ry*ry)/norm;
  }
  else
  {
    tk::real norm = std::sqrt(2*(1+rz*(ry-rx)+rx*ry));
    rotMat[0][0] = rx;
    rotMat[0][1] = ry;
    rotMat[0][2] = rz;
    rotMat[1][0] = (ry+rz)/norm;
    rotMat[1][1] = (rz-rx)/norm;
    rotMat[1][2] = (-rx-ry)/norm;
    rotMat[2][0] = (rx*(rz-ry)-ry*ry-rz*rz)/norm;
    rotMat[2][1] = (ry*(rx+rz)+rx*rx+rz*rz)/norm;
    rotMat[2][2] = (rz*(rx-ry)-rx*rx-ry*ry)/norm;
  }
  return rotMat;
}

//! Multiply two 3x3 matrices in flat row-major storage: C = A*B
KOKKOS_INLINE_FUNCTION
void mat3mulDev( const tk::real A[9], const tk::real B[9], tk::real C[9] )
{
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) {
      tk::real s = 0.0;
      for (std::size_t k=0; k<3; ++k) s += A[i*3+k] * B[k*3+j];
      C[i*3+j] = s;
    }
}

//! Rotate a second order tensor: returns R*mat*R^T
//! \details Mirror of tk::rotateTensorFast in Riemann/TensorFast.hpp, which is
//!   itself a BLAS-free mirror of tk::rotateTensor in Base/Vector.hpp.
KOKKOS_INLINE_FUNCTION
Mat3Dev rotateTensorDev( const Mat3Dev& mat, const Vec3Dev& r )
{
  if (is_matrix_emptyDev(mat)) return zeroMat3Dev();

  Mat3Dev rotMatrix = getRotationMatrixDev(r);

  tk::real M[9], R[9], Rt[9], T[9], O[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) {
      M[i*3+j] = mat[i][j];
      R[i*3+j] = rotMatrix[i][j];
      Rt[j*3+i] = rotMatrix[i][j];
    }

  mat3mulDev( M, Rt, T );   // T = mat*R^T
  mat3mulDev( R, T, O );    // O = R*T

  Mat3Dev out;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) out[i][j] = O[i*3+j];
  return out;
}

//! Un-rotate a second order tensor: returns R^T*mat*R
//! \details Mirror of tk::unrotateTensorFast in Riemann/TensorFast.hpp.
KOKKOS_INLINE_FUNCTION
Mat3Dev unrotateTensorDev( const Mat3Dev& mat, const Vec3Dev& r )
{
  if (is_matrix_emptyDev(mat)) return zeroMat3Dev();

  Mat3Dev rotMatrix = getRotationMatrixDev(r);

  tk::real M[9], R[9], Rt[9], T[9], O[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) {
      M[i*3+j] = mat[i][j];
      R[i*3+j] = rotMatrix[i][j];
      Rt[j*3+i] = rotMatrix[i][j];
    }

  mat3mulDev( M, R, T );    // T = mat*R
  mat3mulDev( Rt, T, O );   // O = R^T*T

  Mat3Dev out;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) out[i][j] = O[i*3+j];
  return out;
}

//! Reflect a second order tensor: returns R^T*mat*R
//! \details Mirror of tk::reflectTensor in Base/Vector.hpp, which issues two
//!   cblas_dgemm calls with m=n=k=3. Same operation ordering, so results are
//!   bitwise identical. Note this takes the reflection matrix directly, not a
//!   normal vector, matching the host signature: the symmetry boundary state
//!   function builds reflectMat itself as I - 2*n_i*n_j.
KOKKOS_INLINE_FUNCTION
Mat3Dev reflectTensorDev( const Mat3Dev& mat, const Mat3Dev& reflectMat )
{
  if (is_matrix_emptyDev(mat)) return zeroMat3Dev();

  tk::real M[9], R[9], Rt[9], T[9], O[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) {
      M[i*3+j] = mat[i][j];
      R[i*3+j] = reflectMat[i][j];
      Rt[j*3+i] = reflectMat[i][j];
    }

  mat3mulDev( M, R, T );    // T = mat*R
  mat3mulDev( Rt, T, O );   // O = R^T*T

  Mat3Dev out;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j) out[i][j] = O[i*3+j];
  return out;
}

//! Rotate a vector: returns R*v
//! \details Mirror of tk::rotateVector in Base/Vector.hpp
KOKKOS_INLINE_FUNCTION Vec3Dev rotateVectorDev( const Vec3Dev& v,
                                                const Vec3Dev& r )
{
  Mat3Dev rotMat = getRotationMatrixDev(r);
  return matvecDev( rotMat, v );
}

//! Un-rotate a vector: returns R^T*v
//! \details Mirror of tk::unrotateVector in Base/Vector.hpp
KOKKOS_INLINE_FUNCTION Vec3Dev unrotateVectorDev( const Vec3Dev& v,
                                                  const Vec3Dev& r )
{
  Mat3Dev rotMat = getRotationMatrixDev(r);
  return matvecDev( transpose3by3Dev(rotMat), v );
}

} // tk::

namespace inciter {

//! Index of a Cauchy stress component in the 6-component packed storage
//! \details Closed form of the stressCmp table in MultiMatIndexing.hpp:
//!   {{0,3,4},{3,1,5},{4,5,2}}. Diagonal gives i, off-diagonal gives 2+i+j.
KOKKOS_INLINE_FUNCTION std::size_t stressCmpDev( std::size_t i, std::size_t j )
{
  return (i == j) ? i : 2+i+j;
}

//! Inverse deformation gradient tensor for a material
//! \details Mirror of inciter::getDeformGrad in MiscMultiMatFns.cpp. solidx is
//!   passed in rather than read from g_inputdeck.
template< class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION tk::Mat3Dev
getDeformGradDev( std::size_t nmat,
                  std::size_t k,
                  const SolidxT& solidx,
                  const StateT& state )
{
  tk::Mat3Dev gk = tk::zeroMat3Dev();

  if (solidx[k] > 0) {
    // deformation gradient for solids
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j)
        gk[i][j] = state[deformIdx(nmat,solidx[k],i,j)];
    }
  }
  // zero tensor for fluids

  return gk;
}

//! Elastic Cauchy stress tensor for a material
//! \details Mirror of inciter::getCauchyStress in MiscMultiMatFns.cpp. solidx
//!   is passed in rather than read from g_inputdeck.
template< class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION tk::Mat3Dev
getCauchyStressDev( std::size_t nmat,
                    std::size_t k,
                    std::size_t ncomp,
                    const SolidxT& solidx,
                    const StateT& state )
{
  tk::Mat3Dev asigk = tk::zeroMat3Dev();

  // elastic Cauchy stress for solids
  if (solidx[k] > 0) {
    for (std::size_t i=0; i<3; ++i) {
      for (std::size_t j=0; j<3; ++j)
        asigk[i][j] = state[ncomp +
          stressIdx(nmat, solidx[k], stressCmpDev(i,j))];
    }
  }

  return asigk;
}

} // inciter::

#endif // RiemannDevice_h
