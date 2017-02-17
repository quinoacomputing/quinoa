/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <app/UseCase_14_Common.hpp>

typedef int		Int;
typedef double		Real;


//
//  Defined numerical constants
//
const double ONE12TH = (1.0/12.0);
//
//  Indexes into a 4 by 3 matrix of hourglass constants
//
const int HG_X1 = 0;
const int HG_Y1 = 1;
const int HG_Z1 = 2;
const int HG_X2 = 3;
const int HG_Y2 = 4;
const int HG_Z2 = 5;
const int HG_X3 = 6;
const int HG_Y3 = 7;
const int HG_Z3 = 8;
const int HG_X4 = 9;
const int HG_Y4 = 10;
const int HG_Z4 = 11;
//
//  Indexes into a 3 by 3 symmetric tensor stored as a length 6 vector
//
const int K_S_XX = 0;
const int K_S_YY = 1;
const int K_S_ZZ = 2;
const int K_S_XY = 3;
const int K_S_YZ = 4;
const int K_S_ZX = 5;
const int K_S_YX = 3;
const int K_S_ZY = 4;
const int K_S_XZ = 5;
//
//  Indexes into a full 3 by 3 tensor stored as a length 9 vector
//
const int K_F_XX = 0;
const int K_F_YY = 1;
const int K_F_ZZ = 2;
const int K_F_XY = 3;
const int K_F_YZ = 4;
const int K_F_ZX = 5;
const int K_F_YX = 6;
const int K_F_ZY = 7;
const int K_F_XZ = 8;
//
//  Indexes into a 3 by 3 skew symmetric tensor stored as a length 9 vector
//
const int K_V_XY = 0;
const int K_V_YZ = 1;
const int K_V_ZX = 2;
//
//  Indexes into a 2 by 2 symmetric tensor stored as a length 3 vector
//
const int K_S_21_XX = 0;
const int K_S_21_YY = 1;
const int K_S_21_XY = 2;
const int K_S_21_YX = 2;
//
//  Indexes into a full 2 by 2 tensor stored as a length 3 vector
//
const int K_F_22_XX = 0;
const int K_F_22_YY = 1;
const int K_F_22_XY = 2;
const int K_F_22_YX = 3;

namespace APS {
namespace Hex {

//***************************************************************************************************************
//
//  Additive decomposition for a gradient operator
//
//    Input:
//      gradient(9)     Input velocity gradient
//
//    Output:
//      stretching(6)
//      vorticity(3)
//
//***************************************************************************************************************
TEMPLATE_SCALAR
INLINE void additive_decomp36(const Scalar *const gradient, Scalar *const stretching, Scalar *const vorticity) {
  //
  //  Symmetric part
  //
  stretching[K_S_XX] = gradient[K_F_XX];
  stretching[K_S_YY] = gradient[K_F_YY];
  stretching[K_S_ZZ] = gradient[K_F_ZZ];
  stretching[K_S_XY] = 0.5*(gradient[K_F_XY] + gradient[K_F_YX]);
  stretching[K_S_YZ] = 0.5*(gradient[K_F_YZ] + gradient[K_F_ZY]);
  stretching[K_S_ZX] = 0.5*(gradient[K_F_ZX] + gradient[K_F_XZ]);
  //
  //  Skew Symmetric part
  //
  vorticity[K_V_XY] = 0.5*(gradient[K_F_XY] - gradient[K_F_YX]);
  vorticity[K_V_YZ] = 0.5*(gradient[K_F_YZ] - gradient[K_F_ZY]);
  vorticity[K_V_ZX] = 0.5*(gradient[K_F_ZX] - gradient[K_F_XZ]);
}

//***************************************************************************************************************
//
//  Polar decomposition of a stretching tensor
//
//    Input
//      dt                   Time multiplication factor
//      stretching(6)         Current element stretching tensor
//      vorticity(3)          The skew symmetric portion of the velocity gradient
//      rotation_old(9)       Old element rotation tensor
//    Output
//      stretch(6)            New calculated element stretch
//      rotation_new(9)       New element rotation tensor
//
//***************************************************************************************************************
TEMPLATE_SCALAR
INLINE void polar_decomp33r2(const Scalar dt,
                             const Scalar *const stretching,
                             const Scalar *const vorticity,
                             const Scalar *const rotation_old,
                             Scalar *const stretch,
                             Scalar *const rotation_new) 
{
  //
  // calculate the rates of rotation via gauss elimination.
  //
  Scalar z1 = stretching[K_S_XY] * stretch[K_S_ZX] - stretching[K_S_ZX] * stretch[K_S_XY] + stretching[K_S_YY] * stretch[K_S_YZ] -
    stretching[K_S_YZ] * stretch[K_S_YY] + stretching[K_S_YZ] * stretch[K_S_ZZ] - stretching[K_S_ZZ] * stretch[K_S_YZ];
  Scalar z2 = stretching[K_S_ZX] * stretch[K_S_XX] - stretching[K_S_XX] * stretch[K_S_ZX] + stretching[K_S_YZ] * stretch[K_S_XY] -
    stretching[K_S_XY] * stretch[K_S_YZ] + stretching[K_S_ZZ] * stretch[K_S_ZX] - stretching[K_S_ZX] * stretch[K_S_ZZ];
  Scalar z3 = stretching[K_S_XX] * stretch[K_S_XY] - stretching[K_S_XY] * stretch[K_S_XX] + stretching[K_S_XY] * stretch[K_S_YY] -
    stretching[K_S_YY] * stretch[K_S_XY] + stretching[K_S_ZX] * stretch[K_S_YZ] - stretching[K_S_YZ] * stretch[K_S_ZX];

  //
  // forward elimination
  //
  const Scalar a1inv = 1.0 / (stretch[K_S_YY] + stretch[K_S_ZZ]);

  const Scalar a4BYa1 = -stretch[K_S_XY] * a1inv;
  const Scalar a2inv = 1.0 / (stretch[K_S_ZZ] + stretch[K_S_XX] + stretch[K_S_XY] * a4BYa1);

  const Scalar a5 =  -stretch[K_S_YZ] + stretch[K_S_ZX] * a4BYa1;
  z2 -= z1 * a4BYa1;
  const Scalar a6BYa1 = -stretch[K_S_ZX] * a1inv;
  const Scalar a5BYa2 = a5 * a2inv;
  z3 -= z1 * a6BYa1 - z2 * a5BYa2;
  //
  // backward substitution -
  //
  z3 /= (stretch[K_S_XX] + stretch[K_S_YY] + stretch[K_S_ZX] * a6BYa1 + a5 * a5BYa2);
  z2 = (z2 - a5 * z3) * a2inv;
  z1 = (z1*a1inv - a6BYa1 * z3 -a4BYa1 * z2);

  //
  // calculate rotation rates - recall that spin_rate is an asymmetric tensor,
  // so compute spin rate vector as dual of spin rate tensor,
  // i.e   w_i = e_ijk * spin_rate_jk
  //
  z1 += vorticity[K_V_YZ];
  z2 += vorticity[K_V_ZX];
  z3 += vorticity[K_V_XY];
  //
  // update rotation tensor:
  //   1) premultiply old rotation tensor to get right-hand side.
  //
  const Scalar dt_half = 0.5 * dt;
  Scalar r_XX = rotation_old[K_F_XX] + dt_half*( z3 * rotation_old[K_F_YX] - z2 * rotation_old[K_F_ZX] );
  Scalar r_YX = rotation_old[K_F_YX] + dt_half*( z1 * rotation_old[K_F_ZX] - z3 * rotation_old[K_F_XX] );
  Scalar r_ZX = rotation_old[K_F_ZX] + dt_half*( z2 * rotation_old[K_F_XX] - z1 * rotation_old[K_F_YX] );
  Scalar r_XY = rotation_old[K_F_XY] + dt_half*( z3 * rotation_old[K_F_YY] - z2 * rotation_old[K_F_ZY] );
  Scalar r_YY = rotation_old[K_F_YY] + dt_half*( z1 * rotation_old[K_F_ZY] - z3 * rotation_old[K_F_XY] );
  Scalar r_ZY = rotation_old[K_F_ZY] + dt_half*( z2 * rotation_old[K_F_XY] - z1 * rotation_old[K_F_YY] );
  Scalar r_XZ = rotation_old[K_F_XZ] + dt_half*( z3 * rotation_old[K_F_YZ] - z2 * rotation_old[K_F_ZZ] );
  Scalar r_YZ = rotation_old[K_F_YZ] + dt_half*( z1 * rotation_old[K_F_ZZ] - z3 * rotation_old[K_F_XZ] );
  Scalar r_ZZ = rotation_old[K_F_ZZ] + dt_half*( z2 * rotation_old[K_F_XZ] - z1 * rotation_old[K_F_YZ] );
  //
  //    2) solve for new rotation tensor via gauss elimination.
  // forward elimination -
  //

  Scalar a12 = - dt_half * z3;
  Scalar a13 =   dt_half * z2;
  Scalar b32 = - dt_half * z1;
  Scalar a22inv = 1.0 / (1.0 + a12 * a12);

  Scalar a13a12 = a13*a12;
  Scalar a23 = b32 + a13a12;
  r_YX += r_XX * a12;
  r_YY += r_XY * a12;
  r_YZ += r_XZ * a12;
  b32 = (b32 - a13a12) * a22inv;
  r_ZX += r_XX * a13 + r_YX * b32;
  r_ZY += r_XY * a13 + r_YY * b32;
  r_ZZ += r_XZ * a13 + r_YZ * b32;
  //
  // backward substitution -
  //
  const Scalar a33inv = 1.0 / (1.0 + a13 * a13 + a23 * b32);

  rotation_new[K_F_ZX]  = r_ZX * a33inv;
  rotation_new[K_F_ZY]  = r_ZY * a33inv;
  rotation_new[K_F_ZZ]  = r_ZZ * a33inv;
  rotation_new[K_F_YX]  = ( r_YX - rotation_new[K_F_ZX] * a23 ) * a22inv;
  rotation_new[K_F_YY]  = ( r_YY - rotation_new[K_F_ZY] * a23 ) * a22inv;
  rotation_new[K_F_YZ]  = ( r_YZ - rotation_new[K_F_ZZ] * a23 ) * a22inv;
  rotation_new[K_F_XX]  = r_XX - rotation_new[K_F_ZX] * a13 - rotation_new[K_F_YX] * a12;
  rotation_new[K_F_XY]  = r_XY - rotation_new[K_F_ZY] * a13 - rotation_new[K_F_YY] * a12;
  rotation_new[K_F_XZ]  = r_XZ - rotation_new[K_F_ZZ] * a13 - rotation_new[K_F_YZ] * a12;
  //
  // update stretch tensor in the new configuration -
  //
  const Scalar a1 = stretching[K_S_XY] + vorticity[K_V_XY];
  const Scalar a2 = stretching[K_S_YZ] + vorticity[K_V_YZ];
  const Scalar a3 = stretching[K_S_ZX] + vorticity[K_V_ZX];
  const Scalar b1 = stretching[K_S_ZX] - vorticity[K_V_ZX];
  const Scalar b2 = stretching[K_S_XY] - vorticity[K_V_XY];
  const Scalar b3 = stretching[K_S_YZ] - vorticity[K_V_YZ];

  const Scalar s_XX = stretch[K_S_XX];
  const Scalar s_YY = stretch[K_S_YY];
  const Scalar s_ZZ = stretch[K_S_ZZ];
  const Scalar s_XY = stretch[K_S_XY];
  const Scalar s_YZ = stretch[K_S_YZ];
  const Scalar s_ZX = stretch[K_S_ZX];

  stretch[K_S_XX] += dt * (stretching[K_S_XX] * s_XX + ( a1 + z3 ) * s_XY + ( b1 - z2 ) * s_ZX);
  stretch[K_S_YY] += dt * (stretching[K_S_YY] * s_YY + ( a2 + z1 ) * s_YZ + ( b2 - z3 ) * s_XY);
  stretch[K_S_ZZ] += dt * (stretching[K_S_ZZ] * s_ZZ + ( a3 + z2 ) * s_ZX + ( b3 - z1 ) * s_YZ);
  stretch[K_S_XY] += dt * (stretching[K_S_XX] * s_XY + ( a1 )      * s_YY + ( b1      ) * s_YZ - z3 * s_XX + z1 * s_ZX);
  stretch[K_S_YZ] += dt * (stretching[K_S_YY] * s_YZ + ( a2 )      * s_ZZ + ( b2      ) * s_ZX - z1 * s_YY + z2 * s_XY);
  stretch[K_S_ZX] += dt * (stretching[K_S_ZZ] * s_ZX + ( a3 )      * s_XX + ( b3      ) * s_XY - z2 * s_ZZ + z3 * s_YZ);
}

TEMPLATE_SCALAR
INLINE void polar_decomp33_nr(const Scalar dt,
                              const Scalar *const stretching,
                              const Scalar *const vorticity,
                              Scalar *const stretch) {
  //
  // calculate the rates of rotation via gauss elimination.
  //
  Scalar z1 = stretching[K_S_XY] * stretch[K_S_ZX] - stretching[K_S_ZX] * stretch[K_S_XY] + stretching[K_S_YY] * stretch[K_S_YZ] -
    stretching[K_S_YZ] * stretch[K_S_YY] + stretching[K_S_YZ] * stretch[K_S_ZZ] - stretching[K_S_ZZ] * stretch[K_S_YZ];
  Scalar z2 = stretching[K_S_ZX] * stretch[K_S_XX] - stretching[K_S_XX] * stretch[K_S_ZX] + stretching[K_S_YZ] * stretch[K_S_XY] -
    stretching[K_S_XY] * stretch[K_S_YZ] + stretching[K_S_ZZ] * stretch[K_S_ZX] - stretching[K_S_ZX] * stretch[K_S_ZZ];
  Scalar z3 = stretching[K_S_XX] * stretch[K_S_XY] - stretching[K_S_XY] * stretch[K_S_XX] + stretching[K_S_XY] * stretch[K_S_YY] -
    stretching[K_S_YY] * stretch[K_S_XY] + stretching[K_S_ZX] * stretch[K_S_YZ] - stretching[K_S_YZ] * stretch[K_S_ZX];
  //
  // forward elimination
  //
  Scalar a1inv = 1.0 / (stretch[K_S_YY] + stretch[K_S_ZZ]);
  Scalar a4BYa1 = -stretch[K_S_XY] * a1inv;
  Scalar a2inv = 1.0 / (stretch[K_S_ZZ] + stretch[K_S_XX] + stretch[K_S_XY] * a4BYa1);
  Scalar a5 =  -stretch[K_S_YZ] + stretch[K_S_ZX] * a4BYa1;
  z2 -= z1 * a4BYa1;
  Scalar a6BYa1 = -stretch[K_S_ZX] * a1inv;
  Scalar a5BYa2 = a5 * a2inv;
  z3 -= z1 * a6BYa1 - z2 * a5BYa2;
  //
  // backward substitution -
  //
  z3 /= (stretch[K_S_XX] + stretch[K_S_YY] + stretch[K_S_ZX] * a6BYa1 + a5 * a5BYa2);
  z2 = (z2 - a5 * z3) * a2inv;
  z1 = (z1*a1inv - a6BYa1 * z3 -a4BYa1 * z2);

  //
  // calculate rotation rates - recall that spin_rate is an asymmetric tensor,
  // so compute spin rate vector as dual of spin rate tensor,
  // i.e   w_i = e_ijk * spin_rate_jk
  //
  z1 += vorticity[K_V_YZ];
  z2 += vorticity[K_V_ZX];
  z3 += vorticity[K_V_XY];
  //
  // update rotation tensor:
  //   1) premultiply old rotation tensor to get right-hand side.
  //
  const Scalar dt_half = 0.5 * dt;
  //
  //    2) solve for new rotation tensor via gauss elimination.
  // forward elimination -
  //

  Scalar a12 = - dt_half * z3;
  Scalar a13 =   dt_half * z2;
  Scalar b32 = - dt_half * z1;
  Scalar a22inv = 1.0 / (1.0 + a12 * a12);
  b32 = (b32-a12 * a13) * a22inv;
  //
  // update stretch tensor in the new configuration -
  //
  {
    Scalar a1 = stretching[K_S_XY] + vorticity[K_V_XY];
    Scalar a2 = stretching[K_S_YZ] + vorticity[K_V_YZ];
    Scalar a3 = stretching[K_S_ZX] + vorticity[K_V_ZX];
    Scalar b1 = stretching[K_S_ZX] - vorticity[K_V_ZX];
    Scalar b2 = stretching[K_S_XY] - vorticity[K_V_XY];
    Scalar b3 = stretching[K_S_YZ] - vorticity[K_V_YZ];

    Scalar s_XX = stretch[K_S_XX];
    Scalar s_YY = stretch[K_S_YY];
    Scalar s_ZZ = stretch[K_S_ZZ];
    Scalar s_XY = stretch[K_S_XY];
    Scalar s_YZ = stretch[K_S_YZ];
    Scalar s_ZX = stretch[K_S_ZX];

    stretch[K_S_XX] += dt * (stretching[K_S_XX] * s_XX + ( a1 + z3 ) * s_XY + ( b1 - z2 ) * s_ZX);
    stretch[K_S_YY] += dt * (stretching[K_S_YY] * s_YY + ( a2 + z1 ) * s_YZ + ( b2 - z3 ) * s_XY);
    stretch[K_S_ZZ] += dt * (stretching[K_S_ZZ] * s_ZZ + ( a3 + z2 ) * s_ZX + ( b3 - z1 ) * s_YZ);
    stretch[K_S_XY] += dt * (stretching[K_S_XX] * s_XY + ( a1 )      * s_YY + ( b1      ) * s_YZ - z3 * s_XX + z1 * s_ZX);
    stretch[K_S_YZ] += dt * (stretching[K_S_YY] * s_YZ + ( a2 )      * s_ZZ + ( b2      ) * s_ZX - z1 * s_YY + z2 * s_XY);
    stretch[K_S_ZX] += dt * (stretching[K_S_ZZ] * s_ZX + ( a3 )      * s_XX + ( b3      ) * s_XY - z2 * s_ZZ + z3 * s_YZ);
  }
}

//***************************************************************************************************************
//
//  Use a rotation matrix to rotate a given tensor forward
//
//    Input:
//      rotation(9)           Rotation tensor
//      orig_tensor(6)        Input symmetric tensor
//
//   Output:
//      rot_tensor(6)         Output symmetric tensor
//
//***************************************************************************************************************
TEMPLATE_SCALAR
INLINE void rotate_tensor_forward(const Scalar *const rotation,
				  const Scalar *const orig_tensor,
				  Scalar *const rot_tensor) {

  Scalar t_const[9];

  t_const[0] = orig_tensor[K_S_XX]*rotation[K_F_XX] + orig_tensor[K_S_XY]*rotation[K_F_YX] + orig_tensor[K_S_XZ]*rotation[K_F_ZX];
  t_const[1] = orig_tensor[K_S_YX]*rotation[K_F_XX] + orig_tensor[K_S_YY]*rotation[K_F_YX] + orig_tensor[K_S_YZ]*rotation[K_F_ZX];
  t_const[2] = orig_tensor[K_S_ZX]*rotation[K_F_XX] + orig_tensor[K_S_ZY]*rotation[K_F_YX] + orig_tensor[K_S_ZZ]*rotation[K_F_ZX];
  t_const[3] = orig_tensor[K_S_XX]*rotation[K_F_XY] + orig_tensor[K_S_XY]*rotation[K_F_YY] + orig_tensor[K_S_XZ]*rotation[K_F_ZY];
  t_const[4] = orig_tensor[K_S_YX]*rotation[K_F_XY] + orig_tensor[K_S_YY]*rotation[K_F_YY] + orig_tensor[K_S_YZ]*rotation[K_F_ZY];
  t_const[5] = orig_tensor[K_S_ZX]*rotation[K_F_XY] + orig_tensor[K_S_ZY]*rotation[K_F_YY] + orig_tensor[K_S_ZZ]*rotation[K_F_ZY];
  t_const[6] = orig_tensor[K_S_XX]*rotation[K_F_XZ] + orig_tensor[K_S_XY]*rotation[K_F_YZ] + orig_tensor[K_S_XZ]*rotation[K_F_ZZ];
  t_const[7] = orig_tensor[K_S_YX]*rotation[K_F_XZ] + orig_tensor[K_S_YY]*rotation[K_F_YZ] + orig_tensor[K_S_YZ]*rotation[K_F_ZZ];
  t_const[8] = orig_tensor[K_S_ZX]*rotation[K_F_XZ] + orig_tensor[K_S_ZY]*rotation[K_F_YZ] + orig_tensor[K_S_ZZ]*rotation[K_F_ZZ];

  rot_tensor[K_S_XX] = rotation[K_F_XX]*t_const[0] + rotation[K_F_YX]*t_const[1] + rotation[K_F_ZX]*t_const[2];
  rot_tensor[K_S_YY] = rotation[K_F_XY]*t_const[3] + rotation[K_F_YY]*t_const[4] + rotation[K_F_ZY]*t_const[5];
  rot_tensor[K_S_ZZ] = rotation[K_F_XZ]*t_const[6] + rotation[K_F_YZ]*t_const[7] + rotation[K_F_ZZ]*t_const[8];

  rot_tensor[K_S_XY] = rotation[K_F_XX]*t_const[3] + rotation[K_F_YX]*t_const[4] + rotation[K_F_ZX]*t_const[5];
  rot_tensor[K_S_YZ] = rotation[K_F_XY]*t_const[6] + rotation[K_F_YY]*t_const[7] + rotation[K_F_ZY]*t_const[8];
  rot_tensor[K_S_ZX] = rotation[K_F_XZ]*t_const[0] + rotation[K_F_YZ]*t_const[1] + rotation[K_F_ZZ]*t_const[2];
}

//***************************************************************************************************************
//
//  Use a rotation matrix to rotate a given tensor backward
//
//    Input:
//      rotation(9)           Rotation tensor
//      orig_tensor(6)        Input symmetric tensor
//
//   Output:
//      rot_tensor(6)         Output symmetric tensor
//
//***************************************************************************************************************
TEMPLATE_SCALAR
INLINE void rotate_tensor_backward(const Scalar *const rotation,
                                   const Scalar *const orig_tensor,
                                   Scalar *const rot_tensor) {
  Scalar t_const[9];

  t_const[0] = orig_tensor[K_S_XX]*rotation[K_F_XX] + orig_tensor[K_S_XY]*rotation[K_F_XY] + orig_tensor[K_S_XZ]*rotation[K_F_XZ];
  t_const[1] = orig_tensor[K_S_YX]*rotation[K_F_XX] + orig_tensor[K_S_YY]*rotation[K_F_XY] + orig_tensor[K_S_YZ]*rotation[K_F_XZ];
  t_const[2] = orig_tensor[K_S_ZX]*rotation[K_F_XX] + orig_tensor[K_S_ZY]*rotation[K_F_XY] + orig_tensor[K_S_ZZ]*rotation[K_F_XZ];
  t_const[3] = orig_tensor[K_S_XX]*rotation[K_F_YX] + orig_tensor[K_S_XY]*rotation[K_F_YY] + orig_tensor[K_S_XZ]*rotation[K_F_YZ];
  t_const[4] = orig_tensor[K_S_YX]*rotation[K_F_YX] + orig_tensor[K_S_YY]*rotation[K_F_YY] + orig_tensor[K_S_YZ]*rotation[K_F_YZ];
  t_const[5] = orig_tensor[K_S_ZX]*rotation[K_F_YX] + orig_tensor[K_S_ZY]*rotation[K_F_YY] + orig_tensor[K_S_ZZ]*rotation[K_F_YZ];
  t_const[6] = orig_tensor[K_S_XX]*rotation[K_F_ZX] + orig_tensor[K_S_XY]*rotation[K_F_ZY] + orig_tensor[K_S_XZ]*rotation[K_F_ZZ];
  t_const[7] = orig_tensor[K_S_YX]*rotation[K_F_ZX] + orig_tensor[K_S_YY]*rotation[K_F_ZY] + orig_tensor[K_S_YZ]*rotation[K_F_ZZ];
  t_const[8] = orig_tensor[K_S_ZX]*rotation[K_F_ZX] + orig_tensor[K_S_ZY]*rotation[K_F_ZY] + orig_tensor[K_S_ZZ]*rotation[K_F_ZZ];

  rot_tensor[K_S_XX] = rotation[K_F_XX]*t_const[0] + rotation[K_F_XY]*t_const[1] + rotation[K_F_XZ]*t_const[2];
  rot_tensor[K_S_YY] = rotation[K_F_YX]*t_const[3] + rotation[K_F_YY]*t_const[4] + rotation[K_F_YZ]*t_const[5];
  rot_tensor[K_S_ZZ] = rotation[K_F_ZX]*t_const[6] + rotation[K_F_ZY]*t_const[7] + rotation[K_F_ZZ]*t_const[8];

  rot_tensor[K_S_XY] = rotation[K_F_XX]*t_const[3] + rotation[K_F_XY]*t_const[4] + rotation[K_F_XZ]*t_const[5];
  rot_tensor[K_S_YZ] = rotation[K_F_YX]*t_const[6] + rotation[K_F_YY]*t_const[7] + rotation[K_F_YZ]*t_const[8];
  rot_tensor[K_S_ZX] = rotation[K_F_ZX]*t_const[0] + rotation[K_F_ZY]*t_const[1] + rotation[K_F_ZZ]*t_const[2];
}

//***************************************************************************************************************
//
//  Use a rotation matrix to rotate a given tensor forward
//
//    Input:
//      rotation(4)           Rotation tensor
//      orig_tensor(3)        Input symmetric tensor
//
//   Output:
//      rot_tensor(3)         Output symmetric tensor
//
//***************************************************************************************************************
TEMPLATE_SCALAR
INLINE void rotate_tensor_21_forward(const Scalar *const rotation, const Scalar *const orig_tensor, Scalar *const rot_tensor){
  Scalar t_const[4];

  t_const[0] = orig_tensor[K_S_21_XX]*rotation[K_F_22_XX] + orig_tensor[K_S_21_XY]*rotation[K_F_22_YX];
  t_const[1] = orig_tensor[K_S_21_YX]*rotation[K_F_22_XX] + orig_tensor[K_S_21_YY]*rotation[K_F_22_YX];
  t_const[2] = orig_tensor[K_S_21_XX]*rotation[K_F_22_XY] + orig_tensor[K_S_21_XY]*rotation[K_F_22_YY];
  t_const[3] = orig_tensor[K_S_21_YX]*rotation[K_F_22_XY] + orig_tensor[K_S_21_YY]*rotation[K_F_22_YY];

  rot_tensor[K_S_21_XX] = rotation[K_F_22_XX]*t_const[0] + rotation[K_F_22_YX]*t_const[1];
  rot_tensor[K_S_21_YY] = rotation[K_F_22_XY]*t_const[2] + rotation[K_F_22_YY]*t_const[3];
  rot_tensor[K_S_21_XY] = rotation[K_F_22_XX]*t_const[2] + rotation[K_F_22_YX]*t_const[3];
}

//***************************************************************************************************************
//
//  Use a rotation matrix to rotate a given tensor backward
//
//    Input:
//      rotation(4)           Rotation tensor
//      orig_tensor(3)        Input symmetric tensor
//
//   Output:
//      rot_tensor(3)         Output symmetric tensor
//
//***************************************************************************************************************
TEMPLATE_SCALAR
INLINE void rotate_tensor_21_backward(const Scalar *const rotation, const Scalar *const orig_tensor, Scalar *const rot_tensor)
{
  Scalar t_const[4];

  t_const[0] = orig_tensor[K_S_21_XX]*rotation[K_F_22_XX] + orig_tensor[K_S_21_XY]*rotation[K_F_22_XY];
  t_const[1] = orig_tensor[K_S_21_YX]*rotation[K_F_22_XX] + orig_tensor[K_S_21_YY]*rotation[K_F_22_XY];
  t_const[2] = orig_tensor[K_S_21_XX]*rotation[K_F_22_YX] + orig_tensor[K_S_21_XY]*rotation[K_F_22_YY];
  t_const[3] = orig_tensor[K_S_21_YX]*rotation[K_F_22_YX] + orig_tensor[K_S_21_YY]*rotation[K_F_22_YY];

  rot_tensor[K_S_21_XX] = rotation[K_F_22_XX]*t_const[0] + rotation[K_F_22_XY]*t_const[1];
  rot_tensor[K_S_21_YY] = rotation[K_F_22_YX]*t_const[2] + rotation[K_F_22_YY]*t_const[3];
  rot_tensor[K_S_21_XY] = rotation[K_F_22_XX]*t_const[2] + rotation[K_F_22_XY]*t_const[3];
}

/**
 * Compute the dot product of two length 8 vectors
 *
 *    Input
 *      x1(8)           first vector
 *      x2(8)           second vector
 */
TEMPLATE_SCALAR
inline Scalar ddot8(const Scalar *const x1, const Scalar *const x2) {
  return (x1[0] * x2[0] +
          x1[1] * x2[1] +
          x1[2] * x2[2] +
          x1[3] * x2[3] +
          x1[4] * x2[4] +
          x1[5] * x2[5] +
          x1[6] * x2[6] +
          x1[7] * x2[7]);
}

/**
 * Compute the dot product of two length 3 vectors
 *
 *    Input
 *      x1(3)           first vector
 *      x2(3)           second vector
 */
TEMPLATE_SCALAR
inline Scalar ddot3(const Scalar *const x1, const Scalar *const x2) {
  return (x1[0] * x2[0] +
          x1[1] * x2[1] +
          x1[2] * x2[2]);
}

/**
 * Compute the dot product of two length 8 vectors using vector striding
 *
 *    Input
 *      x1(8)           first vector
 *      x2(8)           second vector
 */
TEMPLATE_SCALAR
inline Scalar ddot8_stride1_stride3(const Scalar *const x1, const Scalar *const x2) {
  return (x1[0] * x2[ 0] +
          x1[1] * x2[ 3] +
          x1[2] * x2[ 6] +
          x1[3] * x2[ 9] +
          x1[4] * x2[12] +
          x1[5] * x2[15] +
          x1[6] * x2[18] +
          x1[7] * x2[21]);
}

/**
 *  Rotate a [3,8] matrix into a [8,3] matrix
 *
 *    Input:
 *      cur_coords        Old matrix
 *
 *    Output:
 *      cordel_ptr        New matrix
 */
TEMPLATE_SCALAR
INLINE void transform_38_matrix_to_83_matrix(const Scalar *cordel_ptr, Scalar *cur_coords) {
  cur_coords[0] = cordel_ptr[0];
  cur_coords[1] = cordel_ptr[3];
  cur_coords[2] = cordel_ptr[6];
  cur_coords[3] = cordel_ptr[9];
  cur_coords[4] = cordel_ptr[12];
  cur_coords[5] = cordel_ptr[15];
  cur_coords[6] = cordel_ptr[18];
  cur_coords[7] = cordel_ptr[21];

  cur_coords[8] = cordel_ptr[1];
  cur_coords[9] = cordel_ptr[4];
  cur_coords[10] = cordel_ptr[7];
  cur_coords[11] = cordel_ptr[10];
  cur_coords[12] = cordel_ptr[13];
  cur_coords[13] = cordel_ptr[16];
  cur_coords[14] = cordel_ptr[19];
  cur_coords[15] = cordel_ptr[22];

  cur_coords[16] = cordel_ptr[2];
  cur_coords[17] = cordel_ptr[5];
  cur_coords[18] = cordel_ptr[8];
  cur_coords[19] = cordel_ptr[11];
  cur_coords[20] = cordel_ptr[14];
  cur_coords[21] = cordel_ptr[17];
  cur_coords[22] = cordel_ptr[20];
  cur_coords[23] = cordel_ptr[23];
}

/**
 *  Add bulk viscosity terms to a set of symmetric 3x3 tensors.  Bulk viscosity terms
 *  will only be added to the terms on the diagonal of the tensor.
 *
 *    Input:
 *      Int num_elements   :  Number of input tensors
 *      Real *bulkq        :  Bulk viscosity term to add
 *      Real *stress_in    :  Input set of tensors
 *
 *    Output:
 *      Real *stress_out   :  Modified output tensors
 */
TEMPLATE_SCALAR
INLINE void add_bulk_visc(const Int num_elements, const Scalar * bulkq, const Scalar * stress_in, Scalar * stress_out) {
  for(Int ielem = 0; ielem < num_elements; ++ielem) {
    stress_out[K_S_XX] = stress_in[K_S_XX] + bulkq[ielem];
    stress_out[K_S_YY] = stress_in[K_S_YY] + bulkq[ielem];
    stress_out[K_S_ZZ] = stress_in[K_S_ZZ] + bulkq[ielem];
    stress_out[K_S_XY] = stress_in[K_S_XY];
    stress_out[K_S_XZ] = stress_in[K_S_XZ];
    stress_out[K_S_YZ] = stress_in[K_S_YZ];
    stress_out += 6;
    stress_in  += 6;
  }
}

/**
 *  Compute the aspect ratio for an element
 *
 *    Input:
 *      gradop12x(24)   12 times the gradient operator
 *      volume12x       12 times the element volume
 */   
TEMPLATE_SCALAR
INLINE Scalar comp_aspect(const Scalar *const gradop12x, const Scalar &volume12x) {
  return 6.0 * volume12x /
    (gradop12x[0] * gradop12x[0]+
     gradop12x[1] * gradop12x[1]+
     gradop12x[2] * gradop12x[2]+
     gradop12x[3] * gradop12x[3]+
     gradop12x[4] * gradop12x[4]+
     gradop12x[5] * gradop12x[5]+
     gradop12x[6] * gradop12x[6]+
     gradop12x[7] * gradop12x[7]+
     gradop12x[8] * gradop12x[8]+
     gradop12x[9] * gradop12x[9]+
     gradop12x[10] * gradop12x[10]+
     gradop12x[11] * gradop12x[11]+
     gradop12x[12] * gradop12x[12]+
     gradop12x[13] * gradop12x[13]+
     gradop12x[14] * gradop12x[14]+
     gradop12x[15] * gradop12x[15]+
     gradop12x[16] * gradop12x[16]+
     gradop12x[17] * gradop12x[17]+
     gradop12x[18] * gradop12x[18]+
     gradop12x[19] * gradop12x[19]+
     gradop12x[20] * gradop12x[20]+
     gradop12x[21] * gradop12x[21]+
     gradop12x[22] * gradop12x[22]+
     gradop12x[23] * gradop12x[23]);
}


/**
 * **********************************************************************************************
 *  
 *  Purpose:  Perform all element calculations after the element stress calculation
 *            
 *    Input:
 *      Int   nelem                Number of elements in workset
 *      Int   dt                   Last step time step
 *      Real  *cordel              Element nodal coordinates
 *      Real  *vel                 Element nodal velocities
 *      Real  *rotation            elemetn rotation tensor
 *      Real  *stress_new          The new element stress
 *      Real  *rotated_stretching  The element stretching rate tensor
 *      Real  *spin_rate           rigid body elemetn spining rate
 *      Real  linBulkVisc          Linear artificial bulk viscosity parameter      
 *      Real  quadBulkVisc         Quadradic bulk viscosity parameter
 *      Real  *elem_mass           The element mass
 *      Real  *elem_dilmod         Element dilatation modulus as computed by the effective modulus routine
 *      Real  *elem_shrmod         Element shear modulus as coputed by the effective modulus routine
 *      Real  hg_stiffness         Element hourglass stiffness scaling constant
 *      Real  hg_viscosity         Element hourglass viscosity scaling constant
 *      
 *   Output:
 *      Real  *rotated_stress      Element rotated stress tensor
 *      Real  min_elem_time_step   New mimumum element based time step
 *      Real  *volume              Element end of step volume
 *      Real  *elem_time_step      Element time steps on per element basis
 *      Real  *hg_resist_old       Old hourglass resitance rate
 *      Real  *hg_resist_new       New hourglass resistance rate 
 *      Real  *force_new           Computed internal force 
 *      Real  *hg_energy           Per element hourglass force energy
 *      Real  *int_energy          Per element internal energy
 *
 * **********************************************************************************************
 */

//***************************************************************************************************************
//
//  Compute a set of 8 components for a gradient operator.
//
//  48*3 = 144 multiplications.
//
//    Input
//      grad_ptr(24)   Array in which to store the current gradient operator components
//      x(8)           X nodal coordinates
//      y(8)           Y nodal coorindates
//      z(8)           Z nodal coordinates
//
//***************************************************************************************************************
TEMPLATE_SCALAR
INLINE 
void comp_grad12x( Scalar *const grad_ptr, const Scalar * const x, const Scalar *const y, const Scalar *const z) {
  Scalar R42=(z[3] - z[1]);
  Scalar R52=(z[4] - z[1]);
  Scalar R54=(z[4] - z[3]);

  Scalar R63=(z[5] - z[2]);
  Scalar R83=(z[7] - z[2]);
  Scalar R86=(z[7] - z[5]);

  Scalar R31=(z[2] - z[0]);
  Scalar R61=(z[5] - z[0]);
  Scalar R74=(z[6] - z[3]);

  Scalar R72=(z[6] - z[1]);
  Scalar R75=(z[6] - z[4]);
  Scalar R81=(z[7] - z[0]);

  Scalar t1 =(R63 + R54);
  Scalar t2 =(R61 + R74);
  Scalar t3 =(R72 + R81);
  Scalar t4 =(R86 + R42);
  Scalar t5 =(R83 + R52);
  Scalar t6 =(R75 + R31);

  grad_ptr[0] = (y[1] *  t1) - (y[2] * R42) - (y[3] *  t5) + (y[4] *  t4) + (y[5] * R52) - (y[7] * R54);
  grad_ptr[1] = (y[2] *  t2) + (y[3] * R31) - (y[0] *  t1) - (y[5] *  t6) + (y[6] * R63) - (y[4] * R61);
  grad_ptr[2] = (y[3] *  t3) + (y[0] * R42) - (y[1] *  t2) - (y[6] *  t4) + (y[7] * R74) - (y[5] * R72);
  grad_ptr[3] = (y[0] *  t5) - (y[1] * R31) - (y[2] *  t3) + (y[7] *  t6) + (y[4] * R81) - (y[6] * R83);
  grad_ptr[4] = (y[5] *  t3) + (y[6] * R86) - (y[7] *  t2) - (y[0] *  t4) - (y[3] * R81) + (y[1] * R61);
  grad_ptr[5] = (y[6] *  t5) - (y[4] *  t3) - (y[7] * R75) + (y[1] *  t6) - (y[0] * R52) + (y[2] * R72);
  grad_ptr[6] = (y[7] *  t1) - (y[5] *  t5) - (y[4] * R86) + (y[2] *  t4) - (y[1] * R63) + (y[3] * R83);
  grad_ptr[7] = (y[4] *  t2) - (y[6] *  t1) + (y[5] * R75) - (y[3] *  t6) - (y[2] * R74) + (y[0] * R54);

  R42=(x[3] - x[1]);
  R52=(x[4] - x[1]);
  R54=(x[4] - x[3]);

  R63=(x[5] - x[2]);
  R83=(x[7] - x[2]);
  R86=(x[7] - x[5]);

  R31=(x[2] - x[0]);
  R61=(x[5] - x[0]);
  R74=(x[6] - x[3]);

  R72=(x[6] - x[1]);
  R75=(x[6] - x[4]);
  R81=(x[7] - x[0]);

  t1 =(R63 + R54);
  t2 =(R61 + R74);
  t3 =(R72 + R81);
  t4 =(R86 + R42);
  t5 =(R83 + R52);
  t6 =(R75 + R31);

  grad_ptr[8 ] = (z[1] *  t1) - (z[2] * R42) - (z[3] *  t5) + (z[4] *  t4) + (z[5] * R52) - (z[7] * R54);
  grad_ptr[9 ] = (z[2] *  t2) + (z[3] * R31) - (z[0] *  t1) - (z[5] *  t6) + (z[6] * R63) - (z[4] * R61);
  grad_ptr[10] = (z[3] *  t3) + (z[0] * R42) - (z[1] *  t2) - (z[6] *  t4) + (z[7] * R74) - (z[5] * R72);
  grad_ptr[11] = (z[0] *  t5) - (z[1] * R31) - (z[2] *  t3) + (z[7] *  t6) + (z[4] * R81) - (z[6] * R83);
  grad_ptr[12] = (z[5] *  t3) + (z[6] * R86) - (z[7] *  t2) - (z[0] *  t4) - (z[3] * R81) + (z[1] * R61);
  grad_ptr[13] = (z[6] *  t5) - (z[4] *  t3) - (z[7] * R75) + (z[1] *  t6) - (z[0] * R52) + (z[2] * R72);
  grad_ptr[14] = (z[7] *  t1) - (z[5] *  t5) - (z[4] * R86) + (z[2] *  t4) - (z[1] * R63) + (z[3] * R83);
  grad_ptr[15] = (z[4] *  t2) - (z[6] *  t1) + (z[5] * R75) - (z[3] *  t6) - (z[2] * R74) + (z[0] * R54);

  R42=(y[3] - y[1]);
  R52=(y[4] - y[1]);
  R54=(y[4] - y[3]);

  R63=(y[5] - y[2]);
  R83=(y[7] - y[2]);
  R86=(y[7] - y[5]);

  R31=(y[2] - y[0]);
  R61=(y[5] - y[0]);
  R74=(y[6] - y[3]);

  R72=(y[6] - y[1]);
  R75=(y[6] - y[4]);
  R81=(y[7] - y[0]);

  t1 =(R63 + R54);
  t2 =(R61 + R74);
  t3 =(R72 + R81);
  t4 =(R86 + R42);
  t5 =(R83 + R52);
  t6 =(R75 + R31);

  grad_ptr[16] = (x[1] *  t1) - (x[2] * R42) - (x[3] *  t5) + (x[4] *  t4) + (x[5] * R52) - (x[7] * R54);
  grad_ptr[17] = (x[2] *  t2) + (x[3] * R31) - (x[0] *  t1) - (x[5] *  t6) + (x[6] * R63) - (x[4] * R61);
  grad_ptr[18] = (x[3] *  t3) + (x[0] * R42) - (x[1] *  t2) - (x[6] *  t4) + (x[7] * R74) - (x[5] * R72);
  grad_ptr[19] = (x[0] *  t5) - (x[1] * R31) - (x[2] *  t3) + (x[7] *  t6) + (x[4] * R81) - (x[6] * R83);
  grad_ptr[20] = (x[5] *  t3) + (x[6] * R86) - (x[7] *  t2) - (x[0] *  t4) - (x[3] * R81) + (x[1] * R61);
  grad_ptr[21] = (x[6] *  t5) - (x[4] *  t3) - (x[7] * R75) + (x[1] *  t6) - (x[0] * R52) + (x[2] * R72);
  grad_ptr[22] = (x[7] *  t1) - (x[5] *  t5) - (x[4] * R86) + (x[2] *  t4) - (x[1] * R63) + (x[3] * R83);
  grad_ptr[23] = (x[4] *  t2) - (x[6] *  t1) + (x[5] * R75) - (x[3] *  t6) - (x[2] * R74) + (x[0] * R54);
}

//***************************************************************************************************************
//
//  Compute the hourglass operator
//
//    Input:
//      x_ptr, y_ptr, z_ptr        The x, y, and z nodal coordinates
//
//      gradop12x_ptr_x            12 times the gradient operator components
//      gradop12x_ptr_y            12 times the gradient operator components
//      gradop12x_ptr_z            12 times the gradient operator components
//
//      volinv12th                 The inverse of the element voluem, times one 12th
//
//***************************************************************************************************************
TEMPLATE_SCALAR
INLINE
void comp_hgop(const Scalar *const x, const Scalar *const y, const Scalar * const z,
               const Scalar *const gradop12x_ptr_x,
               const Scalar *const gradop12x_ptr_y,
               const Scalar *const gradop12x_ptr_z,
               const Scalar volinv12th,
               Scalar *const hgop) {

  //
  // KHP: Perhaps, this static array should be passed in from the calling
  // routine?
  //
  static const Scalar hgop_arr[32] = {
    1.0, 1.0,-1.0,-1.0,-1.0,-1.0, 1.0, 1.0,
    1.0,-1.0,-1.0, 1.0,-1.0, 1.0, 1.0,-1.0,
    1.0,-1.0, 1.0,-1.0, 1.0,-1.0, 1.0,-1.0,
    -1.0, 1.0,-1.0, 1.0, 1.0,-1.0, 1.0,-1.0};

  // KHP: Alternatively, we could have
  // hx0,hx1,hx2,hx3,...,hz0,hz1,hz2,hz3
  Scalar hgconst12th[12];

  Scalar t0 = x[0] - x[1];
  Scalar t1 = x[2] - x[3];
  Scalar t2 = x[4] - x[5];
  Scalar t3 = x[6] - x[7];

  hgconst12th[0] = ( (x[0]+x[1]) - (x[2]+x[3]) - (x[4]+x[5]) + (x[6]+x[7]) ) * volinv12th;
  hgconst12th[1] = (  t0 - t1 - t2 + t3 ) * volinv12th;
  hgconst12th[2] = (  t0 + t1 + t2 + t3 ) * volinv12th;
  hgconst12th[3] = ( -t0 - t1 + t2 + t3 ) * volinv12th;

  t0 = (y[0] - y[1]);
  t1 = (y[2] - y[3]);
  t2 = (y[4] - y[5]);
  t3 = (y[6] - y[7]);

  hgconst12th[4] = ( (y[0]+y[1]) - (y[2]+y[3]) - (y[4]+y[5]) + (y[6]+y[7]) ) * volinv12th;
  hgconst12th[5] = (  t0 - t1 - t2 + t3 ) * volinv12th;
  hgconst12th[6] = (  t0 + t1 + t2 + t3 ) * volinv12th;
  hgconst12th[7] = ( -t0 - t1 + t2 + t3 ) * volinv12th;

  t0 = (z[0] - z[1]);
  t1 = (z[2] - z[3]);
  t2 = (z[4] - z[5]);
  t3 = (z[6] - z[7]);

  hgconst12th[8]  = ( (z[0]+z[1]) - (z[2]+z[3]) - (z[4]+z[5]) + (z[6]+z[7]) ) * volinv12th;
  hgconst12th[9]  = (  t0 - t1 - t2 + t3 ) * volinv12th;
  hgconst12th[10] = (  t0 + t1 + t2 + t3 ) * volinv12th;
  hgconst12th[11] = ( -t0 - t1 + t2 + t3 ) * volinv12th;

  // 8 times 12 = 96 + 12 = 108 multiplications.
  for(int i = 0; i < 8; ++i) {
    hgop[i   ] = hgop_arr[i   ] - (hgconst12th[0] * gradop12x_ptr_x[i]  + hgconst12th[4] * gradop12x_ptr_y[i] + hgconst12th[8 ] * gradop12x_ptr_z[i]);
    hgop[i+8 ] = hgop_arr[i+8 ] - (hgconst12th[1] * gradop12x_ptr_x[i]  + hgconst12th[5] * gradop12x_ptr_y[i] + hgconst12th[9 ] * gradop12x_ptr_z[i]);
    hgop[i+16] = hgop_arr[i+16] - (hgconst12th[2] * gradop12x_ptr_x[i]  + hgconst12th[6] * gradop12x_ptr_y[i] + hgconst12th[10] * gradop12x_ptr_z[i]);
    hgop[i+24] = hgop_arr[i+24] - (hgconst12th[3] * gradop12x_ptr_x[i]  + hgconst12th[7] * gradop12x_ptr_y[i] + hgconst12th[11] * gradop12x_ptr_z[i]);
  }
}

//*************************************************************************************************
//
//  Compute the new internal forces, hg resistances, and element energies
//
//    Input:
//      dt                    last time step
//      spin_rate_ptr         pointer into the spin rate array
//      hgop                  element hourglass operator
//      vel_ptr               element velocity
//      fac1                  hourglass stiffness resitance factor
//      fac2                  hourglass viscosity resitance factor
//      hg_resist_old_ptr     Old hourglass resitances
//
//   Output:
//      hg_resist_new_ptr     New hourglass resitances
//      hg_resist_total_ptr   Total hourglass resistance
//
//*************************************************************************************************

TEMPLATE_SCALAR
INLINE
void  comp_force(const Scalar &dt,
                 const Scalar *const spin_rate_ptr,
                 const Scalar *const hgop_for_resist_calc,
                 const Scalar *const hgop_ptr0,
                 const Scalar *const hgop_ptr1,
                 const Scalar *const hgop_ptr2,
                 const Scalar *const hgop_ptr3,
                 const Scalar *const vel,
                 const Scalar &fac1,
                 const Scalar &fac2,
                 const Scalar *const hg_resist_old,
                 Scalar *const hg_resist_new,
                 const Scalar *const total_stress12th,
                 const Scalar *const gradop12x_ptr_x,
                 const Scalar *const gradop12x_ptr_y,
                 const Scalar *const gradop12x_ptr_z,
                 Scalar *const hg_energy,
                 Scalar *const int_energy,
                 Scalar *const force_new,
                 bool scaleHGRotation) {

  //
  //  NKC, does Presto Scalar need these spin rate terms?  Pronto appears to have dumped them....
  //
  const Scalar dwxy = dt * spin_rate_ptr[0];
  const Scalar dwyz = dt * spin_rate_ptr[1];
  const Scalar dwzx = dt * spin_rate_ptr[2];
  //
  //  Compute new hourglass resitance by the old rotated hourglass resitance plus a hourglass rate term
  //
  Scalar hg_resist_total[12];
  Scalar * hg_resist_total_ptr = hg_resist_total;
  const Scalar * hgop_ptr = hgop_for_resist_calc;
  const Scalar * vel_ptr = vel;
  const Scalar * hg_resist_old_ptr = hg_resist_old;
  Scalar * hg_resist_new_ptr = hg_resist_new;
  Scalar * hg_energy_ptr = hg_energy;
  Scalar * int_energy_ptr = int_energy;
  Scalar * force_new_ptr = force_new;

  if (!scaleHGRotation){
    for(int i = 0; i < 4; ++i) {

      const Scalar hg_resist_old_0 = hg_resist_old_ptr[0];
      const Scalar hg_resist_old_1 = hg_resist_old_ptr[1];
      const Scalar hg_resist_old_2 = hg_resist_old_ptr[2];

      const Scalar hg_rate_0 = ddot8_stride1_stride3(hgop_ptr, vel_ptr);
      hg_resist_new_ptr[0] = hg_resist_old_0 + dwxy*hg_resist_old_1 - dwzx*hg_resist_old_2 + fac1* hg_rate_0 ;
      hg_resist_total_ptr[0]   = hg_resist_new_ptr[0] + fac2* hg_rate_0 ;

      const Scalar hg_rate_1 = ddot8_stride1_stride3(hgop_ptr, vel_ptr+1);
      hg_resist_new_ptr[1] = hg_resist_old_1 - dwxy*hg_resist_old_0 + dwyz*hg_resist_old_2 + fac1* hg_rate_1 ;
      hg_resist_total_ptr[1]   = hg_resist_new_ptr[1] + fac2* hg_rate_1 ;

      const Scalar hg_rate_2 = ddot8_stride1_stride3(hgop_ptr, vel_ptr+2);
      hg_resist_new_ptr[2] = hg_resist_old_2 + dwzx*hg_resist_old_0 - dwyz*hg_resist_old_1 + fac1* hg_rate_2 ;
      hg_resist_total_ptr[2]   = hg_resist_new_ptr[2] + fac2* hg_rate_2 ;

      hgop_ptr += 8;
      hg_resist_new_ptr+= 3;
      hg_resist_old_ptr+=3;
      hg_resist_total_ptr += 3;
    }
  } else {
    for(int i = 0; i < 4; ++i) {
      const Scalar hg_rate_0 = ddot8_stride1_stride3(hgop_ptr, vel_ptr);
      const Scalar hg_rate_1 = ddot8_stride1_stride3(hgop_ptr, vel_ptr+1);
      const Scalar hg_rate_2 = ddot8_stride1_stride3(hgop_ptr, vel_ptr+2);

      const Scalar hg_resist_old_0 = hg_resist_old_ptr[0];
      const Scalar hg_resist_old_1 = hg_resist_old_ptr[1];
      const Scalar hg_resist_old_2 = hg_resist_old_ptr[2];

      const Scalar rot_hg_resist_old_0 = hg_resist_old_0 + dwxy*hg_resist_old_1 - dwzx*hg_resist_old_2;
      const Scalar rot_hg_resist_old_1 = hg_resist_old_1 - dwxy*hg_resist_old_0 + dwyz*hg_resist_old_2;
      const Scalar rot_hg_resist_old_2 = hg_resist_old_2 + dwzx*hg_resist_old_0 - dwyz*hg_resist_old_1;

      Scalar fnorm = rot_hg_resist_old_0 *rot_hg_resist_old_0  + rot_hg_resist_old_1 *rot_hg_resist_old_1  + rot_hg_resist_old_2 *rot_hg_resist_old_2 ;
      if (fnorm > 1.e-30){
        fnorm = sqrt ( (hg_resist_old_0*hg_resist_old_0 +
                        hg_resist_old_1*hg_resist_old_1 +
                        hg_resist_old_2*hg_resist_old_2) / fnorm );
        hg_resist_new_ptr[0] = fnorm*rot_hg_resist_old_0 + fac1*hg_rate_0 ;
        hg_resist_new_ptr[1] = fnorm*rot_hg_resist_old_1 + fac1*hg_rate_1 ;
        hg_resist_new_ptr[2] = fnorm*rot_hg_resist_old_2 + fac1*hg_rate_2 ;
      } else {
        hg_resist_new_ptr[0] = rot_hg_resist_old_0 + fac1*hg_rate_0 ;
        hg_resist_new_ptr[1] = rot_hg_resist_old_1 + fac1*hg_rate_1 ;
        hg_resist_new_ptr[2] = rot_hg_resist_old_2 + fac1*hg_rate_2 ;
      }
      hg_resist_total_ptr[0] = hg_resist_new_ptr[0] + fac2*hg_rate_0 ;
      hg_resist_total_ptr[1] = hg_resist_new_ptr[1] + fac2*hg_rate_1 ;
      hg_resist_total_ptr[2] = hg_resist_new_ptr[2] + fac2*hg_rate_2 ;
      hgop_ptr += 8;
      hg_resist_new_ptr+= 3;
      hg_resist_old_ptr+=3;
      hg_resist_total_ptr += 3;
    }
  }

  Scalar hg_force_0[8];
  Scalar hg_force_1[8];
  Scalar hg_force_2[8];
  
  *hg_energy_ptr = 0.0;
  *int_energy_ptr = 0.0;
  for(int i = 0; i < 8; ++i) {
    hg_force_0[i] =
      (hg_resist_total[HG_X1] * hgop_ptr0[i] +
       hg_resist_total[HG_X2] * hgop_ptr1[i] +
       hg_resist_total[HG_X3] * hgop_ptr2[i] +
       hg_resist_total[HG_X4] * hgop_ptr3[i]);

    hg_force_1[i] =
      (hg_resist_total[HG_Y1] * hgop_ptr0[i] +
       hg_resist_total[HG_Y2] * hgop_ptr1[i] +
       hg_resist_total[HG_Y3] * hgop_ptr2[i] +
       hg_resist_total[HG_Y4] * hgop_ptr3[i]);

    hg_force_2[i] =
      (hg_resist_total[HG_Z1] * hgop_ptr0[i] +
       hg_resist_total[HG_Z2] * hgop_ptr1[i] +
       hg_resist_total[HG_Z3] * hgop_ptr2[i] +
       hg_resist_total[HG_Z4] * hgop_ptr3[i]);
  }
  
  for(int i = 0; i < 8; ++i) {
    force_new_ptr[0] =
      total_stress12th[K_S_XX] * gradop12x_ptr_x[i] +
      total_stress12th[K_S_XY] * gradop12x_ptr_y[i] +
      total_stress12th[K_S_XZ] * gradop12x_ptr_z[i] + hg_force_0[i] ;

    force_new_ptr[1] =
      total_stress12th[K_S_YX] * gradop12x_ptr_x[i] +
      total_stress12th[K_S_YY] * gradop12x_ptr_y[i] +
      total_stress12th[K_S_YZ] * gradop12x_ptr_z[i] + hg_force_1[i] ;

    force_new_ptr[2] =
      total_stress12th[K_S_ZX] * gradop12x_ptr_x[i] +
      total_stress12th[K_S_ZY] * gradop12x_ptr_y[i] +
      total_stress12th[K_S_ZZ] * gradop12x_ptr_z[i] + hg_force_2[i] ;

    *hg_energy_ptr  += hg_force_0[i]      *vel_ptr[0] + hg_force_1[i]      *vel_ptr[1] + hg_force_2[i]      *vel_ptr[2];
    *int_energy_ptr += force_new_ptr[0]*vel_ptr[0] + force_new_ptr[1]*vel_ptr[1] + force_new_ptr[2]*vel_ptr[2];


    force_new_ptr += 3;
    vel_ptr += 3;
  }
  ++hg_energy_ptr ;
  ++int_energy_ptr ;
}

/**
 * **********************************************************************************************
 *  
 *  Purpose:  Perform all element calculations required for the material stress calculations.
 *            
 *    Input:
 *      Int   nelem             Number of elements in workset
 *      Scalar  dt_scale             Should be -0.5*delta_t, describes how to move the coordinates
 *                              from the end of step to the mid-step using the velocities.
 *      Scalar  *cordel              Array of element nodal coordinates
 *      Scalar  *vel                 Array of element nodal velocities
 *      Scalar  *rotation_old        Array of old roation tensors
 *      
 *   Output:
 *      Scalar  *mid_volume          Compute element volume from the mid step coordinates
 *      Scalar  *vorticity           Asymmetric portion of the element stretch tensor
 *      Scalar  *rotation_new        Updated material rotation tensors  
 *      Scalar  *stretch             New element stretch
 *      Real  *rotated_stretching  Rotated element stretching tensor
 *
 * **********************************************************************************************
 */
TEMPLATE_SCALAR
INLINE Int elem_ug3dh8_mi_compute_stretch(const Int nelem,
                                          const Scalar dt,
                                          const Scalar *const cordel,
                                          const Scalar *const vel,
                                          const Scalar *const rotation_old,
                                          Scalar *const mid_vol,
                                          Scalar *const vorticity,
                                          Scalar *const rotation_new,
                                          Scalar *const stretch,
                                          Scalar *const rotated_stretching,
                                          Scalar *const mid_hgop) {

  //
  //  Extract pointers to the input variables, set up any required temporary variables.  The temporaries will
  //  be reused by every element in the workset.
  //
  const Scalar * cordel_ptr  (cordel  );
  const Scalar * vel_ptr     (vel     );
  Scalar *      mid_vol_ptr (mid_vol );
  Scalar *      vorticity_ptr(vorticity);
  const Scalar * rotation_old_ptr(rotation_old);
  Scalar *      rotation_new_ptr(rotation_new);
  Scalar *      stretch_ptr(stretch);
  Scalar *      rotated_stretching_ptr(rotated_stretching);
  Scalar *      mid_hgop_ptr(mid_hgop);
  //
  //  Storage for the midstep gradient operator, 12x denotes that this will actually be the midstep gradient
  //  operator components.
  //
  Scalar mid_gradop12x[24];
  Scalar *const mid_gradop12x_ptr_x(mid_gradop12x);
  Scalar *const mid_gradop12x_ptr_y(mid_gradop12x+8);
  Scalar *const mid_gradop12x_ptr_z(mid_gradop12x+16);
  //
  //  Midstep coordinates
  //
  Scalar mid_coords[24];
  Scalar *const x_ptr(mid_coords);
  Scalar *const y_ptr(mid_coords+8);
  Scalar *const z_ptr(mid_coords+16);
  //
  //  Element stretching tensor
  //
  Scalar stretching_tensor[6];
  //
  //  Element velocity gradient
  //
  Scalar vel_grad[9];
  //
  //  Loop over all elements in the workset
  //
  Scalar dt_scale = -0.5 * dt;

  Int return_value = 0;

  for(Int ielem(0); ielem < nelem; ++ielem) {
    //
    //  Get current element coordinates and update them so that they reference the midstep
    //  configuration. Note, this routine changes the coordinate ordering from
    //
    //  (x1, y1, z1, x2, y2, z2, ........)
    //
    //    to
    //
    //  (x1, x2, x3, x4, x5, x6, x7, x8, y1, y2, y3 ......)
    //
    //  The array reordering allows better variable access later for optimization purposes
    //
    for(int i = 0; i < 8; ++i) {
      x_ptr[i] = cordel_ptr[0] + dt_scale * vel_ptr[0];
      y_ptr[i] = cordel_ptr[1] + dt_scale * vel_ptr[1];
      z_ptr[i] = cordel_ptr[2] + dt_scale * vel_ptr[2];
      //std::cout << "KHP: x= " << x_ptr[i] << ", " << y_ptr[i] << ", " << z_ptr[i] << "\n";
      //std::cout << "KHP: C= " << cordel_ptr[0] << ", " << cordel_ptr[1] << ", " << cordel_ptr[2] << "\n";
      //std::cout << "KHP: V= " << vel_ptr[0] << ", " << vel_ptr[1] << ", " << vel_ptr[2] << "\n";
      cordel_ptr += 3;
      vel_ptr    += 3;
    }

    vel_ptr -= 24;
    //
    //  Do grad op in three seperate steps using the helper routine.
    //  NKC, I'm not positive what the helper is mathematically equivalent too,
    //  however, the combination of these three steps will compute the full gradient.
    //  Note, the comp12x routine actually calculates 12.0 times the gradient, this
    //  12 will be factored out later.
    //
    comp_grad12x(mid_gradop12x, x_ptr, y_ptr, z_ptr);

    //
    //  Compute the midstep volume.  Note that 12 times the volume is being calculated,
    //  The actually output value of volume is then multipled by one twelth
    //
    const Scalar volume12x = ddot8(x_ptr, mid_gradop12x_ptr_x);
    *mid_vol_ptr = volume12x * ONE12TH;
    //
    //  if volume is <= 0.0, report an error, otherwise continue the calculation
    //
    if(volume12x > 0.0) {
      const Scalar volinv12th(1.0/volume12x);
      //
      //  Compute the velocity gradients from the dot products of the input velocities with the
      //  relevant columns of the gradient operator.  Note that the velocity gradient is calculated
      //  as Velocity * 12 * mid_gradop * volinv/12, the twelves cancel.
      //
      // 81 multiplications (9 tensor entries * (8 mults per dot + 1 mult for volinv))
      //
      vel_grad[K_F_XX] = ddot8_stride1_stride3(mid_gradop12x_ptr_x, vel_ptr  ) * volinv12th;
      vel_grad[K_F_XY] = ddot8_stride1_stride3(mid_gradop12x_ptr_y, vel_ptr  ) * volinv12th;
      vel_grad[K_F_XZ] = ddot8_stride1_stride3(mid_gradop12x_ptr_z, vel_ptr  ) * volinv12th;
      vel_grad[K_F_YZ] = ddot8_stride1_stride3(mid_gradop12x_ptr_z, vel_ptr+1) * volinv12th;
      vel_grad[K_F_YY] = ddot8_stride1_stride3(mid_gradop12x_ptr_y, vel_ptr+1) * volinv12th;
      vel_grad[K_F_YX] = ddot8_stride1_stride3(mid_gradop12x_ptr_x, vel_ptr+1) * volinv12th;
      vel_grad[K_F_ZX] = ddot8_stride1_stride3(mid_gradop12x_ptr_x, vel_ptr+2) * volinv12th;
      vel_grad[K_F_ZZ] = ddot8_stride1_stride3(mid_gradop12x_ptr_z, vel_ptr+2) * volinv12th;
      vel_grad[K_F_ZY] = ddot8_stride1_stride3(mid_gradop12x_ptr_y, vel_ptr+2) * volinv12th;
      //
      //  Compute the stretching tensors from the current velocity gradient
      //
      // Let G = Gradient (velocity gradient in this case).
      // Let S = Stretch
      // Let V = Vorticity
      // 
      // Then:
      //
      // S = 0.5 * ( G + G'); // S = 0.5 * ( G + transpose(G));
      // V = 0.5 * ( G - G'); // V = 0.5 * ( G - transpose(G));
      //
      APS::Hex::additive_decomp36(vel_grad, stretching_tensor, vorticity_ptr);

      polar_decomp33r2(dt, stretching_tensor, vorticity_ptr, rotation_old_ptr, stretch_ptr, rotation_new_ptr);

      rotate_tensor_forward(rotation_new_ptr, stretching_tensor, rotated_stretching_ptr);

      //
      //  Compute the midstep hourglass operator if needed. Here one 12th the volume is 
      //  multiplied by 12 times the gradient operators, so the true hourglass operator 
      //  is computed.
      //
      if (mid_hgop_ptr) {
        comp_hgop(x_ptr, y_ptr, z_ptr,
		  mid_gradop12x_ptr_x,
		  mid_gradop12x_ptr_y,
		  mid_gradop12x_ptr_z,
		  volinv12th, mid_hgop_ptr);
        mid_hgop_ptr += 32;
      }

    } else {
      return_value = -1;
    }
    ++mid_vol_ptr ;
    vel_ptr += 24;
    vorticity_ptr +=3;
    rotation_old_ptr += 9;
    stretch_ptr += 6;
    rotation_new_ptr += 9;
    rotated_stretching_ptr += 6;
  }
  return return_value;
}

//*********************************************************************n
//
//  STANDARD PRESTO VERSION
//
// description:
//     This routine is implemented for the 8 node brick, uniform strain
//     element. It is the divergence operator, hourglass force operator,
//     and time step calcuator.  Combination of multiple routines yields
//     better cache usage.
//
// Input:
//     nelem                    Number of elements in the workset
//
//     dt                       Previous time step
//
//     cordel(24,nelem)         End of step element coordinates, orgainzied by (x1, y1, z1, x2, y2, z2, ......)
//
//     vel(24,nelem)            Element nodal velocity, orgainzied by (x1, y1, z1, x2, y2, z2, ......)
//
//     rotated_stress(6,nelem)  Element stress in the rotated configuration
//
//     rotated_stretching(6,nelem) Element stretch rate tensor in rotated configuration
//
//     spin_rate(3,nelem)       Asymmetric portion of stretching rate, represents amount element is spinning
//                              about each axis.
//
//     linBulkVisc              Linear bulk viscosity material constant
//
//     quadBulkVisc             Quadradic bulk viscosity material constant
//
//     elem_mass(nelem)         Current element masses
//
//     elem_dilmod(nelem)       Current material dilitational modulus
//
//     elem_shrmod(nelem)       Current material shear modulus
//
//     hg_stiffness             Hourglass stiffness parameter
//
//     hg_viscosity             Hourglass viscosity parameter
//
//     hg_resist_old(12,nelem)  Old hourglass resistance
//
// Output:
//     min_elem_time_step       The computed mimimum stable element time step
//
//     volume(nelem)            The element end of step volume
//
//     elem_time_step(nelem)    Time steps on a per element basis
//
//     hg_resist_new(12,nelem)  New hourglass resistance
//
//     force_new(24,nelem)      New element total internal force
//
//     hg_energy(nelem)         New element total hourglass energy increment
//
//     int_energy(ielem)        New element internal energy increment
//
//     Returns,  0 if no errors, -1 if any negative volume elements were found
//
//***********************************************************************

TEMPLATE_SCALAR
INLINE
Int elem_ug3dh8_mi_compute_divergence_presto(const Int nelem,
                                             const Scalar dt,
                                             const Scalar *const cordel,
                                             const Scalar *const vel,
                                             const Scalar *const rotation,
                                             const Scalar *const stress_new,
                                             const Scalar *const rotated_stretching,
                                             const Scalar *const spin_rate,
                                             const Scalar linBulkVisc,
                                             const Scalar quadBulkVisc,
                                             const Scalar *const elem_mass,
                                             const Scalar *const elem_dilmod,
                                             const Scalar *const elem_shrmod,
                                             const Scalar hg_stiffness,
                                             const Scalar hg_viscosity,
                                             Scalar *const rotated_stress,
                                             Scalar &min_elem_time_step,
                                             Scalar *const volume,
                                             Scalar *const elem_time_step,
                                             Scalar *const hg_resist_old,
                                             Scalar *const hg_resist_new,
                                             Scalar *const force_new,
                                             Scalar *const hg_energy,
                                             Scalar *const int_energy,
                                             Scalar *const mid_hgop,
                                             const bool scaleHGRotation) {
  Int return_value = 0;
  //
  //  Store pointers to the input variables
  //
  const Scalar *  cordel_ptr(cordel);
  const Scalar *  vel_ptr(vel);
  const Scalar *  elem_mass_ptr(elem_mass);
  const Scalar *  elem_dilmod_ptr(elem_dilmod);
  const Scalar *  elem_shrmod_ptr(elem_shrmod);
  const Scalar *  rotated_stretching_ptr(rotated_stretching);
  const Scalar *  spin_rate_ptr(spin_rate);
  const Scalar *  stress_new_ptr(stress_new);
  const Scalar *  rotation_ptr(rotation);

  Scalar *  rotated_stress_ptr(rotated_stress);
  Scalar *  volume_ptr(volume);
  Scalar *  elem_time_step_ptr(elem_time_step);
  Scalar *  hg_resist_old_ptr(hg_resist_old);
  Scalar *  hg_resist_new_ptr(hg_resist_new);
  Scalar *  force_new_ptr(force_new);
  Scalar *  hg_energy_ptr(hg_energy);
  Scalar *  int_energy_ptr(int_energy);
  //
  //  Create temporaries to be reused by every element
  //
  Scalar   gradop12x[24];
  Scalar * gradop12x_ptr_x(gradop12x);
  Scalar * gradop12x_ptr_y(gradop12x+8);
  Scalar * gradop12x_ptr_z(gradop12x+16);

  Scalar  cur_coords[24];
  const Scalar *const x_ptr(cur_coords);
  const Scalar *const y_ptr(cur_coords+8);
  const Scalar *const z_ptr(cur_coords+16);

  Scalar total_stress12th[6];

  Scalar hgop[32];
  const Scalar *const hgop_ptr0(hgop);
  const Scalar *const hgop_ptr1(hgop+8);
  const Scalar *const hgop_ptr2(hgop+16);
  const Scalar *const hgop_ptr3(hgop+24);
  Scalar * hgop_for_resist_calc(mid_hgop);
  if (mid_hgop == NULL){
    hgop_for_resist_calc = hgop;
  }
  //
  //  This factor will be used in hourglass computations, however is a loop constant
  //
  const Scalar fac1_pre(dt * hg_stiffness * 0.0625);

  for(Int ielem(0); ielem < nelem; ++ielem) {
    //
    //  Store local coordintes, note this changes the way coordinates are stored from
    //
    //  (x1, y1, z1, x2, y2, z2 .......)
    //    to
    //  (x1, x2, x3, x4, x5, x6, x7, x8, y1, y2, y3, .....)
    //
    transform_38_matrix_to_83_matrix(cordel_ptr, cur_coords);
    //
    //  Compute 12 times the gradient operator in three seperate steps
    //
    comp_grad12x(gradop12x, x_ptr, y_ptr, z_ptr);
    //
    // calculate 12 times the element volume, store the actual element volume
    //
    const Scalar volume12x = ddot8(x_ptr, gradop12x);
    *volume_ptr = volume12x * ONE12TH;
    //
    //  if volume is <= 0.0, report an error, otherwise continue the calculation
    //
    if(volume12x > 0.0) {
      const Scalar volinv12th(1.0/volume12x);
      //
      //  Compute the aspect ratio.  Aspect ratio is 0.5 * volume / (grad . grad)
      //  With the 12 factors this is actually 6.0 * 12 * volume / (12 * grad . 12 * grad)
      //
      const Scalar aspect = comp_aspect(gradop12x, volume12x);
      const Scalar aspect_inv = 1.0/aspect;
      //
      //  Compute the stable time step and bulk viscosity
      //
      const Scalar dtrial(sqrt(elem_mass_ptr[0] * aspect / elem_dilmod_ptr[0]));

      const Scalar traced(rotated_stretching_ptr[0] + rotated_stretching_ptr[1] + rotated_stretching_ptr[2]);

      const Scalar eps(linBulkVisc - quadBulkVisc * std::min(0.0, traced)*dtrial);

      const Scalar bulkq = eps * elem_dilmod_ptr[0] *dtrial *traced;

      const Scalar cur_time_step = dtrial * ( sqrt( 1.0 + eps * eps) - eps);
      *elem_time_step_ptr = cur_time_step;
      min_elem_time_step = std::min(min_elem_time_step, cur_time_step);
      //
      //  Rotate the stress
      //
      rotate_tensor_backward(rotation_ptr, stress_new_ptr, rotated_stress_ptr);
      //
      //  Compute the total stress (includes bulk viscosity terms)
      //
      total_stress12th[0] = ONE12TH*(rotated_stress_ptr[0] + bulkq);
      total_stress12th[1] = ONE12TH*(rotated_stress_ptr[1] + bulkq);
      total_stress12th[2] = ONE12TH*(rotated_stress_ptr[2] + bulkq);
      total_stress12th[3] = ONE12TH*(rotated_stress_ptr[3]);
      total_stress12th[4] = ONE12TH*(rotated_stress_ptr[4]);
      total_stress12th[5] = ONE12TH*(rotated_stress_ptr[5]);
      //
      //  Compute the hourglass operator. Here one 12th the volume is multiplied by
      //  12 times the gradient operators, so the true hourglass operator is computed.
      //
      comp_hgop(x_ptr, y_ptr, z_ptr,
                gradop12x_ptr_x,
                gradop12x_ptr_y,
                gradop12x_ptr_z,
                volinv12th, hgop);
      //
      //  Compute the hourglass resistance terms.  These are based on the hourglass stiffness and hourglass viscosities
      //
      const Scalar fac1 = fac1_pre * elem_shrmod_ptr[0] * aspect_inv;
      const Scalar fac2 = hg_viscosity * sqrt(elem_shrmod_ptr[0] * elem_mass_ptr[0] * aspect_inv);
      comp_force(dt,
                 spin_rate_ptr,
                 hgop_for_resist_calc,
                 hgop_ptr0, hgop_ptr1, hgop_ptr2, hgop_ptr3,
                 vel_ptr, fac1, fac2, hg_resist_old_ptr, hg_resist_new_ptr,
                 total_stress12th,
                 gradop12x_ptr_x, gradop12x_ptr_y, gradop12x_ptr_z,
                 hg_energy_ptr, int_energy_ptr, force_new_ptr, scaleHGRotation);
    } else {
      return_value = -1;
    }

    cordel_ptr += 24;
    vel_ptr += 24;
    ++volume_ptr ;
    ++elem_mass_ptr ;
    ++elem_dilmod_ptr ;
    ++elem_shrmod_ptr ;
    ++elem_time_step_ptr ;
    rotation_ptr += 9;
    stress_new_ptr += 6;
    rotated_stress_ptr += 6;
    rotated_stretching_ptr += 6;
    spin_rate_ptr += 3;
    hg_resist_old_ptr += 12;
    hg_resist_new_ptr += 12;
    ++hg_energy_ptr ;
    ++int_energy_ptr ;
    force_new_ptr += 24;
    if (mid_hgop){
      hgop_for_resist_calc += 32;
    }
  }
  return return_value;
}

} // end namespace Hex
} // end namespace APS


/**
 * This class defines the boundary conditions based on
 * a deformation gradient.  It puts an element into a
 * state of uniform strain through the following displacement field:
 *
 * \f$
 *     u_{i} = \left(F_{ij} - \delta_{ij}\right) X_{j}
 * \f$
 */
class DefGradBC
{
public:
  explicit DefGradBC(double * dg) : dgrad(dg)
  {}

  ~DefGradBC()
  {}

  void getDisplacment(double * coord, double * displ) {
    for ( int i = 1 ; i < 8 ; i++ ){
      int k = 3*i;

      displ[k  ] = ( dgrad[0] - 1.0 ) * coord[k]
	+ dgrad[1]                    * coord[k+1]
	+ dgrad[2]                    * coord[k+2];

      displ[k+1] = ( dgrad[4] - 1.0 ) * coord[k]
	+ dgrad[5]                    * coord[k+1]
	+ dgrad[3]                    * coord[k+2];

      displ[k+2] = ( dgrad[8] - 1.0 ) * coord[k]
	+ dgrad[6]                    * coord[k+1]
	+ dgrad[7]                    * coord[k+2];

    }
  }

private:
  DefGradBC& operator=(const DefGradBC& );
  DefGradBC(const DefGradBC& );
  double * dgrad;
};

void getVelocity(DefGradBC& dgbc, double dt, double * coord, double * vel) {  
  dgbc.getDisplacment(coord,vel);

  double dt_inv = (dt > 0.0) ? 1.0 / dt : 1.0;

  for ( int i = 0 ; i < 24 ; i++ ) {
    vel[i] *= dt_inv;
  }
}

namespace lame {

int
Material::initialize(
  matParams *           p)
{
  int npoints = p->nelements;
  for( int index = 0 ; index < num_state_vars * npoints ; ++index) {
    p->state_old[index] = 0.0;
    p->state_new[index] = 0.0;
  }
  return 0;
}


double
Material::getMaterialProperty(
  const std::string &   name,
  const MatProps &      props)
{

  MatProps::const_iterator itr = props.find(name);

  if (itr != props.end()){

    const double p = itr->second[0];

    return p ;

  }
  else {

    std::string msg = " Material input ";
    msg += name;
    msg += " does not exist\n";
    msg += "    Setting value of ";
    msg += name;
    msg += " to zero\n";

    int code = 1;

    reportError(code,msg);

    return -1.0;
  }
}


int
Elastic::initialize(
  matParams *            p)
{
  double youngs_modulus = properties[0];
  if ( youngs_modulus < 0.0 ) {
    std::cout << "KHP: Negative young's modulus\n";
  }

  double poissons_ratio = properties[1];
  if ( poissons_ratio < -1.0 ) {
    std::cout << "KHP: PR < -1.0\n";
  }

  if ( poissons_ratio > 0.5 ) {
    std::cout << "KHP: PR > 0.5, PR= " << poissons_ratio << "\n";
  }

  return 0;
}


int
Elastic::getStress(
  matParams *           p )
{
  double        youngs_     = properties[0];
  double        pr          = properties[1];
  double        twomu       = youngs_/(1.0+pr);
  double        alambda     = twomu*pr/(1.0-2.0*pr);
  int           num_element = p->nelements;
  double        dt          = p->dt;
  double        twomu_dt    = twomu * dt;
  double *      strain      = p->strain_rate;
  double *      stress_new  = p->stress_new;
  double *      stress_old  = p->stress_old;

  for (int i = 0; i < num_element; ++i) {
    double traced = dt * ( strain[0] + strain[1] + strain[2]);
    stress_new[0] = stress_old[0] + alambda*traced + twomu_dt*strain[0];
    stress_new[1] = stress_old[1] + alambda*traced + twomu_dt*strain[1];
    stress_new[2] = stress_old[2] + alambda*traced + twomu_dt*strain[2];
    stress_new[3] = stress_old[3] + twomu_dt*strain[3];
    stress_new[4] = stress_old[4] + twomu_dt*strain[4];
    stress_new[5] = stress_old[5] + twomu_dt*strain[5];

    strain     += 6;
    stress_new += 6;
    stress_old += 6;
  }

  return 0;
}
  
} // namespace lame


int
APSHex8ug::internalForce(  const int num_elements, 
                           const double dt,
                           double current_stable_time_step,
                           double element_time_step,
                           lame::Material  &material_model,
                           lame::matParams &materialParameters,
                           lame::MatProps  &materialProperties,
                           std::vector<double> &coordinates,
                           std::vector<double> &velocity,
                           std::vector<double> &rotation_old, std::vector<double> &rotation_new, 
                           std::vector<double> &midstep_volume, 
                           std::vector<double> &vorticity_tensor, 
                           std::vector<double> &stretch,
                           std::vector<double> &strain_rate, 
                           std::vector<double> &mid_hgop, 
                           std::vector<double> &stress_old, std::vector<double> &stress_new,
                           std::vector<double> &rotated_stress,
                           std::vector<double> &material_eff_bulk_mod,
                           std::vector<double> &material_eff_twomu,
                           std::vector<double> &shrmod,
                           std::vector<double> &dilmod,
                           std::vector<double> &element_mass,
                           std::vector<double> &force_new,
                           std::vector<double> &hourglass_energy,
                           std::vector<double> &internal_energy,
                           std::vector<double> &hg_resistance_old, std::vector<double> &hg_resistance_new
                           ) const
{

  compute_stretch(num_elements, dt, &coordinates[0], &velocity[0], &rotation_old[0], &midstep_volume[0], &vorticity_tensor[0], &rotation_new[0], &stretch[0], &strain_rate[0], &mid_hgop[0]);

  Int err = material_model.getStress(&materialParameters);

  double dilatationalHGParam = 0.05;
  double deviatoricHGParam   = 0.0;

  for(int k=0; k<num_elements; ++k) {
    shrmod[k] = material_eff_twomu[k];
    dilmod[k] = ( material_eff_bulk_mod[k] + 2.0*shrmod[k] )*(1.0/3.0);
  }

  const bool scaleHGRotation            = false; // default?
  const double linear_bulk_viscosity    = 0.0; // default?
  const double quadratic_bulk_viscosity = 0.0; // default?

  err = APS::Hex::elem_ug3dh8_mi_compute_divergence_presto(num_elements,
                                                           dt,
                                                           &coordinates[0],
                                                           &velocity[0],
                                                           &rotation_new[0],
                                                           &stress_new[0],
                                                           &strain_rate[0],
                                                           &vorticity_tensor[0],
                                                           linear_bulk_viscosity,
                                                           quadratic_bulk_viscosity,
                                                           &element_mass[0],
                                                           &dilmod[0],
                                                           &shrmod[0],
                                                           dilatationalHGParam,
                                                           deviatoricHGParam,
                                                           &rotated_stress[0],
                                                           current_stable_time_step,
                                                           &midstep_volume[0],
                                                           &element_time_step,
                                                           &hg_resistance_old[0],
                                                           &hg_resistance_new[0],
                                                           &force_new[0],
                                                           &hourglass_energy[0],
                                                           &internal_energy[0],
                                                           &mid_hgop[0],
                                                           scaleHGRotation);

  if (err) { std::cerr << "Volume error" << std::endl; }
  return err;
}

  int APSHex8ug::internalForce(const int num_elements, 
                               const double dt,
                               double current_stable_time_step,
                               double*  const element_time_step,
                               lame::Material  &material_model,
                               lame::matParams &materialParameters,
                               lame::MatProps  &/*materialProperties*/,
                               double * const coordinates,
                               double * const velocity,
                               double * const rotation_old, double*  const rotation_new, 
                               double * const midstep_volume, 
                               double * const vorticity_tensor, 
                               double * const stretch,
                               double * const strain_rate, 
                               double * const mid_hgop, 
                               double * const /*stress_old*/, double * const stress_new,
                               double * const rotated_stress,
                               double * const material_eff_bulk_mod,
                               double * const material_eff_twomu,
                               double * const shrmod,
                               double * const dilmod,
                               double * const element_mass,
                               double * const force_new,
                               double * const hourglass_energy,
                               double * const internal_energy,
                               double * const hg_resistance_old, double*  const hg_resistance_new
                               ) const
  {
    compute_stretch(num_elements, dt, coordinates, velocity, rotation_old, midstep_volume, vorticity_tensor, rotation_new, stretch, strain_rate, mid_hgop);

    Int err = material_model.getStress(&materialParameters);

    const double one_third = 1.0 / 3.0 ;
    double dilatationalHGParam = 0.05;
    double deviatoricHGParam   = 0.0;

    for(int k=0; k<num_elements; ++k) {
      shrmod[k] = material_eff_twomu[k];
      dilmod[k] = ( material_eff_bulk_mod[k] + 2.0*shrmod[k] ) * one_third ;
    }

    const bool scaleHGRotation            = false; // default?
    const double linear_bulk_viscosity    = 0.0; // default?
    const double quadratic_bulk_viscosity = 0.0; // default?

    err = APS::Hex::elem_ug3dh8_mi_compute_divergence_presto(num_elements,
                                                             dt,
                                                             coordinates,
                                                             velocity,
                                                             rotation_new,
                                                             stress_new,
                                                             strain_rate,
                                                             vorticity_tensor,
                                                             linear_bulk_viscosity,
                                                             quadratic_bulk_viscosity,
                                                             element_mass,
                                                             dilmod,
                                                             shrmod,
                                                             dilatationalHGParam,
                                                             deviatoricHGParam,
                                                             rotated_stress,
                                                             current_stable_time_step,
                                                             midstep_volume,
                                                             element_time_step,
                                                             hg_resistance_old,
                                                             hg_resistance_new,
                                                             force_new,
                                                             hourglass_energy,
                                                             internal_energy,
                                                             mid_hgop,
                                                             scaleHGRotation);

    if (err) { std::cerr << "Volume error" << std::endl; }
    return err;
  }


//--------------------------------------------------------------------
//--------------------------------------------------------------------

