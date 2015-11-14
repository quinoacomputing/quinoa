// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

//
// 7/4/2002: RAB: I was able to update this class using the
// functions in AbstractLinAlgPack_LinAlgOpPackHack.hpp.  This is wastefull in that
// I am creating temporaries every time any operation is performed
// but this was the easiest way to get things going.
//
// 7/4/2002: RAB : ToDo:  In the future it would be good to create
// some type of temporary vector server so that I could avoid
// creating all of these temporaries.  This will take some thought
// and may not be worth it for now.
//

#include <assert.h>

#include <ostream>
#include <iomanip>
#include <limits>

#include "ConstrainedOptPack_QPSchur.hpp"
#include "ConstrainedOptPack_ComputeMinMult.hpp"
#include "AbstractLinAlgPack_MatrixSymPosDefCholFactor.hpp"
#include "AbstractLinAlgPack_SpVectorOp.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOp.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "AbstractLinAlgPack_GenPermMatrixSliceOut.hpp"
#include "AbstractLinAlgPack_SpVectorOut.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsing.hpp"
#include "AbstractLinAlgPack_MatrixOpOut.hpp"
#include "AbstractLinAlgPack_VectorStdOps.hpp"
#include "AbstractLinAlgPack_EtaVector.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DVectorClassExt.hpp"
#include "DenseLinAlgPack_DVectorOp.hpp"
#include "DenseLinAlgPack_LinAlgOpPack.hpp"
#include "DenseLinAlgPack_DVectorOut.hpp"
#include "DenseLinAlgPack_DMatrixOut.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_Assert.hpp"

namespace LinAlgOpPack {
using AbstractLinAlgPack::Vp_StV;
using AbstractLinAlgPack::Vp_StMtV;
using AbstractLinAlgPack::Vp_StV;
using AbstractLinAlgPack::Vp_StMtV;
}

namespace {

// Some local helper functions.

template< class T >
inline
T my_max( const T& v1, const T& v2 ) { return v1 > v2 ? v1 : v2; }

//
// Print a bnd as a string
//
inline
const char* bnd_str( ConstrainedOptPack::EBounds bnd )
{
  switch(bnd) {
    case ConstrainedOptPack::FREE:
      return "FREE";
    case ConstrainedOptPack::UPPER:
      return "UPPER";
    case ConstrainedOptPack::LOWER:
      return "LOWER";
    case ConstrainedOptPack::EQUALITY:
      return "EQUALITY";
  }
  TEUCHOS_TEST_FOR_EXCEPT(true);	// should never be executed
  return 0;
}

//
// print a bool
//
inline
const char* bool_str( bool b )
{
  return b ? "true" : "false";
}

//
// Return a std::string that has a file name, line number and
// error message.
//
std::string error_msg(
  const char file_name[], const int line_num, const char err_msg[]
  )
{
  std::ostringstream  omsg;
  omsg
    << file_name << ":" << line_num << ":" << err_msg;
  return omsg.str();
}

//
// Deincrement all indices less that k_remove
//
void deincrement_indices(
  DenseLinAlgPack::size_type                k_remove
  ,std::vector<DenseLinAlgPack::size_type>  *indice_vector
  ,size_t                              len_vector
  )
{
  typedef DenseLinAlgPack::size_type				size_type;
  typedef std::vector<DenseLinAlgPack::size_type>	vec_t;
  TEUCHOS_TEST_FOR_EXCEPT( !(  len_vector <= indice_vector->size()  ) );
  for( vec_t::iterator itr = indice_vector->begin(); itr != indice_vector->begin() + len_vector; ++itr ) {
    if( *itr > k_remove )
      --(*itr);
  }
}

//
// Insert the element (r_v,c_v) into r[] and c[] sorted by r[]
//
void insert_pair_sorted(
  DenseLinAlgPack::size_type                r_v
  ,DenseLinAlgPack::size_type               c_v
  ,size_t                              len_vector  // length of the new vector
  ,std::vector<DenseLinAlgPack::size_type>  *r
  ,std::vector<DenseLinAlgPack::size_type>  *c
  )
{
  typedef std::vector<DenseLinAlgPack::size_type> rc_t;
  TEUCHOS_TEST_FOR_EXCEPT( !(  r->size() >= len_vector && c->size() >= len_vector  ) );
  // find the insertion point in r[]
  rc_t::iterator
    itr = std::lower_bound( r->begin(), r->begin() + len_vector-1, r_v );
  const DenseLinAlgPack::size_type p = itr - r->begin();
  // Shift all of the stuff out of the way to make room for the insert
  {for( rc_t::iterator itr_last = r->begin() + len_vector-1;
      itr_last > r->begin() + p; --itr_last )
  {
    *itr_last = *(itr_last-1);
  }}
  {for( rc_t::iterator itr_last = c->begin() + len_vector-1;
      itr_last > c->begin() + p; --itr_last )
  {
    *itr_last = *(itr_last-1);
  }}
  // Insert the new elements
  (*r)[p] = r_v;
  (*c)[p] = c_v;
}

//
// z_hat = inv(S_hat) * ( d_hat - U_hat'*vo )
//
void calc_z(
  const AbstractLinAlgPack::MatrixSymOpNonsing   &S_hat
  ,const DenseLinAlgPack::DVectorSlice           &d_hat
  ,const AbstractLinAlgPack::MatrixOp            &U_hat
  ,const DenseLinAlgPack::DVectorSlice           *vo       // If NULL then assumed zero
  ,DenseLinAlgPack::DVectorSlice                 *z_hat
  )
{
  using LinAlgOpPack::Vp_StMtV;
  using LinAlgOpPack::V_InvMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  Workspace<DenseLinAlgPack::value_type> t_ws(wss,d_hat.dim());
  DenseLinAlgPack::DVectorSlice                t(&t_ws[0],t_ws.size());
  t = d_hat;
  if(vo)
    Vp_StMtV( &t, -1.0, U_hat, BLAS_Cpp::trans, *vo );
  V_InvMtV( z_hat, S_hat, BLAS_Cpp::no_trans, t );
}

//
// v = inv(Ko) * ( fo - U_hat * z_hat )
//
void calc_v(
  const AbstractLinAlgPack::MatrixSymOpNonsing   &Ko
  ,const DenseLinAlgPack::DVectorSlice                         *fo    // If NULL then assumed to be zero
  ,const AbstractLinAlgPack::MatrixOp                &U_hat
  ,const DenseLinAlgPack::DVectorSlice                         &z_hat // Only accessed if U_hat.cols() > 0
  ,DenseLinAlgPack::DVectorSlice                               *v
  )
{
  using DenseLinAlgPack::norm_inf;
  using LinAlgOpPack::Vp_StMtV;
  using LinAlgOpPack::V_InvMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  Workspace<DenseLinAlgPack::value_type> t_ws(wss,v->dim());
  DenseLinAlgPack::DVectorSlice                t(&t_ws[0],t_ws.size());
  if(fo) {	
    t = *fo;
  }
  else {
    t = 0.0;
  }
  if( U_hat.cols() )
    Vp_StMtV( &t, -1.0, U_hat, BLAS_Cpp::no_trans, z_hat );
  if( norm_inf(t) > 0.0 )
    V_InvMtV( v, Ko, BLAS_Cpp::no_trans, t );
  else
    *v = 0.0;
}

//
// mu_D_hat =
// 		- Q_XD_hat' * g
// 		- Q_XD_hat' * G * x
// 		- Q_XD_hat' * A * v(n_R+1:n_R+m)
// 		- Q_XD_hat' * A_bar * P_plus_hat * z_hat
//
void calc_mu_D(
  const ConstrainedOptPack::QPSchur::ActiveSet  &act_set
  ,const DenseLinAlgPack::DVectorSlice                         &x
  ,const DenseLinAlgPack::DVectorSlice                         &v
  ,DenseLinAlgPack::DVectorSlice                               *mu_D
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::Vp_StPtMtV;
  using AbstractLinAlgPack::V_MtV;
  using AbstractLinAlgPack::Vp_MtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const ConstrainedOptPack::QPSchurPack::QP
    &qp = act_set.qp();
  const DenseLinAlgPack::size_type
    n = qp.n(),
    n_R = qp.n_R(),
    m = qp.m();

  const AbstractLinAlgPack::GenPermMatrixSlice &Q_XD_hat = act_set.Q_XD_hat();
  const DenseLinAlgPack::DVectorSlice			 g = qp.g();
  const AbstractLinAlgPack::MatrixSymOp &G = qp.G();
  // mu_D_hat = - Q_XD_hat' * g
  V_StMtV( mu_D, -1.0, Q_XD_hat, trans, g ); 
  // mu_D_hat += - Q_XD_hat' * G * x
  Vp_StPtMtV( mu_D, -1.0, Q_XD_hat, trans, G, no_trans, x );
  // mu_D_hat += - Q_XD_hat' * A * v(n_R+1:n_R+m)
  if( m ) {
    Vp_StPtMtV( mu_D, -1.0, Q_XD_hat, trans, qp.A(), no_trans, v(n_R+1,n_R+m) );
  }
  // p_mu_D_hat += - Q_XD_hat' * A_bar * P_plus_hat * z_hat
  if( act_set.q_plus_hat() && act_set.q_hat() ) {
    const DenseLinAlgPack::DVectorSlice z_hat = act_set.z_hat();
    AbstractLinAlgPack::SpVector P_plus_hat_z_hat;
    V_MtV( &P_plus_hat_z_hat, act_set.P_plus_hat(), no_trans, z_hat ); 
    Vp_StPtMtV( mu_D, -1.0, Q_XD_hat, trans
      , qp.constraints().A_bar(), no_trans, P_plus_hat_z_hat() );
  }
}

//
// p_mu_D_hat =
// 		- Q_XD_hat' * G * Q_R * p_v(1:n_R)
// 		- Q_XD_hat' * G * P_XF_hat * p_z_hat
// 		- Q_XD_hat' * A * p_v(n_R+1:n_R+m)
// 		- Q_XD_hat' * A_bar * (P_plus_hat * p_z_hat + e(ja))
//
void calc_p_mu_D(
  const ConstrainedOptPack::QPSchur::ActiveSet  &act_set
  ,const DenseLinAlgPack::DVectorSlice                         &p_v
  ,const DenseLinAlgPack::DVectorSlice                         &p_z_hat
  ,const DenseLinAlgPack::size_type                           *ja       // If != NULL then we will include the term e(ja)
  ,DenseLinAlgPack::DVectorSlice                               *p_mu_D
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using AbstractLinAlgPack::V_MtV;
  using AbstractLinAlgPack::Vp_MtV;
  using LinAlgOpPack::Vp_StPtMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const ConstrainedOptPack::QPSchurPack::QP
    &qp = act_set.qp();
  const ConstrainedOptPack::QPSchurPack::Constraints
    &constraints = qp.constraints();
  const DenseLinAlgPack::size_type
    n = qp.n(),
    n_R = qp.n_R(),
    m = qp.m();

  const AbstractLinAlgPack::GenPermMatrixSlice &Q_XD_hat = act_set.Q_XD_hat();
  const AbstractLinAlgPack::MatrixSymOp &G = qp.G();
  // p_mu_D_hat = - Q_XD_hat' * G * Q_R * p_v(1:n_R)
  {
    AbstractLinAlgPack::SpVector Q_R_p_v1;
    V_MtV( &Q_R_p_v1, qp.Q_R(), no_trans, p_v(1,n_R) ); 
    Vp_StPtMtV( p_mu_D, -1.0, Q_XD_hat, trans, G, no_trans, Q_R_p_v1(), 0.0 );
  }
  // p_mu_D_hat += - Q_XD_hat' * G * P_XF_hat * p_z_hat
  if( act_set.q_F_hat() ) {
    AbstractLinAlgPack::SpVector P_XF_hat_p_z_hat;
    V_MtV( &P_XF_hat_p_z_hat, act_set.P_XF_hat(), no_trans, p_z_hat ); 
    Vp_StPtMtV( p_mu_D, -1.0, Q_XD_hat, trans, G, no_trans, P_XF_hat_p_z_hat() );
  }
  // p_mu_D_hat += - Q_XD_hat' * A * p_v(n_R+1:n_R+m)
  if( m ) {
    Vp_StPtMtV( p_mu_D, -1.0, Q_XD_hat, trans, qp.A(), no_trans, p_v(n_R+1,n_R+m) );
  }
  // p_mu_D_hat += - Q_XD_hat' * A_bar * ( P_plus_hat * p_z_hat + e(ja) )
  if( act_set.q_plus_hat() || ja ) {
    AbstractLinAlgPack::SpVector p_lambda_bar(
      n+constraints.m_breve(), act_set.q_plus_hat() + (ja ? 1 : 0) );
    if( act_set.q_plus_hat() ) // p_lambda_bar =  P_plus_hat * p_z_hat
      Vp_MtV( &p_lambda_bar, act_set.P_plus_hat(), no_trans, p_z_hat );
    if( ja ) // p_lambda_bar += e(ja) (non-duplicate indices?)
      p_lambda_bar.insert_element(AbstractLinAlgPack::SpVector::element_type(*ja,1.0));
    Vp_StPtMtV( p_mu_D, -1.0, Q_XD_hat, trans
      , constraints.A_bar(), no_trans, p_lambda_bar() );
  }
}

//
// Calculate the residual of the augmented KKT system:
//
// [ ro ] = [   Ko     U_hat ] [   v   ] + [  ao * bo ]
// [ ra ]   [ U_hat'   V_hat ] [ z_hat ]   [  aa * ba ]
//
// Expanding this out we have:
//
// ro = Ko * v + U_hat * z_hat + ao * bo
//
//    = [ Q_R'*G*Q_R   Q_R'*A ] * [   x_R  ] + [ Q_R'*G*P_XF_hat + Q_R'*A_bar*P_plus_hat ] * z_hat + [ ao*boR ]
//      [  A'*Q_R             ]   [ lambda ]   [          A'*P_XF_hat                    ]           [ ao*bom ]
//
//    = [ Q_R'*G*(Q_R*x_R + P_XF_hat*z_hat) + Q_R'*A*lambda + Q_R'*A_bar*P_plus_hat*z_hat + ao*boR ]   
//      [      A'*(Q_R*x_R + P_XF_hat*z_hat) + ao*bom                                              ]
//
// ra = [ U_hat' * v + V_hat * z_hat + aa*ba
//
//    = [ P_XF_hat'*G*Q_R + P_plus_hat'*A_bar'*Q_R , P_XF_hat'*A ] * [ x_R ; lambda ]
//      + [ P_XF_hat'*G*P_XF_hat + P_XF_hat'*A_bar*P_plus_hat + P_plus_hat'*A_bar'*P_XF_hat
//          + P_F_tilde'*P_C_hat + P_C_hat'*P_F_tilde ] * z_hat + aa*ba
//
//    = P_XF_hat'*G*(Q_R*x_R + P_XF_hat*z_hat) + P_plus_hat'*A_bar'*(Q_R*x_R + P_XF_hat*z_hat)
//      + P_XF_hat'*A*lambda + P_XF_hat'*A_bar*P_plus_hat*z_hat
//      + (P_F_tilde'*P_C_hat + P_C_hat'*P_F_tilde)*z_hat + aa*ba
//
// Given the QP matrices G, A, and A_bar and the active set mapping matrices Q_R, P_XF_hat,
// P_plus_hat and P_FC_hat = P_F_tilde'*P_C_hat, we can compute the residual efficiently
// as follows:
//
// x_free = Q_R*x_R + P_XF_hat*z_hat (sparse)
// 
// lambda_bar = P_plus_hat*z_hat (sparse)
//
// t1 = G*x_free
//
// t2 = A*lambda
//
// t3 = A_bar*lambda_bar
//
// roR = Q_R'*t1 + Q_R'*t2 + Q_R'*t3 + ao*boR
//
// rom = A'*t1 + ao*bom
//
// ra = P_XF_hat'*t1 + P_plus_hat'*A_bar'*x_free + P_XF_hat'*t2 + P_XF_hat'*t3
//      + (P_FC_hat + P_FC_hat')*z_hat  + aa*ba
//
// On output we will have set:
//
//   roR_scaling = ||Q_R'*t1||inf + ||Q_R'*t2||inf + ||Q_R'*t3||inf + ||ao*boR||inf
//
//   rom_scaling = ||A'*t1||inf + ||ao*bom||inf
//
//   ra_scaling  = ||P_XF_hat'*t1||inf + ||P_plus_hat'*A_bar'*t1||inf + ||P_XF_hat'*t2||inf
//                 + ||P_XF_hat'*t3||inf + ||(P_FC_hat + P_FC_hat')*z_hat|| + ||aa*ba||inf
// 
// Note: In the future we could make this a little more efficent by combining (Q_R + P_XF_hat) into
// an single permulation matrix and then we could leave out the terms for the variables initially 
// fixed and still fixed rather than computing the terms then throwing them away.
//
// Also, in the future, this could really be implemented for extended precision data types but
// it would require some real work!
//
template<class val_type>
void calc_resid(
  const ConstrainedOptPack::QPSchur::ActiveSet     &act_set
  ,const DenseLinAlgPack::DVectorSlice                            &v
  ,const DenseLinAlgPack::DVectorSlice                            &z_hat        // Only accessed if q_hat > 0
  ,const DenseLinAlgPack::value_type                             ao            // Only accessed if bo != NULL
  ,const DenseLinAlgPack::DVectorSlice                            *bo           // If NULL then considered 0
  ,DenseLinAlgPack::VectorSliceTmpl<val_type>                    *ro
  ,DenseLinAlgPack::value_type                                   *roR_scaling
  ,DenseLinAlgPack::value_type                                   *rom_scaling  // Only set if m > 0
  ,const DenseLinAlgPack::value_type                             aa            // Only accessed if q_hat > 0
  ,const DenseLinAlgPack::DVectorSlice                            *ba           // If NULL then considered 0, Only accessed if q_hat > 0
  ,DenseLinAlgPack::VectorSliceTmpl<val_type>                    *ra           // Only set if q_hat > 0
  ,DenseLinAlgPack::value_type                                   *ra_scaling   // Only set if q_hat > 0
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using DenseLinAlgPack::norm_inf;
  using DenseLinAlgPack::DVectorSlice;
  using DenseLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::V_MtV;
  using AbstractLinAlgPack::SpVector;
  using AbstractLinAlgPack::GenPermMatrixSlice;
  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_MtV;
  using LinAlgOpPack::Vp_StMtV;
  using LinAlgOpPack::Vp_StPtMtV;
  namespace COP = ConstrainedOptPack;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  
  const COP::QPSchurPack::QP
    &qp = act_set.qp();
  const COP::QPSchurPack::Constraints
    &constraints = qp.constraints();
  const GenPermMatrixSlice
    &Q_R        = qp.Q_R(),
    &P_XF_hat   = act_set.P_XF_hat(),
    &P_plus_hat = act_set.P_plus_hat();
  const DenseLinAlgPack::size_type
    n          = qp.n(),
    n_R        = qp.n_R(),
    m          = qp.m(),
    m_breve    = constraints.m_breve(),
    q_hat      = act_set.q_hat(),
    q_F_hat    = act_set.q_F_hat(),
    q_C_hat    = act_set.q_C_hat(),
    q_plus_hat = act_set.q_plus_hat();
  
  const DVectorSlice
    x_R = v(1,n_R),
    lambda = ( m ? v(n_R+1,n_R+m): DVectorSlice() ),
    boR = ( bo ? (*bo)(1,n_R) : DVectorSlice() ),
    bom = ( bo && m ? (*bo)(n_R+1,n_R+m) : DVectorSlice() );

  DenseLinAlgPack::VectorSliceTmpl<val_type>
    roR = (*ro)(1,n_R),
    rom = ( m ? (*ro)(n_R+1,n_R+m) : DenseLinAlgPack::VectorSliceTmpl<val_type>() );

  Workspace<DenseLinAlgPack::value_type>
    x_free_ws(wss,n);
  DenseLinAlgPack::DVectorSlice
    x_free(&x_free_ws[0],x_free_ws.size());
  Workspace<val_type>
    t1_ws(wss,n),
    t2_ws(wss,n),
    t3_ws(wss,n),
    tR_ws(wss,n_R),
    tm_ws(wss,m),
    ta_ws(wss,q_hat);
  DenseLinAlgPack::VectorSliceTmpl<val_type>
    t1(&t1_ws[0],t1_ws.size()),
    t2(&t2_ws[0],t2_ws.size()),
    t3(&t3_ws[0],t3_ws.size()),
    tR(&tR_ws[0],tR_ws.size()),
    tm(tm_ws.size()?&tm_ws[0]:NULL,tm_ws.size()),
    ta(ta_ws.size()?&ta_ws[0]:NULL,ta_ws.size());

  *roR_scaling = 0.0;
  if( m )
    *rom_scaling = 0.0;
  if( q_hat )
    *ra_scaling  = 0.0;
  
  // x_free = Q_R*x_R + P_XF_hat*z_hat (dense for now)
  LinAlgOpPack::V_MtV( &x_free, Q_R, no_trans, x_R );
  if( q_F_hat )
    LinAlgOpPack::Vp_MtV( &x_free, P_XF_hat, no_trans, z_hat );
  // lambda_bar = P_plus_hat*z_hat
  SpVector lambda_bar;
  if( q_plus_hat )
    V_MtV( &lambda_bar, P_plus_hat, no_trans, z_hat );
  // t1 = G*x_free
  Vp_StMtV( &t1, 1.0, qp.G(), no_trans, x_free, 0.0 );
  // t2 = A*lambda
  if( m )
    Vp_StMtV( &t2, 1.0, qp.A(), no_trans, lambda, 0.0 );
  // t3 = A_bar*lambda_bar
  if( q_plus_hat )
    Vp_StMtV( &t3, 1.0, constraints.A_bar(), no_trans, lambda_bar(), 0.0 );
  // roR = Q_R'*t1 + Q_R'*t2 + Q_R'*t3 + ao*boR
  LinAlgOpPack::V_MtV( &tR, Q_R, trans, t1 );  // roR = Q_R'*t1
  *roR_scaling += norm_inf(tR);
  roR = tR;
  if( m ) {                     // roR += Q_R'*t2
    LinAlgOpPack::V_MtV( &tR, Q_R, trans, t2 );
    *roR_scaling += norm_inf(tR);
    LinAlgOpPack::Vp_V(&roR,tR);
  }
  if( q_plus_hat ) {           // roR += Q_R'*t3
    LinAlgOpPack::V_MtV( &tR, Q_R, trans, t3 );
    *roR_scaling += norm_inf(tR);
    LinAlgOpPack::Vp_V(&roR,tR);
  }
  if( bo ) {
    LinAlgOpPack::Vp_StV( &roR, ao, boR );    // roR += ao*boR
    *roR_scaling += std::fabs(ao) * norm_inf(boR);
  }
  // rom = A'*t1 + ao*bom
  if( m ) {
    Vp_StMtV( &tm, 1.0, qp.A(), trans, t1, 0.0 );   // A'*t1
    *rom_scaling += norm_inf(tm);
    LinAlgOpPack::Vp_V(&rom,tm);
    if(bo) {
      LinAlgOpPack::Vp_StV( &rom, ao, bom );        // rom += ao*bom
      *rom_scaling += std::fabs(ao)*norm_inf(bom);
    }
  }
  // ra = P_XF_hat'*t1 + P_plus_hat'*A_bar'*x_free + P_XF_hat'*t2 + P_XF_hat'*t3
  //      +(P_FC_hat + P_FC_hat')*z_hat + aa*ba
  if( q_hat ) {
    if(ba) {              // ra = aa*ba
      V_StV( ra, aa, *ba );
      *ra_scaling += std::fabs(aa) * norm_inf(*ba);
    }
    else {
      *ra = 0.0;
    }
    if( q_F_hat ) {       // ra +=  P_XF_hat'*t1
      LinAlgOpPack::V_MtV( &ta, P_XF_hat, trans, t1 );
      *ra_scaling += norm_inf(ta);
      LinAlgOpPack::Vp_V(ra,ta);
    }
    if( q_plus_hat ) {    // ra += P_plus_hat'*A_bar'*x_free
      Vp_StPtMtV( &ta, 1.0, P_plus_hat, trans, constraints.A_bar(), trans, x_free, 0.0 );
      *ra_scaling += norm_inf(ta);
      LinAlgOpPack::Vp_V(ra,ta);
    }
    if( q_F_hat && m ) {  // ra += P_XF_hat'*t2
      LinAlgOpPack::V_MtV( &ta, P_XF_hat, trans, t2 );
      *ra_scaling += norm_inf(ta);
      LinAlgOpPack::Vp_V(ra,ta);
    }
    if( q_F_hat && q_plus_hat ) { // ra += P_XF_hat'*t3
      LinAlgOpPack::V_MtV( &ta, P_XF_hat, trans, t3 );
      *ra_scaling += norm_inf(ta);
      LinAlgOpPack::Vp_V(ra,ta);
    }
    if( q_C_hat ) {       // ra += (P_FC_hat + P_FC_hat')*z_hat
      const GenPermMatrixSlice
        &P_FC_hat = act_set.P_FC_hat();
      ta = 0.0;
      for( GenPermMatrixSlice::const_iterator itr = P_FC_hat.begin(); itr != P_FC_hat.end(); ++itr ) {
        ta(itr->row_i()) = z_hat(itr->col_j());
        ta(itr->col_j()) = z_hat(itr->row_i());
      }
      *ra_scaling += norm_inf(ta);
      LinAlgOpPack::Vp_V(ra,ta);
    }
  }
}

//
// Correct a nearly degenerate Lagrange multiplier
//
// If < 0 is returned it means that the multiplier could not
// be corrected and this should be considered to be dual infeasible
// In this case the error output is sent to *out if print_dual_infeas
// == true.
//
// If 0 is returned then the multiplier was near degenerate and
// was corrected.  In this case a warning message is printed to
// *out.
//
// If > 0 is returned then the multiplier's sign
// is just fine and no corrections were needed (no output).
//
int correct_dual_infeas(
  const DenseLinAlgPack::size_type                                  j                // for output info only
  ,const ConstrainedOptPack::EBounds                  bnd_j
  ,const DenseLinAlgPack::value_type                                t_P              // (> 0) full step length
  ,const DenseLinAlgPack::value_type                                scale            // (> 0) scaling value
  ,const DenseLinAlgPack::value_type                                dual_infeas_tol
  ,const DenseLinAlgPack::value_type                                degen_mult_val
  ,std::ostream                                                *out             // Can be NULL
  ,const ConstrainedOptPack::QPSchur::EOutputLevel    olevel
  ,const bool                                                  print_dual_infeas
  ,const char                                                  nu_j_n[]         // Name of nu_j
  ,DenseLinAlgPack::value_type                                      *nu_j            // required
  ,DenseLinAlgPack::value_type                                      *scaled_viol     // = scale*nu_j*(bnd_j==UPPER ? 1.0: -1.0 ) (after output)
  ,const char                                                  p_nu_j_n[]    = NULL // Name of p_nu_j (can be NULL if p_nu_j==NULL)
  ,DenseLinAlgPack::value_type                                      *p_nu_j       = NULL // optional (can be NULL)
  ,const char                                                  nu_j_plus_n[] = NULL // Name of nu_j_plus (can be NULL if p_nu_j==NULL)
  ,DenseLinAlgPack::value_type                                      *nu_j_plus    = NULL // optional (can be NULL)
  )
{
  typedef DenseLinAlgPack::value_type value_type;
  namespace COP = ConstrainedOptPack;

  value_type nu_j_max = (*scaled_viol) = scale * (*nu_j) * (bnd_j == COP::UPPER ? +1.0 : -1.0);
  if( nu_j_max > 0.0 || bnd_j == COP::EQUALITY ) // Leave any multiplier value with the correct sign alone!
    return +1; // No correction needed
  // See if we need to correct the multiplier
  nu_j_max = std::fabs(nu_j_max);
  if( nu_j_max < dual_infeas_tol ) {
    // This is a near degenerate multiplier so adjust it
    value_type degen_val = degen_mult_val * ( bnd_j == COP::UPPER ? +1.0 : -1.0 );
    if( out && (int)olevel >= (int)COP::QPSchur::OUTPUT_BASIC_INFO ) {
      *out
        << "\nWarning, the constriant a(" << j << ") currently at its "
        << (bnd_j == COP::UPPER ? "UPPER" : "LOWER") << " bound"
        << " has the wrong Lagrange multiplier value but\n"
        << "scale*|"<<nu_j_n<<"| = " << scale << " * |" << (*nu_j)
        << "| = " << nu_j_max  << " < dual_infeas_tol = " << dual_infeas_tol
        << "\nTherefore, this is considered a degenerate constraint and this "
        << "multiplier is set to " << degen_val << std::endl;
    }
    if(p_nu_j) {
      nu_j_max += std::fabs( t_P * (*p_nu_j) ) * scale;
      if( nu_j_max < dual_infeas_tol ) {
        // The full step is also degenerate so adjust it also
        if( out && (int)olevel >= (int)COP::QPSchur::OUTPUT_BASIC_INFO ) {
          *out
            << "Also, the maximum full step scale*(|"<<nu_j_n<<"|+|t_P*"<<p_nu_j_n<<"|) = "
            << scale << " * (|" << (*nu_j) << "| +  |" << t_P << " * " << (*p_nu_j) << "|) = "
            << nu_j_max << " < dual_infeas_tol = " << dual_infeas_tol
            << "\nTherefore, this is considered degenerate and therefore "
            << "seting "<<p_nu_j_n<<" = 0";
          if(nu_j_plus)
            *out << " and "<< nu_j_plus_n <<" = " << degen_val;
          *out << std::endl;
        }
        *p_nu_j = 0.0; // Don't let it limit the step length
        if(nu_j_plus) {
          *nu_j_plus = degen_val;
        }
      }
    }
    *nu_j = degen_val;  // Now set it
    *scaled_viol = scale * (*nu_j) * (bnd_j == COP::UPPER ? +1.0 : -1.0); // Not violated!
    return 0;
  }
  else {
    if( print_dual_infeas && out && (int)olevel >= (int)COP::QPSchur::OUTPUT_BASIC_INFO ) {
      *out
        << "\nError, the constriant a(" << j << ") currently at its "
        << (bnd_j == COP::UPPER ? "UPPER" : "LOWER") << " bound"
        << " has the wrong Lagrange multiplier value and\n"
        << "scale*|"<<nu_j_n<<"| = " << scale << " * |" << (*nu_j)
        << "| = " << nu_j_max  << " > dual_infeas_tol = " << dual_infeas_tol
        << "\nThis is an indication of instability in the calculations.\n"
        << "The QP algorithm is terminated!\n";
    }
    return -1;
  }
  return 0; // Will never be executed!
}

//
// Calculate the QP objective if it has not been already
//
// qp_grad_norm_inf = ||g + G*x||inf
//
// On input, if qp_grad_norm_inf != 0, then it will be assumed
// that this value has already been computed and the computation will
// be skipped.
//
void calc_obj_grad_norm_inf(
  const ConstrainedOptPack::QPSchurPack::QP     &qp
  ,const DenseLinAlgPack::DVectorSlice                         &x
  ,DenseLinAlgPack::value_type                                *qp_grad_norm_inf
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement this?
}

}	// end namespace

namespace ConstrainedOptPack {

// public member functions for QPSchur::U_hat_t

QPSchur::U_hat_t::U_hat_t()
  :G_(NULL)
  ,A_(NULL)
  ,A_bar_(NULL)
  ,Q_R_(NULL)
  ,P_XF_hat_(NULL)
  ,P_plus_hat_(NULL)
{}

void QPSchur::U_hat_t::initialize( 
  const MatrixSymOp       *G
  ,const MatrixOp         *A
  ,const MatrixOp         *A_bar
  ,const GenPermMatrixSlice   *Q_R
  ,const GenPermMatrixSlice   *P_XF_hat
  ,const GenPermMatrixSlice   *P_plus_hat
  )
{
  G_            = G;
  A_            = A;
  A_bar_        = A_bar;
  Q_R_          = Q_R;
  P_XF_hat_     = P_XF_hat;
  P_plus_hat_   = P_plus_hat;
}

size_type QPSchur::U_hat_t::rows() const
{
  return Q_R_->cols() + ( A_ ? A_->cols() : 0 );
}

size_type QPSchur::U_hat_t::cols() const
{
  return P_plus_hat_->cols();
}

/* 10/25/00: I don't think we need this function!
void QPSchur::U_hat_t::Mp_StM(DMatrixSlice* C, value_type a
  , BLAS_Cpp::Transp M_trans ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using LinAlgOpPack::Mp_StMtP;
  using LinAlgOpPack::Mp_StPtMtP;

  // C += a * op(U_hat)

  LinAlgOpPack::Mp_M_assert_sizes( C->rows(), C->cols(), no_trans
    , rows(), cols(), M_trans );

  const size_type
    n_R	= Q_R_->cols(),
    m 	= A() ? A()->cols()  : 0;

  if( M_trans == no_trans ) {
    //
    // C += a * op(U_hat)
    // 
    //    = a * [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ]
    //          [                  A' * P_XF_hat                   ]
    //          
    // C1 += a * Q_R' * G * P_XF_hat + a * Q_R' * A_bar * P_plus_hat
    // 
    // C2 += a * A' * P_XF_hat
    //
    DMatrixSlice
      C1 = (*C)(1,n_R,1,C->cols()),
      C2 = m ? (*C)(n_R+1,n_R+m,1,C->cols()) : DMatrixSlice();
    // C1 += a * Q_R' * G * P_XF_hat
    if( P_XF_hat().nz() )
      Mp_StPtMtP( &C1, a, Q_R(), trans, G(), no_trans, P_XF_hat(), no_trans );
    // C1 += a * Q_R' * A_bar * P_plus_hat
    if( P_plus_hat().nz() )
      Mp_StPtMtP( &C1, a, Q_R(), trans, A_bar(), no_trans, P_plus_hat(), no_trans );
    // C2 += a * A' * P_XF_hat
    if(m && P_XF_hat().nz())
      Mp_StMtP( &C2, a, *A(), trans, P_plus_hat(), no_trans );
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true);	// Implement this!
  }
}
*/

void QPSchur::U_hat_t::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const DVectorSlice& x, value_type b
  ) const
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using DenseLinAlgPack::Vt_S;
  using AbstractLinAlgPack::V_MtV;
  using AbstractLinAlgPack::Vp_StMtV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_StMtV;
  using LinAlgOpPack::Vp_StPtMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  LinAlgOpPack::Vp_MtV_assert_sizes(y->dim(),rows(),cols(),M_trans,x.dim());

  //
  // U_hat = [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ]
  //         [                  A' * P_XF_hat                   ]
  //         

  const size_type
    n_R	= Q_R_->cols(),
    m 	= A() ? A()->cols()  : 0;

  if( M_trans == BLAS_Cpp::no_trans ) {
    //
    // y =  b*y + a * U_hat * x
    // 
    //   = b*y + a * [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ] * x
    //               [                  A' * P_XF_hat                   ]
    // 
    //  =>
    // 
    // y1 = b * y1 + a * Q_R' * G * P_XF_hat * x + a * Q_R' * A_bar * P_plus_hat * x
    // y2 = b * y2 + a * A' * P_XF_hat * x
    // 
    DVectorSlice
      y1 = (*y)(1,n_R),
      y2 = m ? (*y)(n_R+1,n_R+m) : DVectorSlice();
    SpVector
      P_XF_hat_x,
      P_plus_hat_x;
    // P_XF_hat_x = P_XF_hat * x
    if( P_XF_hat().nz() )
      V_MtV( &P_XF_hat_x, P_XF_hat(), no_trans, x );
    // P_plus_hat_x = P_plus_hat * x
    if(P_plus_hat().nz())
      V_MtV( &P_plus_hat_x, P_plus_hat(), no_trans, x );
    // y1 = b * y1
    if(b==0.0)      y1=0.0;
    else if(b!=1.0) Vt_S(&y1,b);
    // y1 += a * Q_R' * G * P_XF_hat_x
    if(P_XF_hat().nz())
      Vp_StPtMtV( &y1, a, Q_R(), trans, G(), no_trans, P_XF_hat_x() );
    // y1 += a * Q_R' * A_bar * P_plus_hat_x
    if(P_plus_hat().nz())
      Vp_StPtMtV( &y1, a, Q_R(), trans, A_bar(), no_trans, P_plus_hat_x() );
    if(m) {
      // y2 = b * y2
      if(b==0.0)      y2=0.0;
      else if(b!=1.0) Vt_S(&y2,b);
      // y2 +=  a * A' * P_XF_hat_x
      if( P_XF_hat().nz() )
        Vp_StMtV( &y2, a, *A(), trans, P_XF_hat_x() );
    }
  }
  else if( M_trans == BLAS_Cpp::trans ) {
    //
    // y =  b*y + a * U_hat' * x
    // 
    //   = b*y + a * [  P_XF_hat' * G * Q_R + P_plus_hat' * A_bar' * Q_R, P_XF_hat' * A ] * [ x1 ]
    //                                                                                      [ x2 ]
    //  =>
    // 
    // y = b * y + a * P_XF_hat' * G * Q_R * x1 + a * P_plus_hat' * A_bar' * Q_R * x1
    //     + a * P_XF_hat' * A * x2
    // 
    const DVectorSlice
      x1 = x(1,n_R),
      x2 = m ? x(n_R+1,n_R+m) : DVectorSlice();
    SpVector
      Q_R_x1;
    // Q_R_x1 = Q_R * x1
    V_MtV( &Q_R_x1, Q_R(), no_trans, x1 );
    // y = b*y
    if(b==0.0)      *y = 0.0;
    else if(b!=1.0) Vt_S( y, b );
    // y += a * P_XF_hat' * G * Q_R_x1
    if(P_XF_hat().nz())
      Vp_StPtMtV( y, a, P_XF_hat(), trans, G(), no_trans, Q_R_x1() );
    // y += a * P_plus_hat' * A_bar' * Q_R_x1
    if(P_plus_hat().nz())
      Vp_StPtMtV( y, a, P_plus_hat(), trans, A_bar(), trans, Q_R_x1() );
    // y += a * P_XF_hat' * A * x2
    if( m && P_XF_hat().nz() )
      Vp_StPtMtV( y, a, P_XF_hat(), trans, *A(), no_trans, x2 );
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true);	// Invalid value for M_trans
  }
}

void QPSchur::U_hat_t::Vp_StMtV(
  DVectorSlice* y, value_type a, BLAS_Cpp::Transp M_trans
  ,const SpVectorSlice& x, value_type b
  ) const
{
//	// Uncomment to use the default version
//	MatrixOp::Vp_StMtV(y,a,M_trans,x,b); return;

  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using DenseLinAlgPack::Vt_S;
  using LinAlgOpPack::V_MtV;
  using AbstractLinAlgPack::V_MtV;
  using AbstractLinAlgPack::Vp_StMtV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_StMtV;
  using LinAlgOpPack::Vp_StPtMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  LinAlgOpPack::Vp_MtV_assert_sizes(y->dim(),rows(),cols(),M_trans,x.dim());

  //
  // U_hat = [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ]
  //         [                  A' * P_XF_hat                   ]
  //         

  const size_type
    n_R	= Q_R_->cols(),
    m 	= A() ? A()->cols()  : 0;

  if( M_trans == BLAS_Cpp::no_trans ) {
    //
    // y =  b*y + a * U_hat * x
    // 
    //   = b*y + a * [  Q_R' * G * P_XF_hat + Q_R' * A_bar * P_plus_hat ] * x
    //               [                  A' * P_XF_hat                   ]
    // 
    //  =>
    // 
    // y1 = b * y1 + a * Q_R' * G * P_XF_hat * x + a * Q_R' * A_bar * P_plus_hat * x
    // y2 = b * y2 + a * A' * P_XF_hat * x
    // 
    DVectorSlice
      y1 = (*y)(1,n_R),
      y2 = m ? (*y)(n_R+1,n_R+m) : DVectorSlice();
    SpVector
      P_XF_hat_x,
      P_plus_hat_x;
    // P_XF_hat_x = P_XF_hat * x
    if( P_XF_hat().nz() )
      V_MtV( &P_XF_hat_x, P_XF_hat(), no_trans, x );
    // P_plus_hat_x = P_plus_hat * x
    if(P_plus_hat().nz())
      V_MtV( &P_plus_hat_x, P_plus_hat(), no_trans, x );
    // y1 = b * y1
    if(b==0.0)      y1=0.0;
    else if(b!=1.0) Vt_S(&y1,b);
    // y1 += a * Q_R' * G * P_XF_hat_x
    if(P_XF_hat().nz())
      Vp_StPtMtV( &y1, a, Q_R(), trans, G(), no_trans, P_XF_hat_x() );
    // y1 += a * Q_R' * A_bar * P_plus_hat_x
    if(P_plus_hat().nz())
      Vp_StPtMtV( &y1, a, Q_R(), trans, A_bar(), no_trans, P_plus_hat_x() );
    if(m) {
      // y2 = b * y2
      if(b==0.0)      y2=0.0;
      else if(b!=1.0) Vt_S(&y2,b);
      // y2 += a * A' * P_XF_hat_x
      if(P_XF_hat().nz())
        Vp_StMtV( &y2, a, *A(), trans, P_XF_hat_x() );
    }
  }
  else if( M_trans == BLAS_Cpp::trans ) {
    //
    // y =  b*y + a * U_hat' * x
    // 
    //   = b*y + a * [  P_XF_hat' * G * Q_R + P_plus_hat' * A_bar' * Q_R, P_XF_hat' * A ] * [ x1 ]
    //                                                                                      [ x2 ]
    //  =>
    // 
    // y = b * y + a * P_XF_hat' * G * Q_R * x1 + a * P_plus_hat' * A_bar' * Q_R * x1
    //     + a * P_XF_hat' * A * x2
    // 
    const SpVectorSlice
      x1 = x(1,n_R),
      x2 = m ? x(n_R+1,n_R+m) : SpVectorSlice(NULL,0,0,0);
    SpVector
      Q_R_x1;
    // Q_R_x1 = Q_R * x1
    V_MtV( &Q_R_x1, Q_R(), no_trans, x1 );
    // y = b*y
    if(b ==0.0)     *y = 0.0;
    else if(b!=1.0) Vt_S( y, b );
    // y += a * P_XF_hat' * G * Q_R_x1
    if(P_XF_hat().nz())
      Vp_StPtMtV( y, a, P_XF_hat(), trans, G(), no_trans, Q_R_x1() );
    // y += a * P_plus_hat' * A_bar' * Q_R_x1
    if(P_plus_hat().nz())
      Vp_StPtMtV( y, a, P_plus_hat(), trans, A_bar(), trans, Q_R_x1() );
    // y += a * P_XF_hat' * A * x2
    if( m && P_XF_hat().nz() )
      Vp_StPtMtV( y, a, P_XF_hat(), trans, *A(), no_trans, x2 );
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT(true);	// Invalid value for M_trans
  }
}

// public member functions for QPSchur::ActiveSet

QPSchur::ActiveSet::ActiveSet(
  const schur_comp_ptr_t& schur_comp
  ,MSADU::PivotTolerances   pivot_tols
  )
  :schur_comp_(schur_comp)
  ,pivot_tols_(pivot_tols)
  ,initialized_(false)
  ,test_(false)
  ,qp_(NULL)
  ,x_init_(NULL)
  ,n_(0)
  ,n_R_(0)
  ,m_(0)
  ,q_plus_hat_(0)
  ,q_F_hat_(0)
  ,q_C_hat_(0)
{}

void QPSchur::ActiveSet::initialize(
  QP& qp, size_type num_act_change, const int ij_act_change[]
  ,const EBounds bnds[], bool test,  bool salvage_init_schur_comp
  ,std::ostream *out, EOutputLevel output_level )
{
  using LinAlgOpPack::V_mV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::V_InvMtV;
  using LinAlgOpPack::Vp_StPtMtV;
  using AbstractLinAlgPack::V_MtV;
  using AbstractLinAlgPack::Mp_StPtMtP;
  using AbstractLinAlgPack::M_StMtInvMtM;
  using DenseLinAlgPack::sym;
  typedef MatrixSymAddDelUpdateable MSADU;
  namespace GPMSTP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  
  const size_type
    n		= qp.n(),
    n_R		= qp.n_R(),
    n_X		= n - n_R,
    m		= qp.m();
  const QP::x_init_t
    &x_init = qp.x_init();
  const QP::l_x_X_map_t
    &l_x_X_map = qp.l_x_X_map();
  const QP::i_x_X_map_t
    &i_x_X_map = qp.i_x_X_map();
  const DVectorSlice
    b_X = qp.b_X();
  const DVectorSlice
    g = qp.g();
  const MatrixSymOp
    &G = qp.G();
  const QP::Constraints
    &constraints = qp.constraints();
  const size_type
    m_breve	= constraints.m_breve();

  try {

  // Count the number of each type of change
  size_type 
    q_plus_hat		= 0,
    q_F_hat			= 0,
    q_C_hat			= 0;
  if( num_act_change ) {
    for( size_type k = 1; k <= num_act_change; ++k ) {
      const int ij = ij_act_change[k-1];
      const EBounds bnd = bnds[k-1];
      if( ij < 0 ) {
        // Initially fixed variable being freed.
        if( x_init(-ij) == FREE ) {
          std::ostringstream omsg;
          omsg
            << "QPSchur::ActiveSet::initialize(...) : Error, "
            << "The variable x(" << -ij << ") is not initially fixed and can not "
            << "be freed by ij_act_change[" << k-1 << "]\n";
          throw std::invalid_argument( omsg.str() );
        }
        if( x_init(-ij) == EQUALITY ) {
          std::ostringstream omsg;
          omsg
            << "QPSchur::ActiveSet::initialize(...) : Error, "
            << "The variable x(" << -ij << ") is equality fixed and therefore can not "
            << "be freed by ij_act_change[" << k-1 << "]\n";
          throw std::invalid_argument( omsg.str() );
        }
        ++q_F_hat;
      }
      else {
        // Constraint being added to the active-set
        if( ij <= n ) {
          // Fixing a variable to a bound
          EBounds x_init_bnd = x_init(ij);
          if( x_init_bnd == FREE ) {
            // initially free variable being fixed
            ++q_plus_hat;
          }
          else if ( x_init_bnd == EQUALITY ) {
            // ToDo: Throw exception
            TEUCHOS_TEST_FOR_EXCEPT(true);
          }
          else if( x_init_bnd == bnd ) {
            // ToDo: Throw exception
            TEUCHOS_TEST_FOR_EXCEPT(true);
          }
          else {
            // Initially fixed variable being fixed to another bound
            ++q_F_hat;	// We must free the variable first
            ++q_C_hat;	// Now we fix it to a different bound.
          }
        }
        else {
          // Adding a general inequality (or equality) constraint
          if( ij > n + m_breve ) {
            // ToDo: Throw exception
            TEUCHOS_TEST_FOR_EXCEPT(true);
          }		
          ++q_plus_hat;
        }
      }
    } 
  }

  const size_type
    q_D_hat = (n - n_R) - q_F_hat,
    q_D_hat_max = n_X;

  // Now let's set stuff up: ij_map, constr_norm, bnds and part of d_hat
  const size_type
    q_hat = q_plus_hat + q_F_hat + q_C_hat,
    q_hat_max = n_X + n,	// If all the initially fixed variables where freed
                // Then all the degrees of freedom where used up with other constraints.
    q_F_hat_max = n_X,
    q_plus_hat_max = n;
    
  ij_map_.resize(q_hat_max);
  constr_norm_.resize(q_hat_max);
  bnds_.resize(q_hat_max);
  d_hat_.resize(q_hat_max);	// set the terms involving the bounds first.

  if( num_act_change ) {
    size_type s = 0;
    for( size_type k = 1; k <= num_act_change; ++k ) {
      const int ij = ij_act_change[k-1];
      const EBounds bnd = bnds[k-1];
      if( ij < 0 ) {
        // Initially fixed variable being freed.
        ij_map_[s]		= ij;
        constr_norm_[s]	= 1.0;
        bnds_[s]		= FREE;
        d_hat_[s]		= - g(-ij);		// - g^X_{l^{(-)}}
        ++s;
      }
      else {
        // Constraint being added to the active-set
        if( ij <= n ) {
          // Fixing a variable to a bound
          EBounds x_init_bnd = x_init(ij);
          if( x_init_bnd == FREE ) {
            // initially free variable being fixed
            ij_map_[s]		= ij;
            constr_norm_[s]	= 1.0;
            bnds_[s]		= bnd;
            d_hat_[s]		= constraints.get_bnd(ij,bnd);
            ++s;
          }
          else {
            // Initially fixed variable being fixed to another bound
            // Free the variable first
            ij_map_[s]		= ij;
            constr_norm_[s]	= 1.0;
            bnds_[s]		= FREE;
            d_hat_[s]		= - g(ij);		// - g^X_{l^{(-)}}
            ++s;
            // Now fix to a different bound
            ij_map_[s]		= ij;
            constr_norm_[s]	= 1.0;
            bnds_[s]		= bnd;
            d_hat_[s]		= constraints.get_bnd(ij,bnd) - b_X(l_x_X_map(ij));
            ++s;
          }
        }
        else {
          // Adding a general inequality (or equality) constraint
          ij_map_[s]		= ij;
          constr_norm_[s]	= 1.0;	// ToDo: We need to compute this in an efficient way!
          bnds_[s]		= bnd;
          d_hat_[s]		= constraints.get_bnd(ij,bnd);	// \bar{b}_{j^{(+)}}
          ++s;
        }
      }
    }
    TEUCHOS_TEST_FOR_EXCEPT( !( s == q_hat ) );
  }

  // Setup P_XF_hat_ and P_plus_hat_
  P_XF_hat_row_.resize(q_F_hat_max);
  P_XF_hat_col_.resize(q_F_hat_max);
  P_FC_hat_row_.resize(q_F_hat_max);
  P_FC_hat_col_.resize(q_F_hat_max);
  P_plus_hat_row_.resize(q_plus_hat_max);
  P_plus_hat_col_.resize(q_plus_hat_max);
  if(q_hat) {
    // See ConstrainedOptPack_QPSchur.hpp for description of P_XF_hat and P_plus_hat
    size_type
      k_XF_hat = 0,	// zero based
      k_plus_hat = 0;	// zero based
    ij_map_t::const_iterator
      ij_itr 		= ij_map_.begin(),
      ij_itr_end	= ij_itr + q_hat;
    for( size_type s = 1; ij_itr != ij_itr_end; ++ij_itr, ++s ) {
      const int ij = *ij_itr;
      EBounds x_init_ij;
      if( ij < 0 ) {
        const size_type i = -ij;
        TEUCHOS_TEST_FOR_EXCEPT( !(  i <= n  ) );
        // [P_XF_hat](:,s) = e(i)
        P_XF_hat_row_[k_XF_hat] = i;
        P_XF_hat_col_[k_XF_hat] = s;
        ++k_XF_hat;
      }
      else if( !(ij <= n && (x_init_ij = x_init(ij)) != FREE ) ) {
        const size_type j = ij;
        TEUCHOS_TEST_FOR_EXCEPT( !(  0 < j && j <= n + m_breve  ) );
        // [P_plus_hat](:,s) = e(j)
        P_plus_hat_row_[k_plus_hat] = j;
        P_plus_hat_col_[k_plus_hat] = s;
        ++k_plus_hat;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPT( !(  k_XF_hat == q_F_hat  ) );
    TEUCHOS_TEST_FOR_EXCEPT( !(  k_plus_hat == q_plus_hat  ) );
  }
  P_XF_hat_.initialize_and_sort(
    n,q_hat,q_F_hat,0,0,GPMSTP::BY_ROW
    ,q_F_hat ? &P_XF_hat_row_[0] : NULL
    ,q_F_hat ? &P_XF_hat_col_[0] : NULL
    ,test
    );
  P_plus_hat_.initialize_and_sort(
    n+m_breve,q_hat,q_plus_hat,0,0,GPMSTP::BY_ROW
    ,q_plus_hat ? &P_plus_hat_row_[0] : NULL
    ,q_plus_hat ? &P_plus_hat_col_[0] : NULL
    ,test
    );
  
  // Setup P_FC_hat_
  if( q_C_hat ) {
    throw std::logic_error(
      error_msg(__FILE__,__LINE__,"QPSchur::ActiveSet::initialize(...) : "
            "Error, q_C_hat != 0, now supported yet!"));  
    // ToDo: We should implement this but it is unlikely to be needed
  }
  P_FC_hat_.initialize_and_sort(
    q_hat,q_hat,q_C_hat,0,0,GPMSTP::BY_ROW
    ,q_C_hat ? &P_FC_hat_row_[0] : NULL
    ,q_C_hat ? &P_FC_hat_col_[0] : NULL
    ,test
    );

  // Setup Q_XD_hat_
  Q_XD_hat_row_.resize(q_D_hat_max);
  Q_XD_hat_col_.resize(q_D_hat_max);
  if(q_D_hat) {
    // See ConstrainedOptPack_QPSchur.hpp for description of Q_XD_hat
    size_type
      k_XD_hat = 0;	// zero based
    GenPermMatrixSlice::const_iterator
      Q_X_itr = qp.Q_X().begin();	// This is sorted by row already!
    P_row_t::const_iterator
      XF_search 		= P_XF_hat_row_.begin(),	// These are already sorted by row!
      XF_search_end 	= XF_search + q_F_hat;
    for( size_type l = 1; l <= n_X; ++l, ++Q_X_itr ) {
      const size_type i = Q_X_itr->row_i();	// Already sorted by row
      // Look for i in XF
      for( ; XF_search != XF_search_end && *XF_search < i; ++XF_search ) ;
      if( XF_search == XF_search_end || (XF_search < XF_search_end && *XF_search > i) ) {
        // We went right past i and did not find it so
        // this variable has not been freed so lets add it!
        Q_XD_hat_row_[k_XD_hat] = i;
        Q_XD_hat_col_[k_XD_hat] = k_XD_hat + 1;
        ++k_XD_hat;
      }
    }
    TEUCHOS_TEST_FOR_EXCEPT( !(  k_XD_hat == q_D_hat  ) );
  }
  Q_XD_hat_.initialize(
      n,q_D_hat,q_D_hat,0,0,GPMSTP::BY_ROW	// Should already be sorted by row!
    , q_D_hat ? &Q_XD_hat_row_[0] : NULL
    , q_D_hat ? &Q_XD_hat_col_[0] : NULL
    ,test
    );

  // Setup l_fxfx
  l_fxfx_.resize(q_D_hat_max);
  if(q_D_hat) {
    for( size_type k = 0; k < q_D_hat; ++k ) {
      l_fxfx_[k] = l_x_X_map(Q_XD_hat_row_[k]);
      TEUCHOS_TEST_FOR_EXCEPT( !(  l_fxfx_[k] != 0  ) );
    }
  }

  // Set the rest of the terms in d_hat involving matrices
  //
  // d_hat += - P_XF_hat' * G * b_XX - P_plus_hat' * A_bar' * b_XX
  // 
  // where: b_XX = Q_X * b_X
  // 
  if( q_hat ) {
    SpVector b_XX;
    V_MtV( &b_XX, qp.Q_X(), BLAS_Cpp::no_trans, b_X );
    Vp_StPtMtV( &d_hat_(1,q_hat), -1.0, P_XF_hat_, BLAS_Cpp::trans
      , G, BLAS_Cpp::no_trans, b_XX() );
    Vp_StPtMtV( &d_hat_(1,q_hat), -1.0, P_plus_hat_, BLAS_Cpp::trans
      , constraints.A_bar(), BLAS_Cpp::trans, b_XX() );
  }

  // Setup U_hat
  U_hat_.initialize( &G, m ? &qp.A() : NULL, &constraints.A_bar()
    , &qp.Q_R(), &P_XF_hat_, &P_plus_hat_ );

  // Set the rest of the members
  test_		= test;
  qp_			= &qp;
  x_init_		= &x_init;
  n_			= n;
  n_R_		= n_R;
  m_			= m;
  m_breve_	= m_breve;
  q_plus_hat_	= q_plus_hat;
  q_F_hat_	= q_F_hat;
  q_C_hat_	= q_C_hat;

  // Resize storage for z_hat, p_z_hat, mu_D_hat, and p_mu_D_hat and set to zero by default
  z_hat_.resize(q_hat_max);
  p_z_hat_.resize(q_hat_max);
  mu_D_hat_.resize(n_X);
  p_mu_D_hat_.resize(n_X);

  initialized_ = true;	// set to true tenatively so that we can
              // print this stuff.

  if( out && (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
    *out
      << "\nPrint definition of Active-Set before the Schur complement is formed...\n";
    dump_act_set_quantities( *this, *out, false );
  }

  // Initialize and factorize the schur complement
  if( q_hat ) {
    // Temporary storage for S (dense)
    DMatrix S_store(q_hat+1,q_hat+1);
    DMatrixSliceSym S_sym( S_store(2,q_hat+1,1,q_hat), BLAS_Cpp::lower );
    MatrixSymPosDefCholFactor S(&S_store());
    // S = -1.0 * U_hat' * inv(Ko) * U_hat
    M_StMtInvMtM( &S, -1.0, U_hat_, BLAS_Cpp::trans, qp.Ko()
      , MatrixSymNonsing::DUMMY_ARG );
    // Now add parts of V_hat
    if( q_F_hat ) {
      // S += P_XF_hat' * G * P_XF_hat
      Mp_StPtMtP( &S, 1.0, MatrixSymOp::DUMMY_ARG, qp_->G(), P_XF_hat_, BLAS_Cpp::no_trans );
    }
    if( q_F_hat && q_plus_hat ) {
      // S += P_XF_hat' * A_bar * P_plus_hat + P_plus_hat' * A_bar' * P_XF_hat
      AbstractLinAlgPack::syr2k(
        qp_->constraints().A_bar()
        ,BLAS_Cpp::no_trans, 1.0
        ,P_XF_hat_, BLAS_Cpp::no_trans
        ,P_plus_hat_, BLAS_Cpp::no_trans
        ,1.0, &S );
    }
    if( q_F_hat && q_C_hat ) {
      // S += P_FC_hat + P_FC_hat'
      throw std::logic_error(
        error_msg(__FILE__,__LINE__,"QPSchur::ActiveSet::initialize(...) : "
              "Error, q_C_hat != 0, now supported yet!"));  
      // ToDo: We should implement this but it is unlikely to be needed
    }

    if( out && (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
      *out
        << "\nIninitial Schur Complement before it is nonsingular:\n"
        << "\nS_hat =\nLower triangular part (ignore nonzeros above diagonal)\n"
        << S_store();
    }
    // Initialize and factorize the schur complement!
    try {
      schur_comp().update_interface().initialize(
        S_sym,q_hat_max,true,MSADU::Inertia( q_plus_hat + q_C_hat, 0, q_F_hat )
        ,pivot_tols() );
    }
    catch(const MSADU::WarnNearSingularUpdateException& excpt) {
      if( out && (int)output_level >= QPSchur::OUTPUT_BASIC_INFO ) {
        *out
          << "\nActiveSet::initialize(...) : " << excpt.what()
          << std::endl;
      }
    }
    catch(const MSADU::SingularUpdateException& excpt) {
      if( out && (int)output_level >= QPSchur::OUTPUT_BASIC_INFO ) {
        *out
        << "\nActiveSet::initialize(...) : " << excpt.what()
        << std::endl;
      }
      if(salvage_init_schur_comp) {
        if( out && (int)output_level >= QPSchur::OUTPUT_BASIC_INFO ) {
          *out
            << "\nsalvage_init_schur_comp == true\n"
            << "We will attempt to add as many rows/cols of the "
            << "initial Schur complement as possible ...\n";
        }
        // ToDo: We will build the schur complement one row/col at a time
        // skipping those updates that cause it to become singular.  For each
        // update that causes the schur compplement to become singular we
        // will remove the corresponding change.

        throw; // For now just rethrow the exception!
      }
      else {
        if( out && (int)output_level >= QPSchur::OUTPUT_BASIC_INFO ) {
          *out
            << "\nsalvage_init_schur_comp == false\n"
            << "We will just throw this singularity exception out of here ...\n";
        }
        throw;
      }
    }
    if( out && (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
      *out
        << "\nSchur Complement after factorization:\n"
        << "\nS_hat =\n"
         << schur_comp().op_interface();
    }
  }
  else {
    schur_comp().update_interface().set_uninitialized();
  }

  // Success, we are initialized!
  initialized_ = true;
  return;

  }	// end try
  catch(...) {
    initialized_ = false;
    throw;
  }
}

void QPSchur::ActiveSet::refactorize_schur_comp()
{
  // ToDo: Finish Me
  TEUCHOS_TEST_FOR_EXCEPT(true);
}

bool QPSchur::ActiveSet::add_constraint(
  size_type ja, EBounds bnd_ja, bool update_steps
  ,std::ostream *out, EOutputLevel output_level
  ,bool force_refactorization
  ,bool allow_any_cond
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using DenseLinAlgPack::dot;
  using AbstractLinAlgPack::dot;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::Vp_StPtMtV;
  using LinAlgOpPack::V_InvMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  
  typedef AbstractLinAlgPack::EtaVector eta_t;

  assert_initialized();

  bool wrote_output = false;

  const QPSchurPack::QP::Constraints
    &constraints = qp_->constraints();

  if( is_init_fixed(ja) && (*x_init_)(ja) == bnd_ja ) {
    //
    // This is a variable that was initially fixed, then freed and now
    // is being fixed back to its original bound.
    //
    // Here we will shrink the augmented KKT system.
    //
    const size_type q_hat = this->q_hat();
    const size_type sd = s_map(-int(ja));
    const size_type la = qp_->l_x_X_map()(ja);
    TEUCHOS_TEST_FOR_EXCEPT( !( sd ) );
    TEUCHOS_TEST_FOR_EXCEPT( !( la ) );
    wrote_output = remove_augmented_element(
      sd,force_refactorization
      ,MatrixSymAddDelUpdateable::EIGEN_VAL_POS
      ,out,output_level,allow_any_cond);
    // We must remove (ja,sd) from P_XF_hat
    P_row_t::iterator
      itr = std::lower_bound( P_XF_hat_row_.begin(), P_XF_hat_row_.begin()+q_F_hat_, ja );
    TEUCHOS_TEST_FOR_EXCEPT( !(  itr != P_XF_hat_row_.end()  ) );
    const size_type p = itr - P_XF_hat_row_.begin();
    std::copy( P_XF_hat_row_.begin() + p + 1, P_XF_hat_row_.begin()+q_F_hat_,
      P_XF_hat_row_.begin() + p );
    std::copy( P_XF_hat_col_.begin() + p + 1, P_XF_hat_col_.begin()+q_F_hat_,
      P_XF_hat_col_.begin() + p );
    // Deincrement all counters in permutation matrices for removed element
    if( q_F_hat_ > 1 )
      deincrement_indices( sd, &P_XF_hat_col_, q_F_hat_-1 );
    if( q_C_hat_ > 0 )
      deincrement_indices( sd, &P_FC_hat_col_, q_C_hat_ );
    if( q_plus_hat_ > 0 )
      deincrement_indices( sd, &P_plus_hat_col_, q_plus_hat_ );
    //
    // Add the multiplier for mu_D_hat(...)
    //
    const size_type q_D_hat = this->q_D_hat();
    // add l_fxfx(q_D_hat+1) = la
    l_fxfx_[q_D_hat] = la;
    // add mu_D_hat(q_D_hat+1) = 0.0
    mu_D_hat_[q_D_hat] = 0.0;
    // add p_mu_D_hat(q_D_hat+1) = 1.0
    if(update_steps)
      p_mu_D_hat_[q_D_hat] = 1.0;
    // Insert the pair (ja,q_D_hat+1) into Q_XD_hat(...) (sorted by row)
    insert_pair_sorted(ja,q_D_hat+1,q_D_hat+1,&Q_XD_hat_row_,&Q_XD_hat_col_);
    //
    // Update the counts
    //
    --q_F_hat_;
  }
  else {
    //
    // Expand the size of the schur complement to add the new constraint
    //
     
    // Compute the terms for the update
    
    value_type			d_p = 0.0;
    const size_type		q_hat = this->q_hat();
    Workspace<value_type> t_hat_ws(wss,q_hat);
    DVectorSlice t_hat(t_hat_ws.size()?&t_hat_ws[0]:NULL,t_hat_ws.size());
    value_type			alpha_hat = 0.0;
    bool				changed_bounds = false;
    size_type           sd = 0; // Only used if changed_bounds == true
        
    if( ja <= n_ && !is_init_fixed(ja) ) {
      //
      // Fix an initially free variable is being fixed
      // 
      // u_p = [ Q_R' * e(ja) ] <: R^(n_R+m)
      //       [      0       ]
      //
      const size_type
        la = qp_->Q_R().lookup_col_j(ja);
      TEUCHOS_TEST_FOR_EXCEPT( !(  la  ) );
      const eta_t u_p = eta_t(la,n_R_+m_);
      // r = inv(Ko)*u_p
      DVector r;	// ToDo: Make this sparse!
      V_InvMtV( &r, qp_->Ko(), no_trans, u_p() );
      // t_hat = - U_hat' * r
      if(q_hat)
        V_StMtV( &t_hat(), -1.0, U_hat_, trans, r() );
      // alpha_hat = - u_p ' * r
      alpha_hat = - dot( u_p(), r() );
      // d_p = \bar{b}_{j^{(+)}}
      d_p = constraints.get_bnd( ja, bnd_ja );

      changed_bounds = false;
    }
    else if ( is_init_fixed(ja) ) {
      //
      // An intially fixed variable was freed and
      // is now being fixed to the other bound.
      //
      // Here we must expand the augmented KKT system for this
      // simple change.
      //
      // u_p = 0
      //
      // v_p = e(sd) <: R^(q_hat), where sd = s_map(-ja)
      //
      sd = s_map(-int(ja));
      TEUCHOS_TEST_FOR_EXCEPT( !( sd ) );
      const size_type la = qp_->l_x_X_map()(ja);
      TEUCHOS_TEST_FOR_EXCEPT( !( la ) );
      // t_hat = e(sd)
      t_hat = 0.0;
      t_hat(sd) = 1.0;
      // alpha_hat = 0.0
      alpha_hat = 0.0;
      // d_p = \bar{b}_{j^{(+)}} - b_X(la)
      d_p = constraints.get_bnd( ja, bnd_ja ) - qp_->b_X()(la);

      changed_bounds = true;
    }
    else {	// ja > n
      //
      // Add an extra equality or inequality constraint.
      //
      // u_p = [ Q_R' * A_bar * e(ja) ] n_R
      //       [        0             ] m
      const eta_t e_ja = eta_t(ja,n_+m_breve_);
      const MatrixOp &A_bar = constraints.A_bar();
      DVector u_p( n_R_ + m_ );	// ToDo: make this sparse
      Vp_StPtMtV( &u_p(1,n_R_), 1.0, qp_->Q_R(), trans, A_bar, no_trans, e_ja(), 0.0 );
      if( m_ )
        u_p(n_R_+1,n_R_+m_) = 0.0;
      // r = inv(Ko) * u_p
      DVector r;	// ToDo: Make this sparse!
      V_InvMtV( &r, qp_->Ko(), no_trans, u_p() );
      if(q_hat) {
        // t_hat = v_p - U_hat' * r
        //    where: v_p = P_XF_hat' * A_bar * e(ja)
        V_StMtV( &t_hat(), -1.0, U_hat_, trans, r() );
        Vp_StPtMtV( &t_hat(), 1.0, P_XF_hat_, trans, A_bar, no_trans, e_ja() );
      }
      // alpha_hat = - u_p ' * r
      alpha_hat = - dot( u_p(), r() );
      // d_p = \bar{b}_{j^{(+)}} - b_X' * Q_X' * A_bar * e(ja)
      // 
      // d_p = \bar{b}_{j^{(+)}}
      d_p = constraints.get_bnd( ja, bnd_ja );
      if(n_ > n_R_) {
        // d_p += - b_X' * Q_X' * A_bar * e(ja)
        r.resize( n_ - n_R_ );	// reuse storage
        Vp_StPtMtV( &r(), 1.0, qp_->Q_X(), trans, A_bar, no_trans, e_ja(), 0.0 );			
        d_p += - dot( qp_->b_X(), r() );
      }
      
      changed_bounds = false;
    }
    
    // Update the schur complement if nonsingular.  These
    // with throw exceptions if the matrix is singular
    // or has the wrong inertia
    if(q_hat) {
      try {
        schur_comp().update_interface().augment_update(
          &t_hat(), alpha_hat, force_refactorization
          ,MatrixSymAddDelUpdateable::EIGEN_VAL_NEG
          ,MSADU::PivotTolerances(
            pivot_tols().warning_tol
            ,allow_any_cond ? 0.0 : pivot_tols().singular_tol
            ,allow_any_cond ? 0.0 : pivot_tols().wrong_inertia_tol ) );
      }
      catch(const MSADU::WarnNearSingularUpdateException& excpt) {
        if( out && (int)output_level >= QPSchur::OUTPUT_BASIC_INFO ) {
          *out
            << "\nActiveSet::add_constraint(...) : " << excpt.what()
            << std::endl;
          wrote_output = true;
        }
      }
    }
    else {
      schur_comp().update_interface().initialize(
        alpha_hat, (n_-n_R_) + n_-m_ );
    }
    // Update the rest of the augmented KKT system
    if( changed_bounds )
      ++q_C_hat_;
    else
      ++q_plus_hat_;
    const size_type q_hat_new = q_F_hat_ + q_C_hat_ + q_plus_hat_;
    // Add ij_map(q_hat) = ja to ij_map(...)
    ij_map_[q_hat_new - 1]	= ja;
    // Add constr_norm(q_hat) to constr_norm(...)
    constr_norm_[q_hat_new - 1]	= 1.0;	// ToDo: Compute this for real!
    // Add bnds(q_hat)_ = bnd_ja to bnds(...)
    bnds_[q_hat_new - 1] 	= bnd_ja;
    // Augment d_hat = [ d_hat; d_p ]
    d_hat_(q_hat_new) = d_p;
    // Augment z_hat with new (zeroed) multiplier value, z_hat = [ z_hat; 0 ]
    z_hat_(q_hat_new) = 0.0;
    if( update_steps ) {
      // Set the step for this multiplier to 1.0, p_z_hat = [ p_z_hat; 1 ]
      p_z_hat_(q_hat_new) = 1.0;
    }
    if( !changed_bounds ) {
      // Insert (ja, q_hat_new) into P_plus_hat, sorted by row
      insert_pair_sorted(ja,q_hat_new,q_plus_hat_,&P_plus_hat_row_,&P_plus_hat_col_);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT( !( sd ) );
      // Insert (sd,q_hat_new) into P_FC_hat, sorted by row)
      insert_pair_sorted(sd,q_hat_new,q_C_hat_,&P_FC_hat_row_,&P_FC_hat_col_);
    }
  }
  // Update the permutation matrices and U_hat
  reinitialize_matrices(test_);
  return wrote_output;
}

bool QPSchur::ActiveSet::drop_constraint(
  int jd, std::ostream *out, EOutputLevel output_level
  ,bool force_refactorization, bool allow_any_cond
  )
{
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using DenseLinAlgPack::dot;
  using AbstractLinAlgPack::dot;
  using LinAlgOpPack::V_StMtV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_MtV;
  using LinAlgOpPack::Vp_StPtMtV;
  using LinAlgOpPack::V_InvMtV;
  using AbstractLinAlgPack::transVtMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  
  typedef AbstractLinAlgPack::EtaVector eta_t;

  assert_initialized();

  bool wrote_output = false;

  if( jd < 0 ) {
    //
    // A variable initially fixed is being freed.
    // Increase the dimension of the augmented the KKT system!
    //
    size_type
      q_hat      = this->q_hat(),
      q_F_hat    = this->q_F_hat(),
      q_plus_hat = this->q_plus_hat(),
      q_D_hat    = this->q_D_hat();
    // Get indexes
    const size_type id = -jd;
    TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= id && id <= n_  ) );
    const size_type ld = qp_->l_x_X_map()(-jd);
    TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= ld && ld <= n_ - n_R_  ) );
    size_type kd; // Find kd (this is unsorted)
    {for( kd = 1; kd <= q_D_hat; ++kd ) {
      if( l_fxfx_[kd-1] == ld ) break;
    }}
    TEUCHOS_TEST_FOR_EXCEPT( !(  kd <= q_D_hat  ) );
    // Get references
    const MatrixSymOp
      &G           = qp_->G();
    const DVectorSlice
      g            = qp_->g();
    const MatrixOp
      &A_bar       = qp_->constraints().A_bar();
    const MatrixSymOpNonsing
      &Ko          = qp_->Ko();
    const MatrixOp
      &U_hat       = this->U_hat();
    const GenPermMatrixSlice
      &Q_R         = qp_->Q_R(),
      &Q_X         = qp_->Q_X(),
      &P_XF_hat    = this->P_XF_hat(),
      &P_plus_hat  = this->P_plus_hat();
    const DVectorSlice
      b_X          = qp_->b_X();
    //
    // Compute the update quantities to augmented KKT system
    //
    // e_id
    eta_t e_id(id,n_);
    // u_p = [ Q_R'*G*e_id ; A'*e_id ] <: R^(n_R+m)
    DVector u_p(n_R_+m_);
    Vp_StPtMtV( &u_p(1,n_R_), 1.0, Q_R, trans, G, no_trans, e_id(), 0.0 );
    if( m_ )
      V_MtV( &u_p(n_R_+1,n_R_+m_), qp_->A(), trans, e_id() );
    const value_type
      nrm_u_p = DenseLinAlgPack::norm_inf( u_p() );
    // sigma = e_id'*G*e_id <: R
    const value_type
      sigma = transVtMtV( e_id(), G, no_trans, e_id() );
    // d_p = - g(id) - b_X'*(Q_X'*G*e_id) <: R
    DVector Q_X_G_e_id(Q_X.cols());
    Vp_StPtMtV( &Q_X_G_e_id(), 1.0, Q_X, trans, G, no_trans, e_id(), 0.0 );
    const value_type
      d_p = -g(id) - dot( b_X, Q_X_G_e_id() );
    // r = inv(Ko)*u_p <: R^(n_R+m)
    DVector r;
    if( nrm_u_p > 0.0 )
      V_InvMtV( &r, Ko, no_trans, u_p() );
    // t_hat = v_p - U_hat'*r
    // where: v_p = P_XF_hat'*G*e_id + P_plus_hat'*A_bar'*e_id <: R^(q_hat)
    Workspace<value_type>
      t_hat_ws(wss,q_hat);
    DVectorSlice
      t_hat(&t_hat_ws[0],q_hat);
    if(q_hat) {
      t_hat = 0.0;
      // t_hat += v_p
      if(q_F_hat_)
        Vp_StPtMtV( &t_hat(), 1.0, P_XF_hat, trans, G, no_trans, e_id() );
      if(q_plus_hat_)
        Vp_StPtMtV( &t_hat(), 1.0, P_plus_hat, trans, A_bar, trans, e_id() );
      // t_hat += U_hat'*r
      if( nrm_u_p > 0.0 )
        Vp_MtV( &t_hat(), U_hat, trans, r() );
    }
    // alpha_hat = sigma - u_p'*r
    const value_type
      alpha_hat = sigma - ( nrm_u_p > 0.0 ? dot(u_p(),r()) : 0.0 );
    //
    // Update the schur complement if nonsingular.  These
    // with throw exceptions if the matrix is singular
    // or has the wrong inertia
    //
    if(q_hat) {
      try {
        schur_comp().update_interface().augment_update(
          &t_hat(), alpha_hat, force_refactorization
          ,MatrixSymAddDelUpdateable::EIGEN_VAL_POS );
      }
      catch(const MSADU::WarnNearSingularUpdateException& excpt) {
        if( out && (int)output_level >= QPSchur::OUTPUT_BASIC_INFO ) {
          *out
            << "\nActiveSet::drop_constraint(...) : " << excpt.what()
            << std::endl;
          wrote_output = true;
        }
      }
    }
    else {
      schur_comp().update_interface().initialize(
        alpha_hat, (n_-n_R_) + n_-m_ );
    }
    //
    // Remove multiplier from mu_D_hat(...)
    //
    // remove l_fxfx(kd) == ld from l_fxfx(...)
    std::copy( l_fxfx_.begin() + kd, l_fxfx_.begin() + q_D_hat
      , l_fxfx_.begin() + (kd-1) );
    // remove mu_D_hat(kd) from mu_D_hat(...)
    std::copy( mu_D_hat_.begin() + kd, mu_D_hat_.begin() + q_D_hat
      , mu_D_hat_.begin() + (kd-1) );
    // remove Q_XD_hat(id,ld) from Q_XD_hat(...)
    P_row_t::iterator
      itr = std::lower_bound( Q_XD_hat_row_.begin(), Q_XD_hat_row_.begin()+q_D_hat, id );
    TEUCHOS_TEST_FOR_EXCEPT( !(  itr != Q_XD_hat_row_.end()  ) );
    const size_type p = itr - Q_XD_hat_row_.begin();
    std::copy( Q_XD_hat_row_.begin() + p + 1, Q_XD_hat_row_.begin()+q_D_hat,
      Q_XD_hat_row_.begin() + p );
    std::copy( Q_XD_hat_col_.begin() + p + 1, Q_XD_hat_col_.begin()+q_D_hat,
      Q_XD_hat_col_.begin() + p );
    if( q_D_hat > 1 )
      deincrement_indices( kd, &Q_XD_hat_col_, q_D_hat-1 );
    //
    // Update the counts
    //
    ++q_F_hat_;
    q_hat = this->q_hat();
    //
    // Add the elements for newly freed variable
    //
    // add ij_map(q_hat) == -id to ij_map(...)
    ij_map_[q_hat-1] = -id;
    // add s_map(-id) == q_hat to s_map(...)
    // ToDo: implement s_map(...)
    // add bnds(q_hat) == FREE to bnds(...)
    bnds_[q_hat-1] = FREE;
    // add d_hat(q_hat) == d_p to d_hat(...)
    d_hat_[q_hat-1] = d_p;
    // add p_X(ld) == 0 to the end of z_hat(...)
    z_hat_[q_hat-1] = 0.0; // This is needed so that (z_hat + beta*t_D*p_z_hat)(q_hat) == 0
    // Insert (id,q_hat) into P_XF_hat sorted by row
    insert_pair_sorted(id,q_hat,q_F_hat_,&P_XF_hat_row_,&P_XF_hat_col_);
  }
  else {
    //
    // Shrink the dimension of the augmented KKT system to remove the constraint!
    //
    const size_type q_hat = this->q_hat();
    const size_type sd = s_map(jd);
    TEUCHOS_TEST_FOR_EXCEPT( !( sd ) );
    wrote_output = remove_augmented_element(
      sd,force_refactorization
      ,MatrixSymAddDelUpdateable::EIGEN_VAL_NEG
      ,out,output_level,allow_any_cond
      );
    if( is_init_fixed(jd) ) {
      // This must be an intially fixed variable, currently fixed at a different bound.
      // We must remove this element from P_FC_hat(...)
      const size_type sd1 = s_map(-jd); // This is the position in the schur complement where first freed
      TEUCHOS_TEST_FOR_EXCEPT( !( sd1 ) );
      // Remove P_FC_hat(sd1,sd) from P_FC_hat(...)
      P_row_t::iterator
        itr = std::lower_bound( P_FC_hat_row_.begin(), P_FC_hat_row_.begin()+q_C_hat_, sd1 );
      TEUCHOS_TEST_FOR_EXCEPT( !(  itr != P_FC_hat_row_.end()  ) );
      const size_type p = itr - P_FC_hat_row_.begin();
      std::copy( P_FC_hat_row_.begin() + p + 1, P_FC_hat_row_.begin()+q_C_hat_,
        P_FC_hat_row_.begin() + p );
      std::copy( P_FC_hat_col_.begin() + p + 1, P_FC_hat_col_.begin()+q_C_hat_,
        P_FC_hat_col_.begin() + p );
      --q_C_hat_;
    }
    else {
      // We must remove P_plus_hat(jd,sd) from P_plus_hat(...)
      P_row_t::iterator
        itr = std::lower_bound( P_plus_hat_row_.begin(), P_plus_hat_row_.begin()+q_plus_hat_, jd );
      TEUCHOS_TEST_FOR_EXCEPT( !(  itr != P_plus_hat_row_.end()  ) );
      const size_type p = itr - P_plus_hat_row_.begin();
      std::copy( P_plus_hat_row_.begin() + p + 1, P_plus_hat_row_.begin()+q_plus_hat_,
        P_plus_hat_row_.begin() + p );
      std::copy( P_plus_hat_col_.begin() + p + 1, P_plus_hat_col_.begin()+q_plus_hat_,
        P_plus_hat_col_.begin() + p );
      --q_plus_hat_;
    }
    // Deincrement all counters in permutation matrices for removed element
    if( q_F_hat_ > 0 )
      deincrement_indices( sd, &P_XF_hat_col_, q_F_hat_ );
    if( q_C_hat_ > 0 )
      deincrement_indices( sd, &P_FC_hat_col_, q_C_hat_ );
    if( q_plus_hat_ > 0 )
      deincrement_indices( sd, &P_plus_hat_col_, q_plus_hat_ );
  }
  // Update the permutation matrices and U_hat
  reinitialize_matrices(test_);
  return wrote_output;
}

bool QPSchur::ActiveSet::drop_add_constraints(
  int jd, size_type ja, EBounds bnd_ja, bool update_steps
  ,std::ostream *out, EOutputLevel output_level
  )
{
  bool wrote_output = false;
  if( drop_constraint( jd, out, output_level, false, true ) )
    wrote_output = true;
  if( add_constraint( ja, bnd_ja, update_steps, out, output_level, true, true ) )
    wrote_output = true;
  return wrote_output;
}

QPSchur::ActiveSet::QP&
QPSchur::ActiveSet::qp()
{
  assert_initialized();
  return *qp_;
}

const QPSchur::ActiveSet::QP&
QPSchur::ActiveSet::qp() const
{
  assert_initialized();
  return *qp_;
}

size_type QPSchur::ActiveSet::q_hat() const
{
  assert_initialized();
  return q_plus_hat_ + q_F_hat_ + q_C_hat_;
}

size_type QPSchur::ActiveSet::q_plus_hat() const
{
  assert_initialized();
  return q_plus_hat_;
}

size_type QPSchur::ActiveSet::q_F_hat() const
{
  assert_initialized();
  return q_F_hat_;
}

size_type QPSchur::ActiveSet::q_C_hat() const
{
  assert_initialized();
  return q_C_hat_;
}

size_type QPSchur::ActiveSet::q_D_hat() const
{
  assert_initialized();
  return (n_ - n_R_) - q_F_hat_;  // n^{X} - \hat{q}^{F}
}

int QPSchur::ActiveSet::ij_map( size_type s ) const
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= s && s <= this->q_hat()  ) );
  return ij_map_[s-1];
}

size_type QPSchur::ActiveSet::s_map( int ij ) const
{
  ij_map_t::const_iterator
    begin	= ij_map_.begin(),
    end		= begin + q_hat(),
    itr = std::find( begin, end, ij );
  return ( itr != end ? (itr - begin) + 1 : 0 );
}

value_type QPSchur::ActiveSet::constr_norm( size_type s ) const
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= s && s <= this->q_hat()  ) );
  return constr_norm_(s);
}

EBounds QPSchur::ActiveSet::bnd( size_type s ) const
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= s && s <= this->q_hat()  ) );
  return bnds_[s-1];
}

size_type QPSchur::ActiveSet::l_fxfx( size_type k ) const
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  1 <= k && k <= this->q_D_hat()  ) );
  return l_fxfx_[k-1];
}

const QPSchur::U_hat_t& QPSchur::ActiveSet::U_hat() const
{
  assert_initialized();
  return U_hat_;
}

const MatrixSymOpNonsing& QPSchur::ActiveSet::S_hat() const
{
  assert_initialized();
  return schur_comp().op_interface();
}

const GenPermMatrixSlice& QPSchur::ActiveSet::P_XF_hat() const
{
  assert_initialized();
  return P_XF_hat_;
}

const GenPermMatrixSlice& QPSchur::ActiveSet::P_FC_hat() const
{
  assert_initialized();
  return P_FC_hat_;
}

const GenPermMatrixSlice& QPSchur::ActiveSet::P_plus_hat() const
{
  assert_initialized();
  return P_plus_hat_;
}

const GenPermMatrixSlice& QPSchur::ActiveSet::Q_XD_hat() const
{
  assert_initialized();
  return Q_XD_hat_;
}

const DVectorSlice QPSchur::ActiveSet::d_hat() const
{
  assert_initialized();
  return d_hat_(1,q_hat());
}

DVectorSlice QPSchur::ActiveSet::z_hat()
{
  assert_initialized();
  return z_hat_(1,q_hat());
}

const DVectorSlice QPSchur::ActiveSet::z_hat() const
{
  assert_initialized();
  return z_hat_(1,q_hat());
}

DVectorSlice QPSchur::ActiveSet::p_z_hat()
{
  assert_initialized();
  return p_z_hat_(1,q_hat());
}

const DVectorSlice QPSchur::ActiveSet::p_z_hat() const
{
  assert_initialized();
  return p_z_hat_(1,q_hat());
}

DVectorSlice QPSchur::ActiveSet::mu_D_hat()
{
  assert_initialized();
  return mu_D_hat_(1,q_D_hat());
}

const DVectorSlice QPSchur::ActiveSet::mu_D_hat() const
{
  assert_initialized();
  return mu_D_hat_(1,q_D_hat());
}

DVectorSlice QPSchur::ActiveSet::p_mu_D_hat()
{
  assert_initialized();
  return p_mu_D_hat_(1,q_D_hat());
}

const DVectorSlice QPSchur::ActiveSet::p_mu_D_hat() const
{
  assert_initialized();
  return p_mu_D_hat_(1,q_D_hat());
}

bool QPSchur::ActiveSet::is_init_fixed( size_type j ) const
{
  assert_initialized();
  return j <= n_ && (*x_init_)(j) != FREE;
}

bool QPSchur::ActiveSet::all_dof_used_up() const
{
  return n_ - m_ == (n_ - n_R_) - q_F_hat_ + q_C_hat_ + q_plus_hat_;
}

// private member functions for QPSchur::ActiveSet

void QPSchur::ActiveSet::assert_initialized() const
{
  if( !initialized_ )
    throw std::logic_error(
      error_msg(__FILE__,__LINE__,"QPSchur::ActiveSet::assert_initialized() : Error, "
      "The active set has not been initialized yet!") );
}

void QPSchur::ActiveSet::assert_s( size_type s) const
{
  TEUCHOS_TEST_FOR_EXCEPT( !(  s <= q_hat()  ) );	// ToDo: Throw an exception
}

void QPSchur::ActiveSet::reinitialize_matrices(bool test)
{
  namespace GPMSTP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;

  const size_type q_hat = this->q_hat();
  const size_type q_D_hat = this->q_D_hat();

  P_XF_hat_.initialize(
    n_,q_hat,q_F_hat_,0,0,GPMSTP::BY_ROW
    ,q_F_hat_ ? &P_XF_hat_row_[0] : NULL
    ,q_F_hat_ ? &P_XF_hat_col_[0] : NULL
    ,test
    );
  P_FC_hat_.initialize(
    q_hat,q_hat,q_C_hat_,0,0,GPMSTP::BY_ROW
    ,q_C_hat_ ? &P_FC_hat_row_[0] : NULL
    ,q_C_hat_ ? &P_FC_hat_col_[0] : NULL
    ,test
    );
  P_plus_hat_.initialize(
    n_+m_breve_,q_hat,q_plus_hat_,0,0,GPMSTP::BY_ROW
    ,q_plus_hat_ ? &P_plus_hat_row_[0] : NULL
    ,q_plus_hat_ ? &P_plus_hat_col_[0] : NULL
    ,test
    );
  Q_XD_hat_.initialize(
    n_,q_D_hat,q_D_hat,0,0,GPMSTP::BY_ROW
    ,q_D_hat ? &Q_XD_hat_row_[0] : NULL
    ,q_D_hat ? &Q_XD_hat_col_[0] : NULL
    ,test
    );
  U_hat_.initialize( 
    &qp_->G(), m_ ? &qp_->A() : NULL, &qp_->constraints().A_bar()
    ,&qp_->Q_R(), &P_XF_hat_, &P_plus_hat_);
}

bool QPSchur::ActiveSet::remove_augmented_element(
  size_type sd, bool force_refactorization
  ,MatrixSymAddDelUpdateable::EEigenValType eigen_val_drop
  ,std::ostream *out, EOutputLevel output_level
  ,bool allow_any_cond
  )
{
  bool wrote_output = false;
  const size_type q_hat = this->q_hat();
  // Delete the sd row and column for S_hat
  try {
    schur_comp().update_interface().delete_update(
      sd,force_refactorization,eigen_val_drop
      ,MSADU::PivotTolerances(
        pivot_tols().warning_tol
        ,allow_any_cond ? 0.0 : pivot_tols().singular_tol
        ,allow_any_cond ? 0.0 : pivot_tols().wrong_inertia_tol ));
  }
  catch(const MSADU::WarnNearSingularUpdateException& excpt) {
    if( out && (int)output_level >= QPSchur::OUTPUT_BASIC_INFO ) {
      *out
        << "\nActiveSet::drop_constraint(...) : " << excpt.what()
        << std::endl;
      wrote_output = true;
    }
  }
  // Remove the ij_map(s) = jd element from ij_map(...)
  std::copy( ij_map_.begin() + sd, ij_map_.begin() + q_hat
         , ij_map_.begin() + (sd-1) );
  // Remove the constr_norm(s) elment from constr_norm(...)
  std::copy( constr_norm_.begin() + sd, constr_norm_.begin() + q_hat
         , constr_norm_.begin() + (sd-1) );
  // Remove the bnds(s) element from bnds(...)
  std::copy( bnds_.begin() + sd, bnds_.begin() + q_hat
         , bnds_.begin() + (sd-1) );
  // Remove the d_hat(s) element from d_hat(...)
  std::copy( d_hat_.begin() + sd, d_hat_.begin() + q_hat
         , d_hat_.begin() + (sd-1) );
  // Remove the z_hat(s) element from z_hat(...)
  std::copy( z_hat_.begin() + sd, z_hat_.begin() + q_hat
         , z_hat_.begin() + (sd-1) );
  // Remove the p_z_hat(s) element from p_z_hat(...)
  std::copy( p_z_hat_.begin() + sd, p_z_hat_.begin() + q_hat
         , p_z_hat_.begin() + (sd-1) );
  return wrote_output;
}

// public member functions for QPSchur

value_type QPSchur::DEGENERATE_MULT = std::numeric_limits<value_type>::min();

void QPSchur::pivot_tols( MSADU::PivotTolerances pivot_tols )
{
  act_set_.pivot_tols(pivot_tols);
}

QPSchur::MSADU::PivotTolerances QPSchur::pivot_tols() const
{
  return act_set_.pivot_tols();
}

QPSchur::QPSchur(
  const schur_comp_ptr_t&   schur_comp
  ,size_type                max_iter
  ,value_type               max_real_runtime
  ,value_type               feas_tol
  ,value_type               loose_feas_tol
  ,value_type               dual_infeas_tol
  ,value_type               huge_primal_step
  ,value_type               huge_dual_step
  ,value_type               warning_tol
  ,value_type               error_tol
  ,size_type                iter_refine_min_iter
  ,size_type                iter_refine_max_iter
  ,value_type               iter_refine_opt_tol
  ,value_type               iter_refine_feas_tol
  ,bool                     iter_refine_at_solution
  ,bool                     salvage_init_schur_comp
  ,MSADU::PivotTolerances   pivot_tols
  )
  :schur_comp_(schur_comp)
  ,max_iter_(max_iter)
  ,max_real_runtime_(max_real_runtime)
  ,feas_tol_(feas_tol)
  ,loose_feas_tol_(loose_feas_tol)
  ,dual_infeas_tol_(dual_infeas_tol)
  ,huge_primal_step_(huge_primal_step)
  ,huge_dual_step_(huge_dual_step)
  ,warning_tol_(warning_tol)
  ,error_tol_(error_tol)
  ,iter_refine_min_iter_(iter_refine_min_iter)
  ,iter_refine_max_iter_(iter_refine_max_iter)
  ,iter_refine_opt_tol_(iter_refine_opt_tol)
  ,iter_refine_feas_tol_(iter_refine_feas_tol)
  ,iter_refine_at_solution_(iter_refine_at_solution)
  ,salvage_init_schur_comp_(salvage_init_schur_comp)
  ,act_set_(schur_comp,pivot_tols)
{}

QPSchur::ESolveReturn QPSchur::solve_qp(
  QP& qp
  ,size_type num_act_change, const int ij_act_change[], const EBounds bnds[]
  ,std::ostream *out, EOutputLevel output_level, ERunTests test_what
  ,DVectorSlice* x, SpVector* mu, DVectorSlice* lambda, SpVector* lambda_breve
  ,size_type* iter, size_type* num_adds, size_type* num_drops
  )
{
  using std::setw;
  using std::endl;
  using std::right;
  using DenseLinAlgPack::norm_inf;
  using AbstractLinAlgPack::norm_inf;
  using LinAlgOpPack::V_InvMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();
  using StopWatchPack::stopwatch;

  const value_type inf = std::numeric_limits<value_type>::max();

  if( !out )
    output_level = NO_OUTPUT;

  const int dbl_min_w = 20;
  const int dbl_w = (out ? my_max(dbl_min_w,int(out->precision()+8)) : 20 );

  // Set the schur complement
  act_set_.set_schur_comp( schur_comp_ );

  ESolveReturn
    solve_return = SUBOPTIMAL_POINT;

  // Print QPSchur output header
  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\n*** Entering QPSchur::solve_qp(...)\n";
  }

  // Start the timer!
  stopwatch timer;
  timer.reset();
  timer.start();

  // Print the definition of the QP to be solved.
  if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
    *out
      << "\n*** Dump the definition of the QP to be solved ***\n";
    qp.dump_qp(*out);
  }
  
  // Print warm start info.
  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\n*** Warm start info\n"
      << "\nNumber of variables                                          = " << qp.n()
      << "\nNumber of initially fixed variables (not in Ko)              = " << (qp.n() - qp.n_R())
      << "\nNumber of changes to the initial KKT system (num_act_change) = " << num_act_change << endl;
    const size_type
      n = qp.n();
    size_type
      num_var_fixed = 0, num_var_freed = 0, num_gen_equ = 0, num_gen_inequ = 0;
    for( size_type s = 1; s <= num_act_change; ++s ) {
      const int ij = ij_act_change[s-1];
      const EBounds bnd = bnds[s-1];
      if( ij < 0 )
        ++num_var_freed;
      else if( ij < n )
        ++num_var_fixed;
      else if( bnd == EQUALITY )
        ++num_gen_equ;
      else
        ++num_gen_inequ;
    }
    *out
      << "\n    Number of initially fixed variables freed from a bound   = " << num_var_freed
      << "\n    Number of initially free variables fixed to a bound      = " << num_var_fixed
      << "\n    Number of general equality constraints added             = " << num_gen_equ
      << "\n    Number of general inequality constraints added           = " << num_gen_inequ << endl;
  }
  
  if( num_act_change > 0 && (int)output_level >= (int)OUTPUT_ACT_SET ) {
    *out
      << endl
      << right << setw(5) << "s"
      << right << setw(20) << "ij_act_change"
      << right << setw(10) << "bnds" << endl
      << right << setw(5) << "---"
      << right << setw(20) << "-------------"
      << right << setw(10) << "--------" << endl;
    for( size_type s = 1; s <= num_act_change; ++s )
      *out
        << right << setw(5) << s
        << right << setw(20) << ij_act_change[s-1]
        << right << setw(10) << bnd_str(bnds[s-1]) << endl;
  }

  // Initialize the active set.
  try {
    act_set_.initialize( qp, num_act_change, ij_act_change, bnds
      , test_what == RUN_TESTS, salvage_init_schur_comp(), out, output_level );
    // If this throws a WrongInteriaUpdateExecption it will be
    // thrown clean out of here!
  }
  catch( const MatrixSymAddDelUpdateable::SingularUpdateException& excpt ) {
    if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
      *out
        << "\n*** Error in initializing schur complement\n"
        << excpt.what() << std::endl
        << "\nSetting num_act_change = 0 and proceeding with a cold start...\n";
    }
    act_set_.initialize( qp, num_act_change = 0, ij_act_change, bnds
      , test_what == RUN_TESTS, salvage_init_schur_comp(), out, output_level );
  }

  // Compute vo =  inv(Ko) * fo
  Workspace<value_type> vo_ws(wss,qp.n_R()+qp.m());
  DVectorSlice vo(&vo_ws[0],vo_ws.size());
  V_InvMtV( &vo, qp.Ko(), BLAS_Cpp::no_trans, qp.fo() );

  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\nSolution to the initial KKT system, vo = inv(Ko)*fo:\n\n||vo||inf = " << norm_inf(vo) << std::endl;
  }
  if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
    *out
      << "\nvo =\n" << vo;
  }

  // ////////////////////////////////////////////////
  // Remove constraints until we are dual feasible.
  // 
  // Here we are essentially performing a primal algorithm where we are only
  // removing constraints.  If the dual variables are not dual feasible then
  // we will remove the one with the largest scaled dual feasibility
  // violation then compute the dual variables again.  Eventually we
  // will come to a point where we have a dual feasible point.  If
  // we have to, we will remove all of the inequality constraints and then
  // this will by default be a dual feasible point (i.e. we picked all the
  // wrong inequality constraints).
  // 
  // The difficulty here is in dealing with near degenerate constraints.
  // If a constraint is near degenerate then we would like to not drop
  // it since we may have to add it again later.
  // 
  *iter = 0;
  *num_adds = 0;
  *num_drops = 0;
  // Print header for removing constraints
  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\n***"
      << "\n*** Removing constriants until we are dual feasible"
      << "\n***\n"
      << "\n*** Start by removing constraints within the Schur complement first\n";
  }
  // Print summary header for max_viol and jd.
  if( (int)OUTPUT_ITER_SUMMARY <= (int)output_level 
    && (int)output_level < (int)OUTPUT_ITER_QUANTITIES
    && num_act_change > 0  )
  {
    *out     
      << "\nIf max_viol > 0 and jd != 0 then constraint jd will be dropped from the active set\n\n"
      << right << setw(dbl_w)	<< "max_viol"
      << right << setw(5)	<< "sd"
      << right << setw(5)	<< "jd"		<< endl
      << right << setw(dbl_w)	<< "--------------"
      << right << setw(5)	<< "----"
      << right << setw(5)	<< "----"	<< endl;
  }
  for( int k = num_act_change; k > 0; --k, ++(*iter) ) {
    // Check runtime
    if( timeout_return(&timer,out,output_level) )
      return MAX_RUNTIME_EXEEDED_FAIL;
    // Compute z_hat (z_hat = inv(S_hat)*(d_hat - U_hat'*vo))
    DVectorSlice z_hat = act_set_.z_hat();
    calc_z( act_set_.S_hat(), act_set_.d_hat(), act_set_.U_hat(), &vo
      , &z_hat );
    // Determine if we are dual feasible.
    value_type	max_viol = 0.0;	// max scaled violation of dual feasability.
    size_type	jd = 0;			// indice of constraint with max scaled violation.
    DVectorSlice::iterator
      z_itr = z_hat.begin();
    const size_type q_hat = act_set_.q_hat();	// Size of schur complement system.
    // Print header for s, z_hat(s), bnd(s), viol, max_viol and jd
    if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
      *out
        << "\nLooking for a constraint with the maximum dual infeasiblity to drop...\n\n"
        << right << setw(5)		<< "s"
        << right << setw(dbl_w)	<< "z_hat"
        << right << setw(20)	<< "bnd"
        << right << setw(dbl_w)	<< "viol"
        << right << setw(dbl_w)	<< "max_viol"
        << right << setw(5)		<< "jd"	<< endl
        << right << setw(5)		<< "----"
        << right << setw(dbl_w)	<< "--------------"
        << right << setw(20)	<< "--------------"
        << right << setw(dbl_w)	<< "--------------"
        << right << setw(dbl_w)	<< "--------------"
        << right << setw(5)		<< "----"	<< endl;
    }
    for( int s = 1; s <= q_hat; ++s, ++z_itr ) {
      int j = act_set_.ij_map(s);
      if( j > 0 ) {
        // This is for an active constraint not initially in Ko so z_hat(s) = mu(j)
        EBounds bnd = act_set_.bnd(s);
        value_type viol;  // Get the amount of violation and fix near degeneracies
        const int dual_feas_status
          = correct_dual_infeas(
            j,bnd,0.0,act_set_.constr_norm(s),dual_infeas_tol(),DEGENERATE_MULT
            ,out,output_level,false
            ,"z_hat(s)", &(*z_itr), &viol );
        if( dual_feas_status < 0  && viol < max_viol ) {
          max_viol = viol;
          jd = j;
        }
        // Print row for s, z_hat(s), bnd(s), viol, max_viol and jd
        if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
          *out
            << right << setw(5)		<< s
            << right << setw(dbl_w)	<< *z_itr
            << right << setw(20)	<< bnd_str(bnd)
            << right << setw(dbl_w)	<< viol
            << right << setw(dbl_w)	<< max_viol
            << right << setw(5)		<< jd	<< endl;
        }
      }
    }
    // Print row of max_viol and jd
    if( (int)OUTPUT_ITER_SUMMARY <= (int)output_level 
      && (int)output_level < (int)OUTPUT_ITER_QUANTITIES )
    {
      *out
        << right << setw(dbl_w)	<< max_viol
        << right << setw(5)	    << act_set_.s_map(jd)
        << right << setw(5)	    << jd					<< endl;
    }
    if( jd == 0 ) break;	// We have a dual feasible point w.r.t. these constraints
    // Remove the active constraint with the largest scaled violation.
    act_set_.drop_constraint( jd, out, output_level, true, true );
    ++(*iter);
    ++(*num_drops);
    // Print U_hat, S_hat and d_hat.
    if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
      *out
        << "\nPrinting active set after dropping constraint jd = " << jd << " ...\n";
      dump_act_set_quantities( act_set_, *out );
    }
  }

  // Print how many constraints where removed from the schur complement
  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\nThere where " << (*num_drops)
      << " constraints dropped from the schur complement from the initial guess of the active set.\n";
  }

  // Compute v
  Workspace<value_type> v_ws(wss,qp.n_R()+qp.m());
  DVectorSlice v(&v_ws[0],v_ws.size());
  if( act_set_.q_hat() > 0 ) {
    calc_v( qp.Ko(), &qp.fo(), act_set_.U_hat(), act_set_.z_hat(), &v );
    if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
      *out
        << "\nSolution to the system; v = inv(Ko)*(fo - U_hat*z_hat):\n"
        << "\n||v||inf = " << norm_inf(v) << std::endl;
    }
    if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
      *out
        << "\nv =\n" << v;
    }
  }
  else {
    v = vo;
  }

  // Set x
  set_x( act_set_, v, x );
  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\nCurrent guess for unknowns x:\n\n||x||inf = " << norm_inf(*x) << std::endl;
  }
  if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
    *out
      << "\nx =\n" << *x;
  }

  //
  // Determine if any initially fixed variables need to be freed by checking mu_D_hat.
  //
  if( act_set_.q_D_hat() ) {
    if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
      *out << "\n*** Second, free initially fixed variables not in Ko\n";
    }
    const QPSchurPack::QP::i_x_X_map_t&  i_x_X_map = act_set_.qp().i_x_X_map();
    const QPSchurPack::QP::x_init_t&     x_init    = act_set_.qp().x_init();
    // Print summary header for max_viol and jd.
    if( (int)OUTPUT_ITER_SUMMARY <= (int)output_level 
      && (int)output_level < (int)OUTPUT_ITER_QUANTITIES )
    {
      *out     
        << "\nIf max_viol > 0 and id != 0 then the variable x(id) will be freed from its initial bound\n\n"
        << right << setw(dbl_w)	<< "max_viol"
        << right << setw(5)	<< "kd"
        << right << setw(5)	<< "id"		<< endl
        << right << setw(dbl_w)	<< "--------------"
        << right << setw(5)	<< "----"
        << right << setw(5)	<< "----"	<< endl;
    }
    size_type q_D_hat = act_set_.q_D_hat(); // This will be deincremented
    while( q_D_hat > 0 ) {
      // Check runtime
      if( timeout_return(&timer,out,output_level) )
        return MAX_RUNTIME_EXEEDED_FAIL;
      // mu_D_hat = ???
      DVectorSlice mu_D_hat = act_set_.mu_D_hat();
      calc_mu_D( act_set_, *x, v, &mu_D_hat );
      // Determine if we are dual feasible.
      value_type	max_viol = 0.0;	// max scaled violation of dual feasability.
      int			id = 0;			// indice of variable with max scaled violation.
      size_type   kd = 0;
      DVectorSlice::iterator
        mu_D_itr = mu_D_hat.begin();
      // Print header for k, mu_D_hat(k), bnd, viol, max_viol and id
      if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
        *out
          << "\nLooking for a variable bound with the max dual infeasibility to drop...\n\n"
          << right << setw(5)		<< "k"
          << right << setw(dbl_w)	<< "mu_D_hat"
          << right << setw(20)	<< "bnd"
          << right << setw(dbl_w)	<< "viol"
          << right << setw(dbl_w)	<< "max_viol"
          << right << setw(5)		<< "id"	<< endl
          << right << setw(5)		<< "----"
          << right << setw(dbl_w)	<< "--------------"
          << right << setw(20)	<< "--------------"
          << right << setw(dbl_w)	<< "--------------"
          << right << setw(dbl_w)	<< "--------------"
          << right << setw(5)		<< "----"	<< endl;
      }
      for( int k = 1; k <= q_D_hat; ++k, ++mu_D_itr ) {
        int
          i = i_x_X_map(act_set_.l_fxfx(k));
        EBounds
          bnd = x_init(i);
        value_type viol;  // Get the amount of violation and fix near degeneracies
        const int dual_feas_status
          = correct_dual_infeas(
            i,bnd,0.0,1.0,dual_infeas_tol(),DEGENERATE_MULT
            ,out,output_level,false
            ,"mu_D_hat(k)", &(*mu_D_itr), &viol );
        if( dual_feas_status < 0  && viol < max_viol ) {
          max_viol = viol;
          kd = k;
          id = i;
        }
        // Print row for k, mu_D_hat(k), bnd, viol, max_viol and jd
        if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
          *out
            << right << setw(5)		<< k
            << right << setw(dbl_w)	<< *mu_D_itr
            << right << setw(20)	<< bnd_str(bnd)
            << right << setw(dbl_w)	<< viol
            << right << setw(dbl_w)	<< max_viol
            << right << setw(5)		<< id	<< endl;
        }
      }
      // Print row of max_viol and id
      if( (int)OUTPUT_ITER_SUMMARY <= (int)output_level 
        && (int)output_level < (int)OUTPUT_ITER_QUANTITIES )
      {
        *out
          << right << setw(dbl_w)	<< max_viol
          << right << setw(5)	<< kd
          << right << setw(5)	<< id         << endl;
      }
      if( id == 0 ) break;	// We have a dual feasible point w.r.t. these variable bounds
      // Remove the active constraint with the largest scaled violation.
      act_set_.drop_constraint( -id, out, output_level, true, true );
      ++(*iter);
      ++(*num_adds);
      --q_D_hat;
      // Print U_hat, S_hat and d_hat.
      if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
        *out
          << "\nPrinting active set after freeing initially fixed variable id = " << id << " ...\n";
        dump_act_set_quantities( act_set_, *out );
      }
      if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
        *out
          << "\nSolution to the new KKT system; z_hat = inv(S_hat)*(d_hat - U_hat'*vo), v = inv(Ko)*(fo - U_hat*z_hat):\n";
      }
      // Compute z_hat (z_hat = inv(S_hat)*(d_hat - U_hat'*vo))
      calc_z( act_set_.S_hat(), act_set_.d_hat(), act_set_.U_hat(), &vo
          , &act_set_.z_hat() );
      if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
        *out
          << "\n||z_hat||inf = " << norm_inf(act_set_.z_hat()) << std::endl;
      }
      if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
        *out
          << "\nz_hat =\n" << act_set_.z_hat();
      }
      // Compute v
      calc_v( qp.Ko(), &qp.fo(), act_set_.U_hat(), act_set_.z_hat(), &v );
      if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
        *out
          << "\n||v||inf = " << norm_inf(v) << std::endl;
      }
      if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
        *out
          << "\nv =\n" << v;
      }
      // Set x
      set_x( act_set_, v, x );
      if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
        *out
          << "\nCurrent guess for unknowns x:\n\n||x||inf = " << norm_inf(*x) << std::endl;
      }
      if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
        *out
          << "\nx =\n" << *x;
      }
    }
  }

  // Print how many initially fixed variables where freed
  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\nThere where " << (*num_adds)
      << " initially fixed variables not in Ko that were freed and added to the schur complement.\n";
  }

  // Run the primal dual algorithm
  size_type  iter_refine_num_resid = 0, iter_refine_num_solves = 0;
  solve_return = qp_algo(
    PICK_VIOLATED_CONSTRAINT
    ,out, output_level, test_what
    ,vo, &act_set_, &v
    ,x, iter, num_adds, num_drops
    ,&iter_refine_num_resid, &iter_refine_num_solves
    ,&timer
    );

  if( solve_return != OPTIMAL_SOLUTION )
    set_x( act_set_, v, x );

  // Correct the sign of near degenerate multipliers in case it has not been done yet!
  if( solve_return != SUBOPTIMAL_POINT && act_set_.q_hat() ) {
    const size_type q_hat = act_set_.q_hat();
    DVectorSlice z_hat = act_set_.z_hat();
    for( size_type s = 1; s <= q_hat; ++s ) {
      const int       j    = act_set_.ij_map(s);
      value_type      viol = 0.0;
      const EBounds   bnd  = act_set_.bnd(s);
      if(bnd == FREE)
        continue;
      const int dual_feas_status
        = correct_dual_infeas(
          j,bnd,0.0,1.0,dual_infeas_tol(),DEGENERATE_MULT
          ,out,output_level,true,"z_hat(s)",&z_hat(s),&viol );
      if( dual_feas_status < 0 ) {
        solve_return = SUBOPTIMAL_POINT;
        break;
      }
    }
  }
  if( solve_return != SUBOPTIMAL_POINT && act_set_.q_D_hat() ) {
    const GenPermMatrixSlice&          Q_XD_hat = act_set_.Q_XD_hat();
    DVectorSlice                        mu_D_hat = act_set_.mu_D_hat();
    const QPSchurPack::QP::x_init_t&   x_init   = qp.x_init();
    for( GenPermMatrixSlice::const_iterator itr = Q_XD_hat.begin(); itr != Q_XD_hat.end(); ++itr ) {
      const int       i    = itr->row_i();
      value_type      viol = 0.0;
      const EBounds   bnd  = x_init(i);
      TEUCHOS_TEST_FOR_EXCEPT( !(  bnd != FREE  ) );
      const int dual_feas_status
        = correct_dual_infeas(
          i,bnd,0.0,1.0,dual_infeas_tol(),DEGENERATE_MULT
          ,out,output_level,true,"mu_D_hat(k)",&mu_D_hat(itr->col_j()),&viol );
      if( dual_feas_status < 0 ) {
        solve_return = SUBOPTIMAL_POINT;
        break;
      }
    }
  }

  set_multipliers( act_set_, v, mu, lambda, lambda_breve );

  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    switch(solve_return) {
      case OPTIMAL_SOLUTION:
        *out
          << "\n*** Solution found!\n";
        break;
      case MAX_ITER_EXCEEDED:
        *out
          << "\n*** Maximum iterations exceeded!\n";
        break;
      case MAX_RUNTIME_EXEEDED_DUAL_FEAS:
        *out
          << "\n*** Maximum runtime exceeded!\n";
        break;
      case MAX_ALLOWED_STORAGE_EXCEEDED:
        *out
          << "\n*** The maxinum size of the schur complement has been exceeded!\n";
        break;
      case INFEASIBLE_CONSTRAINTS:
        *out
          << "\n*** The constraints are infeasible!\n";
        break;
      case NONCONVEX_QP:
        *out
          << "\n*** The QP appears to be nonconvex but will return current point anyway!\n";
        break;
      case DUAL_INFEASIBILITY:
        *out
          << "\n*** The dual variables are infeasible (numerical roundoff?)!\n";
        break;
      case SUBOPTIMAL_POINT:
        *out
          << "\n*** The current point is suboptimal but we will return it anyway!\n";
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    *out	<< "\nNumber of QP iteratons                                = " << *iter;
    *out	<< "\nNumber of iterative refinement residual calculations  = " << iter_refine_num_resid;
    *out	<< "\nNumber of iterative refinement solves                 = " << iter_refine_num_solves << endl;
    *out	<< "\n||x||inf                = "	<< norm_inf(*x);
    *out	<< "\nmu.nz()                 = "	<< mu->nz();
    *out	<< "\nmax(|mu(i)|)            = " 	<< norm_inf((*mu)());
    *out	<< "\nmin(|mu(i)|)            = " 	<< min_abs((*mu)());
    if(lambda)
      *out
        << "\nmax(|lambda(i)|)        = "	<< norm_inf(*lambda)
        << "\nmin(|lambda(i)|)        = "	<< min_abs(*lambda);
    if(lambda_breve)
      *out
        << "\nlambda_breve.nz()       = "	<< lambda_breve->nz()
        << "\nmax(|lambda_breve(i)|)  = "	<< norm_inf((*lambda_breve)())
        << "\nmin(|lambda_breve(i)|)  = "	<< min_abs((*lambda_breve)());
    *out << std::endl;
  }
  // Print Solution x, lambda and mu
  if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
    *out	<< "\nx =\n" 				<< (*x)();
    *out	<< "\nmu =\n" 				<< (*mu)();
    if(lambda)
      *out << "\nlambda =\n"			<< (*lambda)();
    if(lambda_breve)
      *out << "\nlambda_breve =\n"	<< (*lambda_breve)();
  }
  // Print 'goodby' header.
  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\n*** Leaving QPSchur::solve_qp(...)\n";
  }

  return solve_return;
}

const QPSchur::ActiveSet& QPSchur::act_set() const
{
  return act_set_;
}

// protected member functions for QPSchur

QPSchur::ESolveReturn QPSchur::qp_algo(
  EPDSteps next_step
  , std::ostream *out, EOutputLevel output_level, ERunTests test_what
  , const DVectorSlice& vo, ActiveSet* act_set, DVectorSlice* v
  , DVectorSlice* x, size_type* iter, size_type* num_adds, size_type* num_drops
  , size_type* iter_refine_num_resid, size_type* iter_refine_num_solves
  , StopWatchPack::stopwatch* timer
  )
{
  using std::setw;
  using std::endl;
  using std::right;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using DenseLinAlgPack::dot;
  using DenseLinAlgPack::norm_inf;
  using DenseLinAlgPack::Vt_S;
  using DenseLinAlgPack::V_mV;
  using DenseLinAlgPack::Vp_StV;
  using DenseLinAlgPack::V_VmV;
  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::V_StV;
  using LinAlgOpPack::V_StMtV;
  using AbstractLinAlgPack::dot;
  using AbstractLinAlgPack::norm_inf;
  using AbstractLinAlgPack::V_MtV;
  using AbstractLinAlgPack::EtaVector;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::V_InvMtV;
  using LinAlgOpPack::Vp_StMtV;
  using LinAlgOpPack::Vp_StPtMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  // Print header for "Starting Primal-Dual Iterations"
  if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
    *out
      << "\n*** Starting Primal-Dual Iterations ***\n";
  }

  const int dbl_min_w = 20;
  const int dbl_w = (out ? my_max(dbl_min_w,int(out->precision()+8)): 20 );

  try {

  QPSchurPack::QP
    &qp = act_set->qp();
  const size_type
    n		= qp.n(),
    n_R		= qp.n_R(),
    m		= qp.m(),
    m_breve	= qp.constraints().m_breve();

  Workspace<value_type>
    v_plus_ws(wss,v->dim()),
    z_hat_plus_ws(wss,(n-m)+(n-n_R)),
    p_v_ws(wss,v->dim());
  DVectorSlice
    v_plus(&v_plus_ws[0],v_plus_ws.size()),
    z_hat_plus,
    p_v(&p_v_ws[0],p_v_ws.size());

  const value_type
    inf = std::numeric_limits<value_type>::max();
  size_type itr;	// move to somewhere else?

  // Put these here because they need to be remembered between iterations if a linearly
  // dependent constriant is dropped.
  size_type				ja = 0;		// + indice of violated constraint to add to active set
  size_type				last_ja = 0; // + last_ja to be added
  value_type				con_ja_val;	// value of violated constraint.
  value_type				b_a; // value of the violated bound
  value_type				norm_2_constr;	// norm of violated constraint
  EBounds					bnd_ja;	// bound of constraint ja which is violated.
  bool					can_ignore_ja;	// true if we can ignore a constraint if it is LD.
  bool					assume_lin_dep_ja;
  value_type				gamma_plus;	// used to store the new multipler value for the added
                    // constraint.
  const int				summary_lines_counter_max = 15;
  int						summary_lines_counter = 0;
  long int				jd = 0;	// + indice of constraint to delete from active set.
                  // - indice of intially fixed variable to be freed
  long int				last_jd = 0; // Last jd change to the active set
  value_type				t_P;	// Primal step length (constraint ja made active)
  value_type				t_D;	// Dual step length ( longest step without violating dual
                  // feasibility of currently active constraints ).
  value_type				beta;	// +1 if multiplier of constraint being added is positive
                  // -1 if multiplier of constraint being added is negative.
  bool					warned_degeneracy = false; // Will warn the user if degeneracy
                  // is detected.
  value_type				dual_infeas_scale = 1.0;	// Scaling for determining if a
                  // Lagrange multiplier is near degenerate.
  bool					return_to_init_fixed = false;	// True if the constraint being added
                  // to the active set is a variable returning to its orginally
                                  // fixed variable bound.
  bool                    using_iter_refinement = false; // Will be set to true if instability detected

  for( itr = 0; itr <= max_iter_; ++itr, ++(*iter) ) {
    if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
      *out
        << "\n************************************************"
        << "\n*** qp_iter = " << itr
        << "\n*** q_hat   = " << act_set->q_hat() << std::endl;
    }
    bool schur_comp_update_failed = false;
    switch( next_step ) {	// no break; statements in this switch statement.
      case PICK_VIOLATED_CONSTRAINT: {
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\n*** PICK_VIOLATED_CONSTRAINT\n";
        }
        // Check runtime
        if( timeout_return(timer,out,output_level) )
          return MAX_RUNTIME_EXEEDED_DUAL_FEAS;
        // Save the indice of the last constriant to be added!
        last_ja = ja;
        // Set parts of x that are not currently fixed and may have changed.
        // Also, we want set specifially set those variables that where
        // initially free and then latter fixed to their bounds so that
        // they will not be seen as violated.
        set_x( *act_set, *v, x );
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\n||x||inf = " << norm_inf(*x) << std::endl;
        }
        if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
          *out
            << "\nx =\n" << *x;
        }
        if( test_what == RUN_TESTS ) {
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\nChecking current iteration quantities ...\n";
          }
          const char
            sep_line[] = "\n--------------------------------------------------------------------------------\n";
//					AbstractLinAlgPack::TestingPack::CompareDenseVectors comp_v;

          //
          // Check the optimality conditions of the augmented KKT system!
          //

          // ToDo: Implement
          
          //
          // Check mu_D_hat
          //
          if( act_set->q_D_hat() ) {
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << sep_line
                << "\nComputing: mu_D_hat_calc = -Q_XD_hat'*g - Q_XD_hat'*G*x\n"
                << "    - Q_XD_hat'*A*v(n_R+1:n_R+m) - Q_XD_hat'*A_bar*P_plus_hat*z_hat ...\n";
            }
            Workspace<value_type> mu_D_hat_calc_ws( wss, act_set->q_D_hat() );
            DVectorSlice mu_D_hat_calc( &mu_D_hat_calc_ws[0], mu_D_hat_calc_ws.size() );
            calc_mu_D( *act_set, *x, *v, &mu_D_hat_calc );
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << "\n||mu_D_hat_calc||inf = " << norm_inf(mu_D_hat_calc) << std::endl;
            }
            if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
              *out
                << "\nmu_D_hat_calc =\n" << mu_D_hat_calc;
            }
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << "\nChecking mu_D_hat_calc == mu_D_hat\n";
            }
            DVector mu_D_hat_diff(mu_D_hat_calc.dim());
            LinAlgOpPack::V_VmV( &mu_D_hat_diff(), mu_D_hat_calc(), act_set->mu_D_hat() );
            const value_type
              mu_D_hat_err = norm_inf(mu_D_hat_diff()) / (1.0 + norm_inf(mu_D_hat_calc()));
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << "\n||mu_D_hat_calc-mu_D_hat||inf/(1.0+||mu_D_hat_calc||inf) = "
                << mu_D_hat_err << std::endl;
            }
            TEUCHOS_TEST_FOR_EXCEPTION(
              mu_D_hat_err >= error_tol(), TestFailed
              ,"QPSchur::qp_algo(...) : Error, "
              "||mu_D_hat_calc-mu_D_hat||inf/(1.0+||mu_D_hat_calc||inf) = "
              << mu_D_hat_err << " >= error_tol = " << error_tol()
              );
            if( mu_D_hat_err >= warning_tol() && (int)output_level >= (int)OUTPUT_ACT_SET ) {
              *out
                << "\nWarning! ||mu_D_hat_calc-mu_D_hat||inf/(1.0+||mu_D_hat_calc||inf) = "
                << mu_D_hat_err << " >= warning_tol = " << warning_tol() << std::endl;
            }
          }
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << sep_line;
          }
        }
        act_set->qp().constraints().pick_violated( *x, &ja, &con_ja_val
          , &b_a, &norm_2_constr, &bnd_ja, &can_ignore_ja );
        assume_lin_dep_ja = false;	// Assume this initially.
        if( ja > 0 && act_set->is_init_fixed(ja) && qp.x_init()(ja) == bnd_ja )
          return_to_init_fixed = true;
        else
          return_to_init_fixed = false;
        // Print ja, bnd_ja, can_ignore_ja
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\nja = " << ja	<< endl;
          if(ja) {
            *out
              << "\ncon_ja_val           = "		<< con_ja_val
              << "\nb_a                  = "		<< b_a
              << "\nnorm_2_constr        = "		<< norm_2_constr
              << "\nbnd_ja               = "		<< bnd_str(bnd_ja)
              << "\ncan_ignore_ja        = "		<< bool_str(can_ignore_ja)
              << "\nreturn_to_init_fixed = "		<< bool_str(return_to_init_fixed)
              << endl
              ;
          }
        }
        // Print header for itr, nact, change (ADD, DROP), indice (ja, jb)
        //, bound (LOWER, UPPER, EQUALITY), violation (primal or dual), rank (LD,LI)
        if( (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
          if( summary_lines_counter <= 0 ) {
            summary_lines_counter = summary_lines_counter_max;
            *out
              << endl
              << right << setw(6)	<< "itr"
              << right << setw(6)	<< "qhat"
              << right << setw(6)	<< "q(+)"
              << right << setw(6)	<< "q_F"
              << right << setw(6)	<< "q_C"
              << right << setw(6)	<< "q_D"
              << right << setw(8)	<< "change"
              << right << setw(9)	<< "type"
              << right << setw(6)	<< "j"
              << right << setw(10)	<< "bnd"
              << right << setw(dbl_w)	<< "viol, p_z(jd)"
              << right << setw(6)	<< "rank" << endl
              << right << setw(6)	<< "----"
              << right << setw(6)	<< "----"
              << right << setw(6)	<< "----"
              << right << setw(6)	<< "----"
              << right << setw(6)	<< "----"
              << right << setw(6)	<< "----"
              << right << setw(8)	<< "------"
              << right << setw(9)	<< "-------"
              << right << setw(6)	<< "----"
              << right << setw(10)	<< "--------"
              << right << setw(dbl_w) << "--------------"
              << right << setw(6)	<< "----"	<< endl;
          }
        }
        // Print first part of row for itr, q_hat, q(+), q_D, q_C, q_F, change, type, index, bound, violation
        if( (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
          *out
            << right << setw(6)	<< itr                                   // itr
            << right << setw(6)	<< act_set->q_hat()                      // q_hat
            << right << setw(6)	<< act_set->q_plus_hat()                 // q(+)
            << right << setw(6)	<< act_set->q_F_hat()                    // q_F
            << right << setw(6)	<< act_set->q_C_hat()                    // q_C
            << right << setw(6)	<< act_set->q_D_hat()                    // q_D
            << right << setw(8)  << ( ja ? "ADD" : "-" )                 // change
            << right << setw(9);                                         // type
          if( ja == 0 ) {
            *out << "-";
          }
          else if( act_set->is_init_fixed(ja) ) {
            if( bnd_ja == qp.x_init()(ja) )
              *out << "X_F_X";
            else
              *out << "X_F_C";
          }
          else if( ja <= n ) {
            *out << "R_X";
          }
          else {
            *out << "GEN";
          }
          *out
            << right << setw(6)	<< ja                                    // index
            << right << setw(10) << ( ja ? bnd_str(bnd_ja) : "-" )       // bound
            << right << setw(dbl_w);                                     // violation
          if(ja)
            *out << (con_ja_val - b_a); 
          else
            *out << "-";
          if(!ja)
            *out << right << setw(6) << "-" << endl;                      // rank for last iteration
        }
        bool found_solution = false;
        const size_type sa = act_set->s_map(ja);
        value_type scaled_viol = 0.0;
        if( ja == 0 ) {
          found_solution = true;
        }
        else if( sa != 0 || ( act_set->is_init_fixed(ja) && act_set->s_map(-ja) == 0 ) ) {
          const bool is_most_violated
            = (qp.constraints().pick_violated_policy() == QPSchurPack::Constraints::MOST_VIOLATED);
          if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
            *out
              << "\n\nWarning, we have picked the constriant a(" << ja << ") with violation\n"
              << "(a(ja)'*x - b_a) = (" << con_ja_val << " - " << b_a << ") = " << (con_ja_val - b_a)
              << "\nto add to the active set but it is already part of the active set.\n";
          }
          const EBounds act_bnd = ( sa != 0 ? act_set->bnd(sa) : qp.x_init()(ja) );
          if( act_bnd != bnd_ja ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "However, this is not the same bound as the active bound!\n";
            }
            const value_type
              act_b_a = qp.constraints().get_bnd(ja,act_bnd);
            if( act_bnd == LOWER && act_b_a > b_a || act_bnd == UPPER && act_b_a < b_a ) {
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\nError, c_L_bar(" << ja <<") = " << (act_b_a > b_a ? act_b_a : b_a)
                  << " > c_U_bar(" << ja << ") = " << (act_b_a < b_a ? act_b_a : b_a) << ".\n";
              }
            }
            else {
              TEUCHOS_TEST_FOR_EXCEPT(true); // Should not happen!
            }
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "The constraints are infeasible!  Terminating the QP algorithm!\n";
            }
            return INFEASIBLE_CONSTRAINTS;
          }
          else {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "This is the active bound so this is an indication of instability\n"
                << "in the calculations.\n";
            }
          }
          summary_lines_counter = 0;
          if( !is_most_violated ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\nThis is not the most violated constraint so set\n"
                << "pick_violated_policy = MOST_VIOLATED ...\n";
            }
            summary_lines_counter = 0;
            qp.constraints().pick_violated_policy(QP::Constraints::MOST_VIOLATED);
            next_step = PICK_VIOLATED_CONSTRAINT;
            continue;	// Go back and pick the most violated constraint
          }
          else if( !using_iter_refinement ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\nThis is the most violated constraint but we are currently not using\n"
                << "iterative refinement so let's switch it on ...\n";
            }
            summary_lines_counter = 0;
            EIterRefineReturn status = iter_refine(
              *act_set, out, output_level, -1.0, &qp.fo(), -1.0, act_set->q_hat() ? &act_set->d_hat() : NULL
              ,v, act_set->q_hat() ? &act_set->z_hat() : NULL
              ,iter_refine_num_resid, iter_refine_num_solves
              );
            if( status == ITER_REFINE_IMPROVED || status == ITER_REFINE_CONVERGED ) {
              using_iter_refinement = true; // Use iterative refinement from now on!
              next_step = PICK_VIOLATED_CONSTRAINT;
              continue; // Now go back
            }
            else {
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\nError, iterative refinement was not allowed or failed to improve the residuals.\n"
                  << "Terminating the QP algorithm ...\n";
              }
              return SUBOPTIMAL_POINT;
            }
          }
          else {
            scaled_viol = std::fabs(con_ja_val - b_a)/(1.0 + std::fabs(con_ja_val));
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\n\nThis is the most violated constraint, we are using iterative refinement and\n"
                << "|a("<<ja<<")'*x - b_a| / (1 + |a("<<ja<<")'*x|) = |"<<con_ja_val<<" - "<<b_a
                << "| / (1 + |"<<con_ja_val<<"|) = "<<scaled_viol<< (scaled_viol<loose_feas_tol()?" < ":" > ")
                << "loose_feas_tol = "<<loose_feas_tol();
            }
             if( scaled_viol < loose_feas_tol() ) {
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\nTerminating the algorithm with a near optimal solution!\n";
              }
              found_solution = true;
            }
            else { 
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\nError!  The QP algorithm is terminated!\n";
              }
              return SUBOPTIMAL_POINT;
            }
          }
        }
        else if( act_set->all_dof_used_up()
             && (scaled_viol = std::fabs(con_ja_val - b_a)/(1.0 + std::fabs(con_ja_val))) < feas_tol() )
        {
          if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
            *out
              << "\nWith all the dof used up, the violated inequality constriant "
              << "a("<< ja << ")'*x statisfies:\n"
              << "|a("<<ja<<")'*x - b_a| / (1 + |a("<<ja<<")'*x|) = |" << con_ja_val << " - " << b_a
              << "| / (1 + |"<<con_ja_val<<"|) = "<<scaled_viol<<" < feas_tol = "<<feas_tol();
          }
          if( act_set->qp().constraints().pick_violated_policy() == QP::Constraints::MOST_VIOLATED ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\nThis is the most violated constraint so we will consider this\n"
                << "a near degenerate constraint so we are done!\n";
            }
            found_solution = true;
          }
          else {
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << "\nThis is not the most violated constraint so set pick_violated_policy = "
                << "MOST_VIOLATED and pick another violated constraint ....\n";
            }
            act_set->qp().constraints().pick_violated_policy(QP::Constraints::MOST_VIOLATED);
            next_step = PICK_VIOLATED_CONSTRAINT;
            continue;	// Go back and pick the most violated constraint
          }
        }
        if(found_solution) {
          //
          // Solution found!  All of the inequality constraints are satisfied!
          //
          // Use iterative refinement if we need to
          if( iter_refine_at_solution() && !using_iter_refinement ) {
            // The user has requested iterative refinement at the solution
            // and we have not been using iterative refinement up to this point.
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\nWe think we have found the solution and are not currently using iterative refinement\n"
                << "and iter_refine_at_solution==true so perform iterative refinement ...\n";
            }
            using_iter_refinement = true;
            EIterRefineReturn status = iter_refine(
              *act_set, out, output_level, -1.0, &qp.fo(), -1.0, act_set->q_hat() ? &act_set->d_hat() : NULL
              ,v, act_set->q_hat() ? &act_set->z_hat() : NULL
              ,iter_refine_num_resid, iter_refine_num_solves
              );
            switch(status) {
            case ITER_REFINE_ONE_STEP:
            case ITER_REFINE_IMPROVED:
            case ITER_REFINE_CONVERGED:
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\nIterative refinement may have altered the unknowns so go back and look for another violated constraint ...\n";
              }
              summary_lines_counter = 0;
              next_step = PICK_VIOLATED_CONSTRAINT;
              continue;
            case ITER_REFINE_NOT_PERFORMED:
            case ITER_REFINE_NOT_NEEDED:
            case ITER_REFINE_NOT_IMPROVED:
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\nIterative refinement did not alter the unknowns so exit with this solution...\n";
              }
              set_x( *act_set, *v, x );
              break;
            default:
              TEUCHOS_TEST_FOR_EXCEPT(true); // Local programming error only!
            }
          }
          if( iter_refine_at_solution() || using_iter_refinement ) {
            // The user has requested iterative refinement at the solution
            // or we have been using iterative refinement so we should
            // compute the lagrange multipliers for the initially fixed
            // variables that are still fixed.
            if( act_set->q_D_hat() ) {
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\nRecomputing final mu_D_hat at the solution ...\n";
              }
              calc_mu_D( *act_set, *x, *v, &act_set->mu_D_hat() );
              if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
                *out
                  << "\n||mu_D_hat||inf = " << norm_inf(act_set->mu_D_hat()) << std::endl;
              }
              if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
                *out
                  << "\nmu_D_hat =\n" << act_set->mu_D_hat();
              }
            }
          }
          return OPTIMAL_SOLUTION;	// current point is optimal.
        }
      }
      case UPDATE_ACTIVE_SET: {
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\n*** UPDATE_ACTIVE_SET\n";
        }
        ++(*num_adds);
        if( act_set->all_dof_used_up() || act_set->is_init_fixed(ja) ) {
          // All of the degrees of freedom are currently used up so we must
          // assume that we must remove one of these currently active
          // constraints and replace it with the violated constraint.
          // In this case we know that this will make the schur
          // complement singular so let's just skip the update
          // and set this now.  We also may be here if we are fixing
          // an initially fixed variable to some bound.
          assume_lin_dep_ja = true;
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            if(act_set->all_dof_used_up()) {
              *out
                << "\nAll of the degrees of freedom are used up so "
                << "the constraint ja must be linearly dependent.\n"
                << "The augmented KKT system is definitely singular!\n";
            }
            else {
              *out
                << "\nThis is an initially fixed variable that was freed and "
                << "now is being fixed again.\n"
                << "The augmented KKT system could be singular or nonsingular, "
                << "we don't know at this point.\n";
            }
            *out
              << "\nThe schur complement for the new KKT system will not be "
              << "updated yet in case it is singular!\n";
          }
        }
        else {
          assume_lin_dep_ja = false;
          try {
            if(act_set->add_constraint( ja, bnd_ja, false, out, output_level, true ))
              summary_lines_counter = 0;
            else {
              // Print end of row for rank if the right print level
              if( (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
                *out << right << setw(6) << "LI" << endl;
                out->flush();
                --summary_lines_counter;
              }
            }
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out << "\nNew KKT system is nonsingular! (linearly independent (LI) constraints)\n";
            }
            if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
              *out
                << "\nPrinting active set after adding constraint ja = " << ja
                << " ...\n";
              dump_act_set_quantities( *act_set, *out );
            }
          }
          catch( const MatrixSymAddDelUpdateable::SingularUpdateException& excpt ) {
            // Constraint really is linearly dependent.
            if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
              *out
                << "\n\nSchur complement update failed:\n"
                << "(" << excpt.what() << ")\n"
                << "\nConstraint ja = " << ja << " appears to be linearly dependent!\n\n";
            }
            summary_lines_counter = 0;
            if( !(act_set->q_D_hat() + act_set->q_plus_hat()) ) {
              // Print omsg.
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\nQPSchur::qp_algo(...) : "
                  << "Error, constraint j = "<< ja << " is linearly dependent "
                  << "and there are no other constraints to drop.\n"
                  << "The QP must be infeasible\n";
              }
              return INFEASIBLE_CONSTRAINTS;
            }
            assume_lin_dep_ja = true;
            schur_comp_update_failed = true;
          }
          catch( const MatrixSymAddDelUpdateable::WrongInertiaUpdateException& excpt ) {
            // Reduced Hessian has the wrong inertia
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\n\nSchur complement appears to have the wrong inertia:\n"
                << "(" << excpt.what() << ")\n"
                << "\nThe QP appears to be nonconvex.\n"
                << "\nWe have no choice but to terminate the primal-dual QP algorithm!\n";
            }
            return NONCONVEX_QP;
          }
        }
        if( assume_lin_dep_ja && can_ignore_ja ) {
          act_set->qp().constraints().ignore( ja );
          next_step = PICK_VIOLATED_CONSTRAINT;
          continue;
        }
        // Infer the sign of the multiplier for the new constraint being added
        beta = ( con_ja_val > b_a ? +1.0 : -1.0 );
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\nbeta = " << beta << endl;
        }
      }
      case COMPUTE_SEARCH_DIRECTION: {
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\n*** COMPUTE_SEARCH_DIRECTION\n";
        }
        const EtaVector e_ja( ja, n + m_breve );
        if( assume_lin_dep_ja ) {
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\nThe KKT system for the trial active set is assumed or known to be singular!\n"
              << "\nComputing the steps from:"
              << "\np_z_hat = inv(S_hat)*(-v_a+U_hat'*inv(Ko)*u_a), p_v = inv(Ko)*(-u_a-U_hat*p_z_hat) ...\n";
          }
          //
          // The schur complement is not updated so we must compute
          // p_z_hat and p_v explicitly.
          //
          // If all the degrees of freedom are used up then we know that the step for
          // the primal variables will be zero.  However, if m > 0 then v and p_v also
          // contain terms for the Lagrange multipliers for the equality constriants
          // but we don't need to compute these during the algorithm.
          // Therefore we can just set p_v = 0 and save a solve with Ko.
          // If the QP is feasible then a constraint will be dropped, the
          // KKT system will be updated and then v_plus will be computed
          // at the next iteration and be used to compute v so all is good.
          //
          const bool all_dof_used_up = act_set->all_dof_used_up();
          if( act_set->is_init_fixed(ja) ) {
            //
            // Fix a varaible that was fixed and then freed.
            //
            // [   Ko     U_hat ] [   p_v   ] = [   0  ]
            // [ U_hat'   V_hat ] [ p_z_hat ]   [ -v_a ]
            //
            // v_a = e(sa) <: R^q_hat            : where sa = s_map(-ja)
            //
            //        / 0                        : if x_init(ja) == bnd_ja
            // d_a =  |
            //        \ b_a - b_X(l_x_X_map(ja)) : otherwise
            //
            // p_z_hat = -inv(S_hat)*v_a
            //
            // p_v = inv(Ko)*(-U_hat*p_z_hat)
            //
            // gamma_plus = (d_a - v_a' * z_hat ) / ( v_a' * p_z_hat )
            //
            const size_type
              sa = act_set->s_map(-int(ja)),
              la = act_set->qp().l_x_X_map()(ja);
            TEUCHOS_TEST_FOR_EXCEPT( !( sa ) );
            TEUCHOS_TEST_FOR_EXCEPT( !( la ) );
            // v_a = e(sa) <: R^q_hat
            Workspace<value_type> v_a_ws(wss,act_set->q_hat());
            DVectorSlice v_a(&v_a_ws[0],v_a_ws.size());
            v_a = 0.0;
            v_a(sa) = 1.0;
            // d_a
            const value_type
              d_a = ( bnd_ja == qp.x_init()(ja) 
                  ? 0.0
                  : b_a - qp.b_X()(la) );
            // p_z_hat = -inv(S_hat) * v_a
            V_InvMtV( &act_set->p_z_hat(), act_set->S_hat(), no_trans, v_a() );
            Vt_S( &act_set->p_z_hat(), -1.0 );
            // p_v = inv(Ko)*(-U_hat*p_z_hat)
            if(!all_dof_used_up) {
              calc_v( qp.Ko(), NULL, act_set->U_hat(), act_set->p_z_hat(), &p_v() );
            }
            else {
              p_v = 0.0;
            }
            // Iterative refinement?
            if( using_iter_refinement ) {
              if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
                *out
                  << "\n\nPerforming iterative refinement on p_v, p_z_hat system ...\n";
              }
              summary_lines_counter = 0;
              // [   Ko     U_hat ] [   p_v   ] = [   0  ]
              // [ U_hat'   V_hat ] [ p_z_hat ]   [ -v_a ]
              EIterRefineReturn status = iter_refine(
                *act_set, out, output_level, 0.0, NULL, +1.0, &v_a
                ,&p_v, &act_set->p_z_hat()
                ,iter_refine_num_resid, iter_refine_num_solves
                );
            }
            // gamma_plus = ( d_a - v_a'*z_hat ) / ( v_a'*p_z_hat )
            if(!all_dof_used_up)
              gamma_plus = ( ( d_a - dot(v_a(),act_set->z_hat()) )
                       / ( dot(v_a(),act_set->p_z_hat()) ) );
            else
              gamma_plus = beta * inf;
          }
          else {
            //
            // Add a constraint that is not an initially fixed
            // variable bound.
            // 
            // [   Ko     U_hat ] [   p_v   ] = [ -u_a ]
            // [ U_hat'   V_hat ] [ p_z_hat ]   [ -v_a ]
            //
            // p_z_hat = inv(S_hat) * ( - v_a + U_hat' * inv(Ko) * u_a )
            // 
            // p_v = inv(Ko) * ( -u_a - U_hat * p_z_hat )
            // 
            // gamma_plus = ( d_a - u_a'*v - v_a'*z_hat ) / ( u_a'*p_v + v_a'*p_z_hat )
            // 
            // ToDo: (9/25/00): Make u_a and v_a both sparse and combine the following code.
            // 
            if( ja <= n ) {
              // Fix an initially free variable
              //
              // u_a = [ Q_R' * e(ja) ] <: R^(n_R+m)
              //       [      0       ]
              //
              // v_a = 0     <: R^(q_hat)
              //
              // d_a = b_a   <: R
              // 
              const size_type
                la = act_set->qp().Q_R().lookup_col_j(ja);
              TEUCHOS_TEST_FOR_EXCEPT( !(  la  ) );
              const EtaVector u_a = EtaVector(la,n_R+m);
              const value_type d_a = b_a;
              DVector t1;
              // t1 = inv(Ko) * u_a
              V_InvMtV( &t1, qp.Ko(), no_trans, u_a() );
              if( act_set->q_hat() ) {
                // t2 = U_hat'*t1
                DVector t2;
                V_MtV( &t2, act_set->U_hat(), trans, t1() );
                // p_z_hat = inv(S_hat) * t2
                V_InvMtV( &act_set->p_z_hat(), act_set->S_hat(), no_trans, t2() );
                // t1 = - u_a
                V_StV( &t1, -1.0, u_a() );
                // t1 += - U_hat * p_z_hat
                Vp_StMtV( &t1(), -1.0, act_set->U_hat(), no_trans, act_set->p_z_hat() );
                // p_v = inv(Ko) * t1
                if(!all_dof_used_up)
                  V_InvMtV( &p_v, qp.Ko(), no_trans, t1() );
                else
                  p_v = 0.0;
              }
              else {
                // p_v = -t1
                V_mV( &p_v, t1() );
              }
              // Iterative refinement?
              if( using_iter_refinement ) {
                if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
                  *out
                    << "\n\nPerforming iterative refinement on p_v, p_z_hat system ...\n";
                }
                summary_lines_counter = 0;
                // [   Ko     U_hat ] [   p_v   ] = [ -u_a ]
                // [ U_hat'   V_hat ] [ p_z_hat ]   [   0  ]
                Workspace<value_type> dense_u_a_ws(wss,u_a().dim());
                DVectorSlice dense_u_a(&dense_u_a_ws[0],dense_u_a_ws.size());
                dense_u_a = 0.0; // Make a dense copy of u_a!
                dense_u_a(u_a().begin()->index()+u_a().offset()) = 1.0;
                EIterRefineReturn status = iter_refine(
                  *act_set, out, output_level, +1.0, &dense_u_a, 0.0, NULL
                  ,&p_v, &act_set->p_z_hat()
                  ,iter_refine_num_resid, iter_refine_num_solves
                  );
                summary_lines_counter = 0;
              }
              // gamma_plus = ( d_a - u_a'*v) / ( u_a'*p_v )
              if(!all_dof_used_up)
                gamma_plus = ( d_a - dot(u_a(),*v) ) / dot(u_a(),p_v());
              else
                gamma_plus = beta * inf;
            }
            else {
              // Add a general inequality (or equality) constraint
              //
              // u_a = [ Q_R' * A_bar * e(ja) ] <: R^(n_R + m)
              //       [          0           ]
              // 
              // v_a = P_XF_hat' * A_bar * e_ja <: R^(q_hat)
              //
              // d_a = b_a - b_X' * (Q_X' * A_bar * e_ja) <: R
              //
              Workspace<value_type> u_a_ws( wss, n_R + m );
              DVectorSlice u_a( &u_a_ws[0], u_a_ws.size() );
              Workspace<value_type> v_a_ws( wss, act_set->q_hat() );
              DVectorSlice v_a( &v_a_ws[0], v_a_ws.size() );
              // u_a(1:n_R) =  Q_R' * A_bar * e(ja)
              Vp_StPtMtV( &u_a(1,n_R), 1.0, qp.Q_R(), trans
                    , qp.constraints().A_bar(), no_trans, e_ja(), 0.0 );
              // u_a(n_R+1:n_R+m) = 0.0
              if(m)
                u_a(n_R+1,n_R+m) = 0.0;
              // t0 = Q_X' * A_bar * e_ja
              Workspace<value_type> t0_ws( wss, n-n_R );
              DVectorSlice t0( &t0_ws[0], t0_ws.size() );
              if( n > n_R )
                Vp_StPtMtV( &t0(), 1.0, qp.Q_X(), trans
                      , qp.constraints().A_bar(), no_trans, e_ja(), 0.0 );
              // d_a = b_a - b_X'*t0
              const value_type
                d_a = b_a - ( n > n_R ? dot( qp.b_X(), t0() ) : 0.0 );
              // t1 = inv(Ko) * u_a
              Workspace<value_type> t1_ws( wss, n_R + m );
              DVectorSlice t1( &t1_ws[0], t1_ws.size() );
              V_InvMtV( &t1, qp.Ko(), no_trans, u_a );
              if( act_set->q_hat() ) {
                // t2 = U_hat'*t1
                Workspace<value_type> t2_ws( wss, act_set->q_hat() );
                DVectorSlice t2( &t2_ws[0], t2_ws.size() );
                V_MtV( &t2, act_set->U_hat(), trans, t1() );
                // v_a = P_XF_hat' * A_bar * e_ja
                Vp_StPtMtV( &v_a(), 1.0, act_set->P_XF_hat(), trans
                      , qp.constraints().A_bar(), no_trans, e_ja(), 0.0 );
                // t2 += -v_a
                Vp_StV( &t2(), -1.0, v_a() );
                // p_z_hat = inv(S_hat) * t2
                V_InvMtV( &act_set->p_z_hat(), act_set->S_hat(), no_trans, t2() );
                if(!all_dof_used_up) {
                  // t1 = - u_a
                  V_StV( &t1, -1.0, u_a() );
                    // t1 += - U_hat * p_z_hat
                  Vp_StMtV( &t1(), -1.0, act_set->U_hat(), no_trans, act_set->p_z_hat() );
                    // p_v = inv(Ko) * t1
                  V_InvMtV( &p_v, qp.Ko(), no_trans, t1() );
                }
                else {
                  p_v = 0.0;
                }
                // Iterative refinement?
                if( using_iter_refinement ) {
                  if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
                    *out
                      << "\n\nPerforming iterative refinement on p_v, p_z_hat system ...\n";
                  }
                  summary_lines_counter = 0;
                  // [   Ko     U_hat ] [   p_v   ] = [ -u_a ]
                  // [ U_hat'   V_hat ] [ p_z_hat ]   [ -v_a ]
                  EIterRefineReturn status = iter_refine(
                    *act_set, out, output_level, +1.0, &u_a, +1.0, act_set->q_hat() ? &v_a : NULL
                    ,&p_v, act_set->q_hat() ? &act_set->p_z_hat() : NULL
                    ,iter_refine_num_resid, iter_refine_num_solves
                    );
                  summary_lines_counter = 0;
                }
                // gamma_plus = ( d_a - u_a'*v - v_a'*z_hat ) / ( u_a'*p_v + v_a'*p_z_hat )
                if(!all_dof_used_up)
                  gamma_plus = ( ( d_a - dot(u_a,*v) - dot(v_a(),act_set->z_hat()) )
                           / ( dot(u_a,p_v()) + dot(v_a(),act_set->p_z_hat()) ) );
                else
                  gamma_plus = beta * inf;
              }
              else {
                // p_v = -t1
                if(!all_dof_used_up)
                  V_mV( &p_v, t1() );
                else
                  p_v = 0.0;
                // gamma_plus = ( d_a - u_a'*v) / ( u_a'*p_v )
                if(!all_dof_used_up)
                  gamma_plus = ( d_a - dot(u_a,*v) ) / dot(u_a,p_v());
                else
                  gamma_plus = beta * inf;
              }
            }
          }
          if( schur_comp_update_failed && gamma_plus * beta < 0 ) {
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << "\nThe schur complement update failed and gamma_plus = " << gamma_plus << " is the wrong sign"
                << "\nso we will assume the sign error for (...)/+-0 was due to arbitrary roundoff"
                << "\nand therefore we will set gamma_plus = -gamma_plus\n";
            }
            gamma_plus = -gamma_plus;
          }
          // Print the steps p_v and p_z_hat
          if(act_set->q_hat()) {
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << "\n||p_z_hat||inf = " << norm_inf(act_set->p_z_hat()) << endl;
            }
            if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
              *out << "\np_z_hat =\n" << act_set->p_z_hat();
            }
          }
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\n||p_v||inf = " << norm_inf(p_v()) << endl;
          }
          if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
            *out << "\np_v =\n" << p_v();
          }
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\ngamma_plus = " << gamma_plus << endl;
          }
          // Compute step for mu_D_hat
          if( act_set->q_D_hat() ) {
            // Compute for steps of all the constraints in the current active set
            calc_p_mu_D( *act_set, p_v(), act_set->p_z_hat(), &ja, &act_set->p_mu_D_hat() );
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << "\n||p_mu_D_hat||inf = " << norm_inf(act_set->p_mu_D_hat()) << std::endl;
            }
            if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
              *out
                << "\np_mu_D_hat =\n" << act_set->p_mu_D_hat();
            }
          }
        }
        else {
          // The new schur complement is already updated so compute
          // the solution outright.
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\nThe KKT system for the new active set is or known to be nonsingular!\n"
              << "\nComputing the steps from:"
              << "\nz_hat_plus = inv(S_hat)*(d_hat-U_hat'vo), v_plus = inv(Ko)*(fo-U_hat*z_hat_plus) ...\n";
          }
        
          // Compute z_hat_plus, v_plus

          // z_hat_plus = inv(S_hat) * ( d_hat - U_hat' * vo  )
          const size_type q_hat = act_set->q_hat();
          z_hat_plus.bind(DVectorSlice(&z_hat_plus_ws[0],q_hat));
          calc_z( act_set->S_hat(), act_set->d_hat(), act_set->U_hat(), &vo
            , &z_hat_plus );
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\n||z_hat_plus||inf = " << norm_inf(z_hat_plus()) << std::endl;
          }
          if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
            *out
              << "\nz_hat_plus =\n" << z_hat_plus();
          }
          // v_plus = inv(Ko) * (fo - U_hat * z_hat_plus)
          calc_v( qp.Ko(), &qp.fo(), act_set->U_hat(), z_hat_plus
            , &v_plus );
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\n||v_plus||inf = " << norm_inf(v_plus()) << std::endl;
          }
          if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
            *out
              << "\nv_plus =\n" << v_plus();
          }
          if( using_iter_refinement ) {
            if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
              *out
                << "\n\nPerforming iterative refinement on v_plus, z_hat_plus system ...\n";
            }
            summary_lines_counter = 0;
            // [   Ko     U_hat ] [   v_plus   ] = [   fo  ]
            // [ U_hat'   V_hat ] [ z_hat_plus ]   [ d_hat ]
            EIterRefineReturn status = iter_refine(
              *act_set, out, output_level, -1.0, &qp.fo(), -1.0, &act_set->d_hat()
              ,&v_plus, &z_hat_plus
              ,iter_refine_num_resid, iter_refine_num_solves
              );
          }
          // Compute p_z_hat (change in z_hat w.r.t newly added constriant multiplier)
          DVectorSlice p_z_hat = act_set->p_z_hat();
          // p_z_hat = z_hat_plus - z_hat
          V_VmV( &p_z_hat(), z_hat_plus(), act_set->z_hat() );
          // p_v = v_plus - v
          V_VmV( &p_v(), v_plus(), *v );
          // p_mu_D_hat
          if( act_set->q_D_hat() )
            calc_p_mu_D( *act_set, p_v(), p_z_hat(), NULL, &act_set->p_mu_D_hat() );
          // gamma_plus
          const size_type sa = act_set->s_map(ja);
          if(sa) {
            // This is not an initially fixed variable that returned to its
            // initial.  The multiplier for this constriant may not be the
            // last element if an ADD/DROP was performed on the last iteration
            // in order to get here where the DROP was an initially fixed variable
            // that was freed and therefore the KKT system was augmented so this
            // multiplier is not the last element of z_hat(...).
            gamma_plus = z_hat_plus(sa);
          }
          else {
            // This must be an initially fixed variable that returned to its
            // initial bound.  This will be the last element even if an ADD/DROP
            // was just performed since a drop would only remove elements from
            // p_mu_D_hat, not add them.
            gamma_plus = act_set->p_mu_D_hat()(act_set->q_D_hat());
          }
          // p_z_hat = p_z_hat / gamma_plus
          Vt_S( &p_z_hat(), 1.0 / gamma_plus );
          // p_v = p_v / gamma_plus
          Vt_S( &p_v(), 1.0 / gamma_plus );
          // Print gama_plus, p_z_hat, p_v and p_mu_D_hat
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\ngamma_plus = " << gamma_plus << std::endl;
          }
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\n||p_z_hat||inf = " << norm_inf(p_z_hat()) << std::endl;
          }
          if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
            *out
              << "\np_z_hat =\n" << p_z_hat();
          }
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\n||p_v||inf = " << norm_inf(p_v()) << std::endl;
          }
          if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
            *out
              << "\np_v =\n" << p_v();
          }
          if( act_set->q_D_hat() ) {
            // p_mu_D_hat = p_mu_D_hat / gamma_plus
            Vt_S( &act_set->p_mu_D_hat(), 1.0 / gamma_plus ); 
            if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
              *out
                << "\n||p_mu_D_hat||inf =\n" << norm_inf(act_set->p_mu_D_hat()) << std::endl;
            }
            if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES ) {
              *out
                << "\np_mu_D_hat =\n" << act_set->p_mu_D_hat();
            }
          }
        }
      }
      case COMPUTE_STEP_LENGTHS: {

        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\n*** COMPUTE_STEP_LENGTHS\n";
        }
        // Compute the dual infeasibility scaling
        const size_type q_hat = act_set->q_hat();
        dual_infeas_scale = 1.0;
//				if( q_hat )
//					dual_infeas_scale = my_max( dual_infeas_scale, norm_inf( act_set->z_hat() ) );
//				if( m )
//					dual_infeas_scale = my_max( dual_infeas_scale, norm_inf( (*v)(n_R+1,n_R+m) ) );
//				if( act_set->q_D_hat() )
//					dual_infeas_scale = my_max( dual_infeas_scale, norm_inf( act_set->mu_D_hat() ) );
        
        // Primal step length, t_P = beta * gamma_plus, z_plus = [ z_hat_plus; gama_plus ].
        // Or constraint ja is linearly dependent in which case p_x is zero so
        // t_P is infinite.
        t_P = beta * gamma_plus;	// Could be < 0
        if( t_P < 0.0 ) {
          if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
            *out
              << "\nWarning, A near degenerate inequality constraint ja = " << ja
              << " is being added that has the wrong sign with:\n"
              << "    t_P                     = " << t_P 					<< std::endl
              << "    dual_infeas_scale       = " << dual_infeas_scale	<< std::endl
              << "    norm_2_constr           = " << norm_2_constr  		<< std::endl
              << "    |t_P/(norm_2_constr*dual_infeas_scale)| = "
              << std::fabs(t_P/(norm_2_constr*dual_infeas_scale))
              << " <= dual_infeas_tol = " << dual_infeas_tol()			<< std::endl;
          }
          summary_lines_counter = 0;
          if( !using_iter_refinement ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out << "We are not using iterative refinement yet so turn it on"
                 << "\nthen recompute the steps ...\n";
            }
            using_iter_refinement = true;
            next_step = COMPUTE_SEARCH_DIRECTION;
            continue;
          }
          else {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "We are already using iterative refinement so the QP algorithm is terminated!\n";
            }
            return DUAL_INFEASIBILITY;
          }
        }
        t_P = beta * gamma_plus;	// Now guaranteed to be > 0
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\nt_P = " << t_P << endl;
        }

        /////////////////////////////////////////////////////////////////////////
        // Dual step length.  Largest step t that does not cause violation in
        // dual feasibility (i.e. lagrange multipliers for inequalities are
        // dual feasible, or primal optimal ).
        // lambda_hat_new = lambda_hat + beta * t_D * p_lambda_hat must be dual feasible.
        t_D = inf;
        jd = 0;
        value_type max_feas_viol = 0.0; // Remember the amount of violation.
        int j_degen = 0;	// remember which (if any) constraint was near
                  // degenerate and had an incorrect sign.
        EBounds	    bnd_jd;	// The bound of the constraint to be dropped.

        // Search through Lagrange multipliers in z_hat
        if( act_set->q_hat() ) {
          DVectorSlice z_hat = act_set->z_hat();
          DVectorSlice p_z_hat = act_set->p_z_hat();
          DVectorSlice::iterator
            z_itr		= z_hat.begin(),
            p_z_itr		= p_z_hat.begin();
          const size_type
            qq = assume_lin_dep_ja || (!assume_lin_dep_ja && return_to_init_fixed)
              ? q_hat : q_hat - 1;
          // Print header for s, j, z_hat(s), p_z_hat(s), bnds(s), t, t_D, jd
          if( qq > 0 && (int)output_level >= (int)OUTPUT_ACT_SET ) {
              *out
              << "\nComputing the maximum step for multiplers for dual feasibility\n\n"
              << right << setw(5)	<< "s"
              << right << setw(5)	<< "j"
              << right << setw(dbl_w)	<< "z_hat"
              << right << setw(dbl_w)	<< "p_z_hat"
              << right << setw(20)	<< "bnd"
              << right << setw(dbl_w)	<< "t"
              << right << setw(dbl_w)	<< "t_D"
              << right << setw(5)	<< "jd"	<< endl
              << right << setw(5)	<< "----"
              << right << setw(5)	<< "----"
              << right << setw(dbl_w)	<< "--------------"
              << right << setw(dbl_w)	<< "--------------"
              << right << setw(20)	<< "--------------"
              << right << setw(dbl_w)	<< "--------------"
              << right << setw(dbl_w)	<< "--------------"
              << right << setw(5)	<< "----"	<< endl;
          }
          for( int s = 1; s <= qq; ++s, ++z_itr, ++p_z_itr) {
            int j = act_set->ij_map(s);
            if( j > 0 ) {
              namespace ns = QPSchurPack;
              EBounds bnd = act_set->bnd(s);
              // Print first part of row for s, j, z_hat(s), p_z_hat(s), bnds(s) ....
              if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
                *out
                  << right << setw(5)	<< s
                  << right << setw(5)	<< j
                  << right << setw(dbl_w)	<< *z_itr
                  << right << setw(dbl_w)	<< *p_z_itr
                  << right << setw(20)	<< bnd_str(bnd);
              }
              value_type t = inf;
              // Lookout for degeneracy.
              bool j_is_degen = false;
              value_type viol;
              const int dual_feas_status
                = correct_dual_infeas(
                  j,bnd,t_P,1.0,dual_infeas_tol(),DEGENERATE_MULT
                  ,out,output_level,true,"z_hat(s)",&(*z_itr),&viol
                  ,"p_z_hat(s)",&(*p_z_itr),"z_hat_plus(s)"
                  , (assume_lin_dep_ja ? NULL: &z_hat_plus(s) ) );
              if( dual_feas_status < 0 ) {
                if( !using_iter_refinement ) {
                  if( (int)output_level >= (int)OUTPUT_BASIC_INFO )
                    *out << "We are not using iterative refinement yet so turn it on"
                       << "\nthen recompute the steps ...\n";
                  using_iter_refinement = true;
                  next_step = COMPUTE_SEARCH_DIRECTION;
                  continue;
                }
                else {
                  if( (int)output_level >= (int)OUTPUT_BASIC_INFO )
                    *out << "We are already using iterative refinement so the QP algorithm is terminated!\n";
                  return DUAL_INFEASIBILITY;
                }
              }
              else if( dual_feas_status == 0 ) {
                j_is_degen = true;
              }
              // If we get here either the dual variable was feasible or it
              // was near degenerate and was corrected!
              const value_type feas_viol = beta*(*p_z_itr);
              if( bnd == EQUALITY )
                ; // We don't care
              else if( bnd == LOWER && feas_viol <= 0.0 )
                ;	// dual feasible for all t > 0
              else if( bnd == UPPER && feas_viol >= 0.0 )
                ;	// dual feasible for all t > 0
              else {
                // finite t.
                t = -beta*(*z_itr)/(*p_z_itr);
                if( t < t_D ) {	// remember minimum step length
                  t_D = t;
                  jd = j;
                  if(j_is_degen) j_degen = j;
                  max_feas_viol = feas_viol;
                  bnd_jd = bnd;
                }
              }
              // Print rest of row for ... t, t_D, jd
              if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
                *out
                  << right << setw(dbl_w)	<< t
                  << right << setw(dbl_w)	<< t_D
                  << right << setw(5)	<< jd	<< endl;
              }
            }
          }
        }
          
        // Search through Lagrange multipliers in mu_D_hat
        if( act_set->q_D_hat() ) {
          const QPSchurPack::QP::x_init_t     &x_init    = qp.x_init();
          const QPSchurPack::QP::i_x_X_map_t  &i_x_X_map = qp.i_x_X_map();
          const size_type q_D_hat = act_set->q_D_hat();
          DVectorSlice mu_D_hat = act_set->mu_D_hat();
          DVectorSlice p_mu_D_hat = act_set->p_mu_D_hat();
          const size_type
            qD = assume_lin_dep_ja && return_to_init_fixed ? q_D_hat-1 : q_D_hat;
          // Print header for k, i, mu_D_hat(k), p_mu_D_hat(k), x_init(k), t, t_D, jd
          if( qD > 0 && (int)output_level >= (int)OUTPUT_ACT_SET ) {
            *out
              << "\nComputing the maximum step for multiplers for dual feasibility\n\n"
              << right << setw(5)	<< "k"
              << right << setw(5)	<< "i"
              << right << setw(dbl_w)	<< "mu_D_hat"
              << right << setw(dbl_w)	<< "p_mu_D_hat"
              << right << setw(20)	<< "x_init"
              << right << setw(dbl_w)	<< "t"
              << right << setw(dbl_w)	<< "t_D"
              << right << setw(5)	<< "jd"	<< endl
              << right << setw(5)	<< "----"
              << right << setw(5)	<< "----"
              << right << setw(dbl_w)	<< "--------------"
              << right << setw(dbl_w)	<< "--------------"
              << right << setw(20)	<< "--------------"
              << right << setw(dbl_w)	<< "--------------"
              << right << setw(dbl_w)	<< "--------------"
              << right << setw(5)	<< "----"	<< endl;
          }
          GenPermMatrixSlice::const_iterator
            Q_XD_itr = act_set->Q_XD_hat().begin(),
            Q_XD_end = Q_XD_itr + qD;
          for( ; Q_XD_itr != Q_XD_end; ++Q_XD_itr ) {
            const size_type k = Q_XD_itr->col_j();
            const size_type i = Q_XD_itr->row_i();
            DVectorSlice::iterator
              mu_D_itr		= mu_D_hat.begin() + (k-1),
              p_mu_D_itr		= p_mu_D_hat.begin() + (k-1);
            const size_type l = act_set->l_fxfx(k);
            EBounds bnd = qp.x_init()(i);
            // Print first part of row for s, j, z_hat(s), p_z_hat(s), bnds(s) ....
            if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
              *out
                << right << setw(5)	<< k
                << right << setw(5)	<< i
                << right << setw(dbl_w)	<< *mu_D_itr
                << right << setw(dbl_w)	<< *p_mu_D_itr
                << right << setw(20)	<< bnd_str(bnd);
            }
            value_type t = inf;
            // Lookout for degeneracy.
            bool j_is_degen = false;
            value_type viol;
            const int dual_feas_status
              = correct_dual_infeas(
                i,bnd,t_P,1.0,dual_infeas_tol(),DEGENERATE_MULT
                ,out,output_level,true,"mu_D_hat(k)",&(*mu_D_itr),&viol
                ,"p_mu_D_hat(k)",&(*p_mu_D_itr) );
            if( dual_feas_status < 0 ) {
              if( !using_iter_refinement ) {
                if( (int)output_level >= (int)OUTPUT_BASIC_INFO )
                  *out << "We are not using iterative refinement yet so turn it on"
                     << "\nthen recompute the steps ...\n";
                using_iter_refinement = true;
                next_step = COMPUTE_SEARCH_DIRECTION;
                continue;
              }
              else {
                if( (int)output_level >= (int)OUTPUT_BASIC_INFO )
                  *out << "We are already using iterative refinement so the QP algorithm is terminated!\n";
                return DUAL_INFEASIBILITY;
              }
            }
            else if( dual_feas_status == 0 ) {
              j_is_degen = true;
            }
            // If we get here either the dual variable was feasible or it
            // was near degenerate and was corrected!
            const value_type feas_viol = beta*(*p_mu_D_itr);
            if( bnd == EQUALITY )
              ; // We don't care
            else if( bnd == LOWER && feas_viol <= 0.0 )
              ;	// dual feasible for all t > 0
            else if( bnd == UPPER && feas_viol >= 0.0 )
              ;	// dual feasible for all t > 0
            else {
              // finite t.
              t = -beta*(*mu_D_itr)/(*p_mu_D_itr);
              if( t < t_D ) {	// remember minimum step length
                t_D = t;
                jd = -i;
                if(j_is_degen) j_degen = jd;
                max_feas_viol = feas_viol;
                bnd_jd = bnd;
                }
            }
            // Print rest of row for ... t, t_D, jd
            if( (int)output_level >= (int)OUTPUT_ACT_SET ) {
              *out
                << right << setw(dbl_w)	<< t
                << right << setw(dbl_w)	<< t_D
                << right << setw(5)	    << jd   << endl;
            }
          }
        }
        // Print t_D, jd
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\nt_D = "	<< t_D	<< endl
            << "jd = "		<< jd	<< endl;
        }
        if( jd == j_degen && jd != 0 && t_D < t_P ) {
          if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
            *out
              << "\nWarning, the near degenerate constraint j = "
              << jd << " which had the incorrect sign\nand was adjusted "
              << "was selected to be dropped from the active set.\n";
          }
        }
        // Print end of row for rank if the right print level
        if( assume_lin_dep_ja && !schur_comp_update_failed && (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
          if( t_P < huge_primal_step() )
            *out << right << setw(6) << "LI" << endl;
          else
            *out << right << setw(6) << "LD" << endl;
          out->flush();
          --summary_lines_counter;
        }
        // Print start of row for itr, q_hat, q(+), q_D, q_C, q_F, change, type, index, bound, violation
        if( t_D < t_P && (int)output_level == (int)OUTPUT_ITER_SUMMARY ) {
          *out
            << right << setw(6)	<< itr                   // itr
            << right << setw(6)	<< act_set->q_hat()      // q_hat
            << right << setw(6)	<< act_set->q_plus_hat() // q(+)
            << right << setw(6)	<< act_set->q_F_hat()    // q_F
            << right << setw(6)	<< act_set->q_C_hat()    // q_C
            << right << setw(6)	<< act_set->q_D_hat()    // q_D
            << right << setw(8)	<< "DROP"                // change
            << right << setw(9);                         // type
          if( jd < 0 ) {
            *out << "X_F";
          }
          else if( jd <= n ) {
            if( bnd_jd == qp.x_init()(jd) )
              *out << "X_F_C_F";
            else
              *out << "R_X_R";
          }
          else {
            *out << "GEN";
          }
          *out
            << right << setw(6)		<< jd                   // index
            << right << setw(10)	<< bnd_str(bnd_jd)      // bound
            << right << setw(dbl_w)	<< max_feas_viol        // violation
            << right << setw(6)		<< "LI" << endl;        // rank (this should be true!)
        }
      }
      case TAKE_STEP: {
        if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
          *out
            << "\n*** TAKE_STEP\n";
        }
        if( t_P >= huge_primal_step() && t_D >= huge_dual_step() ) {
          if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
            *out
              << "Error, QP is infeasible, inconsistent constraint a("<<ja<<") detected\n";
          }
          if( using_iter_refinement ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "We are already using iterative refinement so the QP algorithm is terminated!\n";
            }
            return INFEASIBLE_CONSTRAINTS;
          }
          else {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out << "We are not using iterative refinement yet so turn it on";
              if(assume_lin_dep_ja)
                *out << "\nthen pick another violated constriant to add ... \n";
              else
                *out << "\nthen recompute the steps ...\n";
            }
            summary_lines_counter = 0;
            last_jd = 0;  // erase this memory!
            last_ja = 0;  // ..
            using_iter_refinement = true;
            if(assume_lin_dep_ja) {
              EIterRefineReturn status = iter_refine(
                *act_set, out, output_level, -1.0, &qp.fo(), -1.0, act_set->q_hat() ? &act_set->d_hat() : NULL
                ,v, act_set->q_hat() ? &act_set->z_hat() : NULL
                ,iter_refine_num_resid, iter_refine_num_solves
                );
              next_step = PICK_VIOLATED_CONSTRAINT;
            }
            else {
              // Iterative refinement will be performed there
              next_step = COMPUTE_SEARCH_DIRECTION;
            }
            continue;
          }
        }
        else if( t_P > t_D ) {
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            if( t_P >= huge_primal_step() ) {
              *out
                << "\n*** (b) Dual Step (t_P = " << t_P << " >= huge_primal_step = "
                  << huge_primal_step() << ")\n";
            }
            else {
              *out
                << "\n*** (b) Partial Primal-Dual Step\n";
            }
          }
          // Check for cycling
          if( ja == last_jd && jd == last_ja ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\n\nQPSchur::qp_algo(...) : Error, the constraint "
                << "a(" << ja << ") with violation\n"
                << "(a(ja)'*x - b_a) = (" << con_ja_val
                << " - " << b_a << ") = " << (con_ja_val - b_a) << "\n"
                << "we are adding to the active set and the constraint constriant\n"
                << "a(" << jd << ") we are dropping were just dropped and added respectively!\n"
                << "The algorithm is cycling!\n";
            }
            if( using_iter_refinement ) {
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "We are already using iterative refinement so the QP algorithm is terminated!\n";
              }
              return SUBOPTIMAL_POINT;
            }
            else {
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out << "We are not using iterative refinement yet so turn it on";
                if(assume_lin_dep_ja)
                  *out << "\nthen pick another violated constriant to add ... \n";
                else
                  *out << "\nthen recompute the steps ...\n";
              }
              summary_lines_counter = 0;
              last_jd = 0;  // erase this memory!
              last_ja = 0;  // ..
              using_iter_refinement = true;
              if(assume_lin_dep_ja) {
                EIterRefineReturn status = iter_refine(
                  *act_set, out, output_level, -1.0, &qp.fo(), -1.0, act_set->q_hat() ? &act_set->d_hat() : NULL
                  ,v, act_set->q_hat() ? &act_set->z_hat() : NULL
                  ,iter_refine_num_resid, iter_refine_num_solves
                  );
                next_step = PICK_VIOLATED_CONSTRAINT;
              }
              else {
                // Iterative refinement will be performed there
                next_step = COMPUTE_SEARCH_DIRECTION;
              }
              continue;
            }
          }
          // Update the augmented KKT system
          try {
            if( assume_lin_dep_ja ) {
              if(act_set->drop_add_constraints( jd, ja, bnd_ja, true, out, output_level ))
                summary_lines_counter = 0;
            }
            else {
              if(act_set->drop_constraint( jd, out, output_level, true, true ))
                summary_lines_counter = 0;
            }
          }
          catch( const MatrixSymAddDelUpdateable::SingularUpdateException& excpt ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\n\nSchur complement appears to be singular and should not be:\n"
                << excpt.what()
                << "\nThe QP appears to be nonconvex and we therefore terminate the primal-dual QP algorithm!\n";
            }
            return NONCONVEX_QP;
          }
          catch( const MatrixSymAddDelUpdateable::WrongInertiaUpdateException& excpt ) {
            if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
              *out
                << "\n\nSchur complement appears to have the wrong inertia:\n"
                << excpt.what()
                << "\nThe QP appears to be nonconvex and we therefore terminate the primal-dual QP algorithm!\n";
            }
            return NONCONVEX_QP;
          }
          // z_hat = z_hat + beta * t_D * p_z_hat
          if(act_set->q_hat())
            Vp_StV( &act_set->z_hat(), beta * t_D, act_set->p_z_hat() );
          // v = v + beta * t_D * p_v
          Vp_StV( v, beta * t_D, p_v() );
          // mu_D_hat = mu_D_hat + beta * t_D * p_mu_D_hat
          if(act_set->q_D_hat())
            Vp_StV( &act_set->mu_D_hat(), beta * t_D, act_set->p_mu_D_hat() );

          ++(*num_drops);

          if( (int)output_level >= (int)OUTPUT_ITER_STEPS )
          {
            *out
              << "\nUpdated primal and dual variables:\n"
              << "\n||v||inf           = " << norm_inf(*v) << endl;
            if(act_set->q_hat()) {
              *out
                << "||z_hat||inf       = " << norm_inf(act_set->z_hat()) << endl;
            }
            if(act_set->q_D_hat()) {
              *out
                << "max(|mu_D_hat(i)|) = " << norm_inf(act_set->mu_D_hat()) << endl
                << "min(|mu_D_hat(i)|) = " << min_abs(act_set->mu_D_hat()) << endl;
            }
          }
          if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES )
          {
            *out << "\nv = \n" << *v << endl;
            if(assume_lin_dep_ja) {
              *out
                << "\nPrinting active set after dropping constraint jd = " << jd
                << " and adding constraint ja = " << ja << " ...\n";
            }
            else {
              *out
                << "\nPrinting active set after dropping constraint jd = " << jd
                << " ...\n";
            }
            dump_act_set_quantities( *act_set, *out );
          }
          last_jd = jd;
          assume_lin_dep_ja = false;  // If we get here then we know these are true!
          next_step = COMPUTE_SEARCH_DIRECTION;
          continue;
        }
        else {	// t_P < t_D
          if( (int)output_level >= (int)OUTPUT_ITER_STEPS ) {
            *out
              << "\n*** (c) Full Primal-Dual Step\n";
          }
          ++(*num_adds);
          if( !assume_lin_dep_ja ) {
            act_set->z_hat() 	= z_hat_plus;
            *v 					= v_plus;
          }
          else {
            bool threw_exception = false;
            try {
              if(act_set->add_constraint( ja, bnd_ja, true, out, output_level, true, true ))
                summary_lines_counter = 0;
            }
            catch( const MatrixSymAddDelUpdateable::SingularUpdateException& excpt ) {
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\n\nSchur complement appears to be singular and should not be:\n"
                  << excpt.what() << std::endl;
              }
              threw_exception = true;
            }
            catch( const MatrixSymAddDelUpdateable::WrongInertiaUpdateException& excpt ) {
              if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                *out
                  << "\n\nSchur complement appears to have the wrong inertia:\n"
                  << excpt.what() << std::endl;
              }
              threw_exception = true;
            }
            if( threw_exception ) {
              if( !using_iter_refinement ) {
                if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                  *out << "We are not using iterative refinement yet so turn it on and\n"
                     << "go back and pick a new violated constraint to add to the active set ...\n";
                }
                using_iter_refinement = true;
                if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
                  *out
                    << "\n\nPerforming iterative refinement on v, z_hat system ...\n";
                }
                summary_lines_counter = 0;
                // [   Ko     U_hat ] [   v   ] = [   fo  ]
                // [ U_hat'   V_hat ] [ z_hat ]   [ d_hat ]
                EIterRefineReturn status = iter_refine(
                  *act_set, out, output_level, -1.0, &qp.fo(), -1.0, &act_set->d_hat()
                  ,v, &act_set->z_hat()
                  ,iter_refine_num_resid, iter_refine_num_solves
                  );
                next_step = PICK_VIOLATED_CONSTRAINT;
                continue;
              }
              else {
                if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
                  *out << "Darn, we are already using iterative refinement!"
                     << "\nThe QP appears to be nonconvex and we therefore terminate the primal-dual QP algorithm!\n";
                }
                return NONCONVEX_QP;
              }
            }
            // z_hat = z_hat + beta * t_P * p_z_hat
            if(act_set->q_hat())
              Vp_StV( &act_set->z_hat(), beta * t_P, act_set->p_z_hat() );
            // v = v + beta * t_P * p_v
            Vp_StV( v, beta * t_P, p_v() );
          }
          // mu_D_hat = mu_D_hat + beta * t_P * p_mu_D_hat
          if(act_set->q_D_hat())
            Vp_StV( &act_set->mu_D_hat(), beta * t_P, act_set->p_mu_D_hat() );


          if( (int)output_level >= (int)OUTPUT_ITER_STEPS )
          {
            *out
              << "\nUpdated primal and dual variables:\n"
              << "\n||v||inf           = " << norm_inf(*v) << endl;
            if(act_set->q_hat()) {
              *out
                << "||z_hat||inf       = " << norm_inf(act_set->z_hat()) << endl;
            }
            if(act_set->q_D_hat()) {
              *out
                << "max(|mu_D_hat(i)|) = " << norm_inf(act_set->mu_D_hat()) << endl
                << "min(|mu_D_hat(i)|) = " << min_abs(act_set->mu_D_hat()) << endl;
            }
          }
          if( (int)output_level >= (int)OUTPUT_ITER_QUANTITIES )
          {
            *out << "\nv = \n" << *v << endl;
            if( assume_lin_dep_ja ) {
              *out
                << "\nPrinting active set after adding constraint ja = " << ja
                << " ...\n";
              dump_act_set_quantities( *act_set, *out );
            }
            else {
              if(act_set->q_hat())
                *out << "\nz_hat =\n" << act_set->z_hat();
              if(act_set->q_D_hat())
                *out << "\nmu_D_hat =\n" << act_set->mu_D_hat();
            }
          }
          assume_lin_dep_ja = false;  // If we get here then we know these are true!
          next_step = PICK_VIOLATED_CONSTRAINT;
          continue;
        }
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true);	// only a local programming error
    }
  }

  } // end try
  catch( std::exception& excpt ) {
    if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
      *out
        << "\n\n*** Caught a standard exception :\n"
        << excpt.what()
        << "\n*** Rethrowing the exception ...\n";
    }
    throw;
  }
  catch(...) {
    if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
      *out
        << "\n\n*** Caught an unknown exception.  Rethrowing the exception ...\n";
    }
    throw;
  }
  // If you get here then the maximum number of QP iterations has been exceeded
  return MAX_ITER_EXCEEDED;
}

void QPSchur::set_x( const ActiveSet& act_set, const DVectorSlice& v, DVectorSlice* x )
{
  using BLAS_Cpp::no_trans;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_MtV;
  
  // x = Q_R * v(1:n_R) + Q_X * b_X + P_XF_hat * z_hat
  V_MtV( x, act_set.qp().Q_R(), no_trans, v(1,act_set.qp().n_R()) );
  if( act_set.qp().n() > act_set.qp().n_R() )
    Vp_MtV( x, act_set.qp().Q_X(), no_trans, act_set.qp().b_X() );
  if( act_set.q_F_hat() )
    Vp_MtV( x, act_set.P_XF_hat(), no_trans, act_set.z_hat() );
}

void QPSchur::set_multipliers(
  const ActiveSet& act_set, const DVectorSlice& v
  ,SpVector* mu, DVectorSlice* lambda, SpVector* lambda_breve
  )
{
  using BLAS_Cpp::no_trans;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_MtV;
  using LinAlgOpPack::Vp_StMtV;
  using LinAlgOpPack::V_MtV;
  using LinAlgOpPack::Vp_MtV;
  using AbstractLinAlgPack::V_MtV;
  using AbstractLinAlgPack::Vp_MtV;
  namespace GPMSTP = AbstractLinAlgPack::GenPermMatrixSliceIteratorPack;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const size_type
    n       = act_set.qp().n(),
    n_R     = act_set.qp().n_R(),
    m       = act_set.qp().m(),
    m_breve = act_set.qp().constraints().m_breve(),
    q_hat   = act_set.q_hat();
  const QPSchurPack::QP::x_init_t
    &x_init = act_set.qp().x_init();

  //
  // mu = P_plus_hat(1:n,:) * z_hat + Q_XD_hat * mu_D + (steps for initially fixed
  // 		variables fixed to the other bounds)
  //
  // lambda_breve = P_plus_hat(n+1:n+m_breve,:) * z_hat
  //
  typedef SpVector::element_type ele_t;
  mu->resize( n, n-m );	                // Resize for the maxinum number
  lambda_breve->resize( m_breve, n-m );   // of active constraints possible.
  // mu += Q_XD_hat * mu_D_hat
  if( act_set.q_D_hat() )
    Vp_MtV( mu, act_set.Q_XD_hat(), no_trans, act_set.mu_D_hat() );
  // Set all the multipliers in z_hat
  if(q_hat){
    const DVectorSlice
      z_hat = act_set.z_hat();
    for( size_type s = 1; s <= q_hat; ++s ) {
      const int ij = act_set.ij_map(s);
      if(ij > 0) {
        const size_type j = ij;
        if( j <= n )
          mu->add_element(ele_t(j,z_hat(s)));
        else
          lambda_breve->add_element(ele_t(j-n,z_hat(s)));
      }
    }
  }
  mu->sort();
  lambda_breve->sort();
  // lambda = v(n_R+1,n_R+m)
  if( m ) {
    *lambda = v(n_R+1,n_R+m);
  }
}

bool QPSchur::timeout_return( StopWatchPack::stopwatch* timer, std::ostream *out, EOutputLevel output_level ) const
{
  const value_type minutes = timer->read() / 60;
  if( minutes >= max_real_runtime() ) {
    if( (int)output_level >= (int)OUTPUT_BASIC_INFO ) {
      *out
        << "\n*** Runtime = " << minutes << " min >= max_real_runtime = " << max_real_runtime() << " min!\n"
        << "Must terminite the algorithm!\n";
    }
    return true;
  }
  return false;
}

QPSchur::EIterRefineReturn
QPSchur::iter_refine(
  const ActiveSet      &act_set
  ,std::ostream        *out
  ,EOutputLevel        output_level
  ,const value_type    ao
  ,const DVectorSlice  *bo
  ,const value_type    aa
  ,const DVectorSlice  *ba
  ,DVectorSlice        *v
  ,DVectorSlice        *z
  ,size_type           *iter_refine_num_resid
  ,size_type           *iter_refine_num_solves
  )
{
  using std::endl;
  using std::setw;
  using std::left;
  using std::right;
  using BLAS_Cpp::no_trans;
  using BLAS_Cpp::trans;
  using DenseLinAlgPack::norm_inf;
  using DenseLinAlgPack::Vp_StV;
  using LinAlgOpPack::Vp_V;
  using LinAlgOpPack::Vp_StMtV;
  using LinAlgOpPack::V_InvMtV;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  typedef DenseLinAlgPack::value_type           extra_value_type;

  const value_type small_num = std::numeric_limits<value_type>::min();

  const int int_w = 8;
  const char int_ul[] = "------";
  const int dbl_min_w = 20;
  const int dbl_w = ( out ? my_max(dbl_min_w,int(out->precision()+8)): 20 );
  const char dbl_ul[] = "------------------";

  const QPSchurPack::QP
    &qp    = act_set.qp();
  const MatrixSymOpNonsing
    &Ko    = qp.Ko(),
    &S_hat = act_set.S_hat();
  const MatrixOp
    &U_hat = act_set.U_hat();
  const DenseLinAlgPack::size_type
    n          = qp.n(),
    n_R        = qp.n_R(),
    m          = qp.m(),
    q_hat      = act_set.q_hat();
  const DVectorSlice
    fo    = qp.fo(),
    d_hat = (q_hat ? act_set.d_hat() : DVectorSlice());

  Workspace<extra_value_type>
    ext_ro_ws(wss,n_R+m),
    ext_ra_ws(wss,q_hat);
  DenseLinAlgPack::VectorSliceTmpl<extra_value_type>
    ext_ro(&ext_ro_ws[0],ext_ro_ws.size()),
    ext_ra(ext_ra_ws.size()?&ext_ra_ws[0]:NULL,ext_ra_ws.size());
  Workspace<value_type>
    ro_ws(wss,n_R+m),
    ra_ws(wss,q_hat),
    t1_ws(wss,n_R+m),
    del_v_ws(wss,n_R+m),
    del_z_ws(wss,q_hat),
    v_itr_ws(wss,n_R+m),
    z_itr_ws(wss,q_hat);
  DVectorSlice
    ro(&ro_ws[0],ro_ws.size()),
    ra(ra_ws.size()?&ra_ws[0]:NULL,ra_ws.size()),
    t1(&t1_ws[0],t1_ws.size()),
    del_v(&del_v_ws[0],del_v_ws.size()),
    del_z(del_z_ws.size()?&del_z_ws[0]:NULL,del_z_ws.size()),
    v_itr(&v_itr_ws[0],v_itr_ws.size()),
    z_itr(z_itr_ws.size()?&z_itr_ws[0]:NULL,z_itr_ws.size());
  
  // Accumulate into temporary variables
  v_itr = *v;
  if(q_hat)
    z_itr = *z;

  // Print summary header
  if( out && (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
    *out
      << "\nBeginning iterative refinement ..."
      << "\niter_refine_opt_tol = " << iter_refine_opt_tol()
      << ", iter_refine_feas_tol = " << iter_refine_feas_tol()
      << "\niter_refine_min_iter = " << iter_refine_min_iter()
      << ", iter_refine_max_iter = " << iter_refine_max_iter() << "\n\n";
    //
    *out
      << right << setw(int_w) << "ir_itr"
      << right << setw(dbl_w) << "roR_scaling"
      << right << setw(dbl_w) << "||roR||s"
      << left  << setw(1)     << " ";
    if(m) {
      *out
        << right << setw(dbl_w) << "rom_scaling"
        << right << setw(dbl_w) << "||rom||s"
        << left  << setw(1)     << " ";
    }
    if(q_hat) {
      *out
        << right << setw(dbl_w) << "ra_scaling"
        << right << setw(dbl_w) << "||ra||s"
        << left  << setw(1)     << " ";
    }
    *out
      << right << setw(dbl_w) << "||del_v||/||v||inf"
      << right << setw(dbl_w) << "||del_z||/||z||inf"
      << endl;
    //
    *out
      << right << setw(int_w) << int_ul
      << right << setw(dbl_w) << dbl_ul
      << right << setw(dbl_w) << dbl_ul
      << left  << setw(1)     << " ";
    if(m) {
      *out
      << right << setw(dbl_w) << dbl_ul
      << right << setw(dbl_w) << dbl_ul
      << left  << setw(1)     << " ";
    }
    if(q_hat) {
      *out
      << right << setw(dbl_w) << dbl_ul
      << right << setw(dbl_w) << dbl_ul
      << left  << setw(1)     << " ";
    }
    *out
      << right << setw(dbl_w) << dbl_ul
      << right << setw(dbl_w) << dbl_ul
      << endl;
  }
  //
  // Perform iterative refinement iterations
  //
  EIterRefineReturn return_status = ITER_REFINE_NOT_PERFORMED;
  value_type
    roR_nrm_o,   rom_nrm_o,   ra_nrm_o,
    roR_nrm,     rom_nrm,     ra_nrm;
  for( size_type iter_refine_k = 0;
     iter_refine_k < iter_refine_max_iter() && return_status != ITER_REFINE_CONVERGED;
     ++iter_refine_k)
  {
    //
    // Compute the residual (in extended precision?)
    //
    // [ ro ] = [   Ko     U_hat ] [ v ] + [ ao*bo ]
    // [ ra ]   [ U_hat'   V_hat ] [ z ]   [ aa*ba ]
    //
    value_type
      roR_scaling = 0.0,
      rom_scaling = 0.0,
      ra_scaling  = 0.0;
    ++(*iter_refine_num_resid);
    calc_resid(
      act_set
      ,v_itr, z_itr
      ,ao
      ,bo
      ,&ext_ro
      ,&roR_scaling
      ,m ? &rom_scaling : NULL
      ,aa
      ,ba
      ,q_hat ? &ext_ra : NULL
      ,q_hat ? &ra_scaling : NULL
      );
    std::copy(ext_ro.begin(),ext_ro.end(),ro.begin());  // Convert back to standard precision
    if(q_hat) std::copy(ext_ra.begin(),ext_ra.end(),ra.begin());
    //
    // Calcuate convergence criteria
    //
    roR_nrm  = norm_inf(ro(1,n_R));
    rom_nrm  = (m ? norm_inf(ro(n_R+1,n_R+m)) : 0.0);
    ra_nrm   = (q_hat ? norm_inf(ra) : 0.0);
    if( iter_refine_k == 0 ) {
      roR_nrm_o = roR_nrm;
      rom_nrm_o = rom_nrm;
      ra_nrm_o  = rom_nrm;
    }
    const bool
      is_roR_conv = roR_nrm / (1.0 + roR_scaling) < iter_refine_opt_tol(),
      is_rom_conv = (m ? rom_nrm / (1.0 + rom_scaling) < iter_refine_feas_tol() : true ),
      is_ra_conv  = (q_hat ?  ra_nrm / (1.0 + ra_scaling) < iter_refine_feas_tol() : true ),
      is_conv     = is_roR_conv && is_rom_conv && is_ra_conv;
    //
    // Print beginning of summary line for residuals
    //
    if( out && (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
      *out
        << right << setw(int_w) << iter_refine_k
        << right << setw(dbl_w) << roR_scaling
        << right << setw(dbl_w) << (roR_nrm / (1.0 + roR_scaling))
        << left  << setw(1)     << (is_roR_conv ? "*" : " ");
      if(m) {
        *out
          << right << setw(dbl_w) << rom_scaling
          << right << setw(dbl_w) << (rom_nrm /(1.0 + rom_scaling))
          << left  << setw(1)     << (is_rom_conv ? "*" : " ");
      }
      if(q_hat) {
        *out
          << right << setw(dbl_w) << ra_scaling
          << right << setw(dbl_w) << (ra_nrm /(1.0 + ra_scaling))
          << left  << setw(1)     << (is_ra_conv ? "*" : " ");
      }
    }
    //
    // Check for convergence
    //
    if( iter_refine_k + 1 < iter_refine_min_iter() ) {
      // Keep on going even if we have converged to desired tolerances!
    }
    else if( is_conv ) {
      if( out && (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
        *out
          << right << setw(dbl_w) << "-"
          << right << setw(dbl_w) << "-"
          << endl;
      }
      if( iter_refine_k == 0 )
        return_status = ITER_REFINE_NOT_NEEDED;
      else
        return_status = ITER_REFINE_CONVERGED;
      break;
    }
    // Make sure we have made progress
    if( roR_nrm_o < roR_nrm && rom_nrm_o < rom_nrm && ra_nrm_o < rom_nrm ) {
      return_status = ITER_REFINE_NOT_IMPROVED;
      break; // No progress was make in converging the equations!
    }
    //
    // Solve for the steps
    //
    // [   Ko     U_hat ] [ del_v ] = [ ro ]
    // [ U_hat'   V_hat ] [ del_z ]   [ ra ]
    //
    ++(*iter_refine_num_solves);
    if( q_hat ) {
      // del_z = inv(S_hat)*(ra - U_hat'*inv(Ko)*ro)
      V_InvMtV( &t1, Ko, no_trans, ro );
      calc_z( act_set.S_hat(), ra, act_set.U_hat(), &t1, &del_z );
    }
    calc_v( Ko, &ro, U_hat, del_z, &del_v );
    //
    // Compute steps:
    //
    // v += -del_v
    // z += -del_z
    //
    Vp_StV( &v_itr, -1.0, del_v );
    if( q_hat )
      Vp_StV( &z_itr, -1.0, del_z );
    //
    // Print rest of summary line for steps
    //
    if( out && (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
      *out
        << right << setw(dbl_w) << norm_inf(del_v) / (norm_inf(v_itr) + small_num)
        << right << setw(dbl_w) << norm_inf(del_z) / (norm_inf(z_itr) + small_num)
        << endl;
    }
  }
  if( iter_refine_max_iter() == 0 ) {
    if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
      *out
        << "\nWarning, iter_refine_max_iter == 0.  Iterative refinement was not performed."
        << "\nLeaving the original solution intact ...\n";
    }
    return_status = ITER_REFINE_NOT_PERFORMED;
  }
  else {
    if( iter_refine_max_iter() == 1 ) {
      if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
        *out
          << "\nWarning, iter_refine_max_iter == 1.  Only one step of iterative refinement"
          << "was performed and the step is taken which out checking the residual ...\n";
      }
      *v = v_itr;
      if(q_hat)
        *z = z_itr;
      return_status = ITER_REFINE_ONE_STEP;
    }
    else if( roR_nrm_o < roR_nrm && rom_nrm_o < rom_nrm && ra_nrm_o < rom_nrm ) {
      if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
        *out
          << "\nNo progress was made in reducing the residuals."
          << "\nLeaving the original solution intact ...\n";
      }
      return_status = ITER_REFINE_NOT_IMPROVED;
    }
    else {
      // The residuals were at least not increased so let's take the new solution
      *v = v_itr;
      if(q_hat)
        *z = z_itr;
      if( return_status != ITER_REFINE_CONVERGED && return_status != ITER_REFINE_NOT_NEEDED ) {
        if( (int)output_level >= (int)OUTPUT_ITER_SUMMARY ) {
          *out
            << "\nThe residuals were not converged but they were not increased either."
            << "\nTake the new point anyway ...\n";
        }
        return_status = ITER_REFINE_IMPROVED;
      }
    }
  }
  return return_status;
}

// private static member functions for QPSchur

void QPSchur::dump_act_set_quantities(
  const ActiveSet& act_set, std::ostream& out
  ,bool print_S_hat
  )
{
  using std::endl;
  using std::setw;
  using std::left;
  using std::right;

  const QPSchurPack::QP
    &qp = act_set.qp();
  const QPSchurPack::Constraints
    &constraints = qp.constraints();

  const int  int_w = 10;
  const char int_ul[] = "--------";
  const int  dbl_min_w = 20;
  const int  dbl_w = my_max(dbl_min_w,int(out.precision()+8));
  const char dbl_ul[] = "------------------";

    out << "\n*** Dumping the current active set ***\n"
    << "\nDimensions of the current active set:\n"
    << "\nn           = " << right << setw(int_w) << qp.n()					<< " (Number of unknowns)"
    << "\nn_R         = " << right << setw(int_w) << qp.n_R()				<< " (Number of initially free variables in Ko)"
    << "\nm           = " << right << setw(int_w) << qp.m()					<< " (Number of initially fixed variables not in Ko)"
    << "\nm_breve     = " << right << setw(int_w) << constraints.m_breve()	<< " (Number of extra general equality/inequality constriants)"
    << "\nq_hat       = " << right << setw(int_w) << act_set.q_hat()		<< " (Number of augmentations to the initial KKT system Ko)"
    << "\nq_plus_hat  = " << right << setw(int_w) << act_set.q_plus_hat()	<< " (Number of added variable bounds and general constraints)"
    << "\nq_F_hat     = " << right << setw(int_w) << act_set.q_F_hat()		<< " (Number of initially fixed variables not at their initial bound)"
    << "\nq_C_hat     = " << right << setw(int_w) << act_set.q_C_hat()		<< " (Number of initially fixed variables at the other bound)"
    << "\nq_D_hat     = " << right << setw(int_w) << act_set.q_D_hat()		<< " (Number of initially fixed variables still fixed at initial bound)"
    << endl;

  // Print table of quantities in augmented KKT system
  out	<< "\nQuantities for augmentations to the initial KKT system:\n";
  const size_type q_hat = act_set.q_hat();
  out	<< endl
    << right << setw(int_w) << "s"
    << right << setw(int_w) << "ij_map(s)"
    << right << setw(int_w) << "bnd(s)"
    << right << setw(dbl_w) << "constr_norm(s)"
    << right << setw(dbl_w) << "d_hat(s)"
    << right << setw(dbl_w) << "z_hat(s)"
    << right << setw(dbl_w) << "p_z_hat(s)"
    << endl;
  out	<< right << setw(int_w) << int_ul
    << right << setw(int_w) << int_ul
    << right << setw(int_w) << int_ul
    << right << setw(dbl_w) << dbl_ul
    << right << setw(dbl_w) << dbl_ul
    << right << setw(dbl_w) << dbl_ul
    << right << setw(dbl_w) << dbl_ul
    << endl;
  {for( size_type s = 1; s <= q_hat; ++s ) {
    out	<< right << setw(int_w) << s
      << right << setw(int_w) << act_set.ij_map(s)
      << right << setw(int_w) << bnd_str(act_set.bnd(s))
      << right << setw(dbl_w) << act_set.constr_norm(s)
      << right << setw(dbl_w) << act_set.d_hat()(s)
      << right << setw(dbl_w) << act_set.z_hat()(s)
      << right << setw(dbl_w) << act_set.p_z_hat()(s)
      << endl;
  }}
  
  // Print P_XF_hat, P_FC_hat, P_plus_hat, U_hat and S_hat
  out	<< "\nP_XF_hat =\n" 	<< act_set.P_XF_hat();
  out	<< "\nP_FC_hat =\n" 	<< act_set.P_FC_hat();
  out	<< "\nP_plus_hat =\n" 	<< act_set.P_plus_hat();
  out	<< "\nU_hat =\n" 		<< act_set.U_hat();
  if(print_S_hat)
    out	<< "\nS_hat =\n" 		<< act_set.S_hat();
  
  // Print table of multipliers for q_D_hat
  out	<< "\nQuantities for initially fixed variables which are still fixed at their initial bound:\n";
  const size_type q_D_hat = act_set.q_D_hat();
  out	<< endl
    << right << setw(int_w) << "k"
    << right << setw(int_w) << "l_fxfx(k)"
    << right << setw(dbl_w) << "mu_D_hat(k)"
    << right << setw(dbl_w) << "p_mu_D_hat(s)"
    << endl;
  out	<< right << setw(int_w) << int_ul
    << right << setw(int_w) << int_ul
    << right << setw(dbl_w) << dbl_ul
    << right << setw(dbl_w) << dbl_ul
    << endl;
  {for( size_type k = 1; k <= q_D_hat; ++k ) {
    out	<< right << setw(int_w) << k
      << right << setw(int_w) << act_set.l_fxfx(k)
      << right << setw(dbl_w) << act_set.mu_D_hat()(k)
      << right << setw(dbl_w) << act_set.p_mu_D_hat()(k)
      << endl;
  }}
  
  // Print Q_XD_hat
  out	<< "\nQ_XD_hat =\n" << act_set.Q_XD_hat();

  out << "\n*** End dump of current active set ***\n";
}

// QPSchurPack::QP

void QPSchurPack::QP::dump_qp( std::ostream& out )
{
  using std::endl;
  using std::setw;
  using std::left;
  using std::right;
  using Teuchos::Workspace;
  Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

  const Constraints
    &constraints = this->constraints();

  const size_type
    n = this->n(),
    n_R = this->n_R(),
    m = this->m(),
    m_breve = constraints.m_breve();

  out	<< "\n*** Original QP ***\n"
    << "\nn       = " << n
    << "\nm       = " << m
    << "\nm_breve = " << m_breve
    << endl;
  out	<< "\ng =\n" << g();
  out	<< "\nG =\n" << G();
  if(m) {
    out	<< "\nA =\n" << A();
    // Le'ts recover c from fo(n_R+1:n_R+m) = c - A' * Q_X * b_x
    throw std::logic_error(
      error_msg(__FILE__,__LINE__,"QPSchurPack::QP::dump_qp(...) : Error, "
            "m != not supported yet!"));  
    // ToDo: Implement this when needed!
  }
  out	<< "\nA_bar =\n" << constraints.A_bar();
  // Get c_L_bar and c_U_bar
  DVector c_L_bar(n+m_breve), c_U_bar(n+m_breve);
  {for( size_type j = 1; j <= n+m_breve; ++j ){
    c_L_bar(j) = constraints.get_bnd(j,LOWER);
    c_U_bar(j) = constraints.get_bnd(j,UPPER);
  }}
  out	<< "\nc_L_bar =\n" << c_L_bar();
  out	<< "\nc_U_bar =\n" << c_U_bar();
  
  out	<< "\n*** Initial KKT system (fixed and free variables) ***\n"
    << "\nn_R = " << n_R
    << endl;
  out	<< "\nb_X =\n" << b_X();
  out	<< "\nQ_R =\n" << Q_R();
  out	<< "\nQ_X =\n" << Q_X();
  out	<< "\nKo =\n" << Ko();
  out	<< "\nfo =\n" << fo();
}

}	// end namespace ConstrainedOptPack
