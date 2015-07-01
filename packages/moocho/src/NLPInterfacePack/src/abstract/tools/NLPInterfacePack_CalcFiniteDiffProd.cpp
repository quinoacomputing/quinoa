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

#include <assert.h>
#include <math.h>

#include <typeinfo>
#include <iomanip>
#include <sstream>
#include <limits>

#include "NLPInterfacePack_CalcFiniteDiffProd.hpp"
#include "NLPInterfacePack_NLP.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_assert_print_nan_inf.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Assert.hpp"

namespace NLPInterfacePack {

CalcFiniteDiffProd::CalcFiniteDiffProd(
  EFDMethodOrder              fd_method_order
  ,EFDStepSelect              fd_step_select
  ,value_type                 fd_step_size
  ,value_type                 fd_step_size_min
  ,value_type                 fd_step_size_f
  ,value_type                 fd_step_size_c
  )
  :fd_method_order_(fd_method_order)
  ,fd_step_select_(fd_step_select)
  ,fd_step_size_(fd_step_size)
  ,fd_step_size_min_(fd_step_size_min)
  ,fd_step_size_f_(fd_step_size_f)
  ,fd_step_size_c_(fd_step_size_c)
{}

bool CalcFiniteDiffProd::calc_deriv_product(
  const Vector       &xo
  ,const Vector      *xl
  ,const Vector      *xu
  ,const Vector      &v
  ,const value_type  *fo
  ,const Vector      *co
  ,bool              check_nan_inf
  ,NLP               *nlp
  ,value_type        *Gf_prod
  ,VectorMutable     *Gc_prod
  ,std::ostream      *out_arg
  ,bool              trace
  ,bool              dump_all
  ) const
{

  using std::setw;
  using std::endl;
  using std::right;

  using BLAS_Cpp::rows;
  using BLAS_Cpp::cols;
  
  typedef VectorSpace::vec_mut_ptr_t  vec_mut_ptr_t;
  using AbstractLinAlgPack::Vt_S;
  using AbstractLinAlgPack::Vp_StV;
  using AbstractLinAlgPack::max_near_feas_step;
  using AbstractLinAlgPack::assert_print_nan_inf;
  using LinAlgOpPack::V_StV;

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::getFancyOStream(Teuchos::rcp(out_arg,false));
  Teuchos::OSTab tab(out);

  //
  // The gradient of the contraints is defined as the matrix Gc as:
  //
  // Gc= [ Gc1, Gc2, ..., Gcm ]
  //
  //     [  dc1/dx(1)    dc2/dx(1)  ...    dcm/dx(1)  ]
  //     [  dc1/dx(2)    dc2/dx(2)  ...    dcm/dx(2)  ]
  // Gc= [  .            .          ...    .          ]
  //     [  dc1/dx(n)    dc2/dx(n)  ...    dcm/dx(n)  ]
  //
  //     [  (dc/dx(1))'  ]
  //     [  (dc/dx(2))'  ]
  // Gc= [  .            ]
  //     [  (dc/dx(n))'  ]
  //
  // The gradient of the objective function is defined as the
  // vector Gf as:
  //
  //     [  (df/dx(1))'  ]
  //     [  (df/dx(2))'  ]
  // Gf= [  .            ]
  //     [  (df/dx(n))'  ]
  //
  // To illustrate the theory behind this implementation consider
  // the generic multi-variable function g(x) <: R^n -> R.  Now let's
  // consider we have the base point xo and the vector v to
  // perturb g(x) along.  First form the function g(xo+a*v) and then
  // let's compute dg/da at a = 0:
  // 
  // (1) d(g(xo+a*v))/d(a) at a = 0
  //         = sum( dg/dx(i) * dx(i)/da, i = 1...n)
  //         = sum( dg/dx(i) * v(i), i = 1...n)
  //         = Gf'*v
  //
  // Now we can approximate (1) using central differences as:
  // 
  // (2) d(g(xo+a*v))/d(a) at a = 0
  //          \approx ( g(xo+h*v) - g(xo+h*v) ) / (2*h)
  //
  // If we equate (1) and (2) we have the approximation:
  // 
  // (3) Gg' * v \approx ( g(xo+h*v) - g(xo+h*v) ) / (2*h)
  // 
  // It should be clear how this applies to computing Gf'*v and Gc'*v.
  // 

  const size_type
    n = nlp->n(),
    m = nlp->m();

  const value_type
    max_bnd_viol = nlp->max_var_bounds_viol();

  // /////////////////////////////////////////
  // Validate the input

  TEUCHOS_TEST_FOR_EXCEPTION(
    m==0 && Gc_prod, std::invalid_argument
    ,"CalcFiniteDiffProd::calc_deriv(...) : "
    "Error, if nlp->m() == 0, then Gc_prod must equal NULL" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    Gc_prod && !Gc_prod->space().is_compatible(*nlp->space_c())
    ,std::invalid_argument
    ,"CalcFiniteDiffProd::calc_deriv(...) : "
    "Error, Gc_prod (type \' "<<typeName(*Gc_prod)<<"\' "
    "is not compatible with the NLP" );
  TEUCHOS_TEST_FOR_EXCEPTION(
    (xl && !xu) || (!xl && xu), std::invalid_argument
    ,"CalcFiniteDiffProd::calc_deriv(...) : "
    "Error, both xl = "<<xl<<" and xu = "<<xu
    <<" must be NULL or not NULL" );

  assert_print_nan_inf(xo,"xo",true,out.get()); 

  switch(this->fd_method_order()) {
    case FD_ORDER_ONE:
      if(out.get()&&trace) *out<<"\nUsing one-sided, first-order finite differences ...\n";
      break;
    case FD_ORDER_TWO:
      if(out.get()&&trace) *out<<"\nUsing one-sided, second-order finite differences ...\n";
      break;
    case FD_ORDER_TWO_CENTRAL:
      if(out.get()&&trace) *out<<"\nUsing second-order central finite differences ...\n";
      break;
    case FD_ORDER_TWO_AUTO:
      if(out.get()&&trace) *out<<"\nUsing auto selection of some second-order finite difference method ...\n";
      break;
    case FD_ORDER_FOUR:
      if(out.get()&&trace) *out<<"\nUsing one-sided, fourth-order finite differences ...\n";
      break;
    case FD_ORDER_FOUR_CENTRAL:
      if(out.get()&&trace) *out<<"\nUsing fourth-order central finite differences ...\n";
      break;
    case FD_ORDER_FOUR_AUTO:
      if(out.get()&&trace) *out<<"\nUsing auto selection of some fourth-order finite difference method ...\n";
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here!
  }

  // ////////////////////////
  // Find the step size

  //
  // Get defaults for the optimal step sizes
  //

  const value_type
    sqrt_epsilon = ::pow(std::numeric_limits<value_type>::epsilon(),1.0/2.0),
    u_optimal_1  = sqrt_epsilon,
    u_optimal_2  = ::pow(sqrt_epsilon,1.0/2.0),
    u_optimal_4  = ::pow(sqrt_epsilon,1.0/4.0),
    xo_norm_inf  = xo.norm_inf();

  value_type
    uh_opt = 0.0;
  switch(this->fd_method_order()) {
    case FD_ORDER_ONE:
      uh_opt = u_optimal_1 * ( fd_step_select() == FD_STEP_ABSOLUTE ? 1.0 : xo_norm_inf + 1.0 );
      break;
    case FD_ORDER_TWO:
    case FD_ORDER_TWO_CENTRAL:
    case FD_ORDER_TWO_AUTO:
      uh_opt = u_optimal_2 * ( fd_step_select() == FD_STEP_ABSOLUTE ? 1.0 : xo_norm_inf + 1.0 );
      break;
    case FD_ORDER_FOUR:
    case FD_ORDER_FOUR_CENTRAL:
    case FD_ORDER_FOUR_AUTO:
      uh_opt = u_optimal_4 * ( fd_step_select() == FD_STEP_ABSOLUTE ? 1.0 : xo_norm_inf + 1.0 );
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here!
  }

  if(out.get()&&trace) *out<<"\nDefault optimal step length uh_opt = " << uh_opt << " ...\n";

  //
  // Set the step sizes used.
  //

  value_type
    uh      = this->fd_step_size(),
    uh_f    = this->fd_step_size_f(),
    uh_c    = this->fd_step_size_c(),
    uh_min  = this->fd_step_size_min();

  // uh
  if( uh < 0 )
    uh = uh_opt;
  else if(fd_step_select() == FD_STEP_RELATIVE)
    uh *= (xo_norm_inf + 1.0);
  // uh_f
  if( uh_f < 0 )
    uh_f = uh;
  else if(fd_step_select() == FD_STEP_RELATIVE)
    uh_f *= (xo_norm_inf + 1.0);
  // uh_c
  if( uh_c < 0 )
    uh_c = uh;
  else if(fd_step_select() == FD_STEP_RELATIVE)
    uh_c *= (xo_norm_inf + 1.0);

  if(out.get()&&trace) *out<<"\nIndividual step sizes initally set: uh="<<uh<<",uh_f="<<uh_f<<",uh_c="<<uh_c<<"\n";

  //
   // Determine the maximum step size that can be used and
  // still stay in the relaxed bounds.
  //
  // ToDo: Consider cramped bounds, one sided differences!
  //

  value_type max_u_feas = std::numeric_limits<value_type>::max();
  if( xl ) {
    std::pair<value_type,value_type>
      u_pn
      = max_near_feas_step(
        xo
        ,v
        ,*xl
        ,*xu
        ,max_bnd_viol
        );
    if( u_pn.first < -u_pn.second )
      max_u_feas = u_pn.first;
    else
      max_u_feas = u_pn.second;
    const value_type abs_max_u_feas = ::fabs(max_u_feas);
    if( abs_max_u_feas < uh ) {
      if( abs_max_u_feas < uh_min ) {
        if(out.get())
          *out
            << "Warning, the size of the maximum finite difference step length\n"
            << "that does not violate the relaxed variable bounds uh = "
            << max_u_feas << " is less than the mimimum allowable step length\n"
            << "uh_min = " << uh_min << " and the finite difference "
            << "derivatives are not computed!\n";
        return false;
      }
      if(out.get())
        *out
          << "Warning, the size of the maximum finite difference step length\n"
          << "that does not violate the relaxed variable bounds uh = "
          << max_u_feas << " is less than the desired step length\n"
          << "uh = " << uh << " and the finite difference "
          << "derivatives may be much less accurate!\n";
    }
  }

  //
  // Set the actual method being used
  //
  // ToDo: Consider cramped bounds and method order.
  //
  
  EFDMethodOrder  fd_method_order = this->fd_method_order();
  switch(fd_method_order) {
    case FD_ORDER_TWO_AUTO:
      fd_method_order = FD_ORDER_TWO_CENTRAL;
      break;
    case FD_ORDER_FOUR_AUTO:
      fd_method_order = FD_ORDER_FOUR_CENTRAL;
      break;
  }

  // Compute the actual individual step size so as to stay in bounds
  const value_type
    abs_max_u_feas = ::fabs(max_u_feas);
  value_type
     num_u_i = 0;
  switch(fd_method_order) {
    case FD_ORDER_ONE:
      num_u_i = 1.0;
      break;
    case FD_ORDER_TWO:
      num_u_i = 2.0;
      break;
    case FD_ORDER_TWO_CENTRAL:
      num_u_i = 1.0;
      break;
    case FD_ORDER_FOUR:
      num_u_i = 4.0;
      break;
    case FD_ORDER_FOUR_CENTRAL:
      num_u_i = 2.0;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here!
  }

  uh   = ( abs_max_u_feas/num_u_i < uh   ? max_u_feas/num_u_i : uh   ); // This can be a negative number!
  uh_f = ( abs_max_u_feas/num_u_i < uh_f ? max_u_feas/num_u_i : uh_f ); //""
  uh_c = ( abs_max_u_feas/num_u_i < uh_c ? max_u_feas/num_u_i : uh_c ); //""

  if( uh_min < 0 ) {
    uh_min = uh / 100.0;
  }

  if(out.get()&&trace) *out<<"\nIndividual step sizes to fit in bounds: uh="<<uh<<",uh_f="<<uh_f<<",uh_c="<<uh_c<<"\n";

  //
  // Remember some stuff
  //
  
  value_type        *f_saved = NULL;
  VectorMutable     *c_saved = NULL;

  f_saved = nlp->get_f();
  if(m)  c_saved = nlp->get_c();

  int p_saved;
  if(out.get())
    p_saved = out->precision();

  // ///////////////////////////////////////////////
  // Compute the directional derivatives

  try {

  value_type
    f;
  vec_mut_ptr_t
    x = nlp->space_x()->create_member();
  vec_mut_ptr_t
    c = m  && Gc_prod ? nlp->space_c()->create_member() : Teuchos::null;
  
  // Set the quanitities used to compute with

  nlp->set_f(&f);
  if(m)  nlp->set_c( c.get() );

  const int dbl_p = 15;
  if(out.get())
    *out << std::setprecision(dbl_p);

  //
  // Compute the weighted sum of the terms
  //

  int          num_evals  = 0;
  value_type   dwgt       = 0.0;
  switch(fd_method_order) {
    case FD_ORDER_ONE: // may only need one eval if f(xo) etc is passed in
      num_evals = 2;
      dwgt      = 1.0;
      break;
    case FD_ORDER_TWO: // may only need two evals if c(xo) etc is passed in
      num_evals = 3;
      dwgt      = 2.0;
      break;
    case FD_ORDER_TWO_CENTRAL:
      num_evals = 2;
      dwgt      = 2.0;
      break;
    case FD_ORDER_FOUR:
      num_evals = 5;
      dwgt      = 12.0;
      break;
    case FD_ORDER_FOUR_CENTRAL:
      num_evals = 5;
      dwgt      = 12.0;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true); // Should not get here!
  }
  if(Gc_prod) *Gc_prod = 0.0;
  if(Gf_prod) *Gf_prod = 0.0;
  for( int eval_i = 1; eval_i <= num_evals; ++eval_i ) {
    // Set the step constant and the weighting constant
    value_type
      uh_i   = 0.0,
      wgt_i  = 0.0;
    switch(fd_method_order) {
      case FD_ORDER_ONE: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -1.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +1.0;
            break;
        }
        break;
      }
      case FD_ORDER_TWO: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -3.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +4.0;
            break;
          case 3:
            uh_i  = +2.0;
            wgt_i = -1.0;
            break;
        }
        break;
      }
      case FD_ORDER_TWO_CENTRAL: {
        switch(eval_i) {
          case 1:
            uh_i  = -1.0;
            wgt_i = -1.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +1.0;
            break;
        }
        break;
      }
      case FD_ORDER_FOUR: {
        switch(eval_i) {
          case 1:
            uh_i  = +0.0;
            wgt_i = -25.0;
            break;
          case 2:
            uh_i  = +1.0;
            wgt_i = +48.0;
            break;
          case 3:
            uh_i  = +2.0;
            wgt_i = -36.0;
            break;
          case 4:
            uh_i  = +3.0;
            wgt_i = +16.0;
            break;
          case 5:
            uh_i  = +4.0;
            wgt_i = -3.0;
            break;
        }
        break;
      }
      case FD_ORDER_FOUR_CENTRAL: {
        switch(eval_i) {
          case 1:
            uh_i  = -2.0;
            wgt_i = +1.0;
            break;
          case 2:
            uh_i  = -1.0;
            wgt_i = -8.0;
            break;
          case 3:
            uh_i  = +1.0;
            wgt_i = +8.0;
            break;
          case 4:
            uh_i  = +2.0;
            wgt_i = -1.0;
            break;
        }
        break;
      }
    }

    if(out.get()&&dump_all) {
      *out<<"\nxo =\n" << xo;
      *out<<"\nv =\n" << v;
      if(fo) *out << "\nfo = " << *fo << "\n";
      if(co) *out << "\nco =\n" << *co;
    }
    // Compute the weighted term and add it to the sum
    bool new_point = true;
    if(Gc_prod) {
      if( co && uh_i == 0.0 ) {
        if(out.get()&&trace) *out<<"\nBase c = co ...\n";
        *c = *co;
      }
      else {
        if( new_point || uh_c != uh ) {
          *x = xo; Vp_StV( x.get(), uh_i * uh_c, v ); // x = xo + uh_i*uh_c*v
        }
        if(out.get()&&trace) *out<<"\nComputing c = c(xo+"<<(uh_i*uh_c)<<"*v) ...\n";
        if(out.get()&&dump_all) *out<<"\nxo+"<<(uh_i*uh_c)<<"*v =\n" << *x;
        nlp->calc_c(*x,new_point);
      }
      new_point = false;
      if(out.get() && dump_all) *out << "\nc =\n" << *c;
      if(check_nan_inf)
        assert_print_nan_inf(*c,"c(xo+u*v)",true,out.get());
      if(out.get()&&trace) *out<<"\nGc_prod += "<<wgt_i<<"*c ...\n";
      Vp_StV( Gc_prod, wgt_i, *c );
      if(out.get() && dump_all) *out<<"\nGc_prod =\n" << *Gc_prod;
    }
    
    if(Gf_prod) {
      if( fo && uh_i == 0.0 ) {
        if(out.get()&&trace) *out<<"\nBase f = fo ...\n";
        f = *fo;
      }
      else {
        if( new_point || uh_f != uh ) {
          *x = xo; Vp_StV( x.get(), uh_i * uh_f, v ); // x = xo + uh_i*uh_f*v
          new_point = true;
        }
        if(out.get()&&trace) *out<<"\nComputing f = f(xo+"<<(uh_i*uh_f)<<"*v) ...\n";
        if(out.get() && dump_all) *out<<"\nxo+"<<(uh_i*uh_f)<<"*v =\n" << *x;
        nlp->calc_f(*x,new_point);
      }
      new_point = false;
      if(out.get() && dump_all) *out<<"\nf = " << f << "\n";
      if(check_nan_inf)
        assert_print_nan_inf(f,"f(xo+u*v)",true,out.get());
      if(out.get()&&trace) *out<<"\nGf_prod += "<<wgt_i<<"*f ...\n";
      *Gf_prod += wgt_i * f;
      if(out.get() && dump_all) *out<<"\nGf_prod = " << *Gf_prod << "\n";
    }

  }

  //
  // Multiply by the scaling factor!
  //

  if(Gc_prod) {
    if(out.get()&&trace) *out<<"\nGc_prod *= "<<(1.0 / (dwgt * uh_c))<<" ...\n";
    Vt_S( Gc_prod, 1.0 / (dwgt * uh_c) );
    if(out.get() && dump_all)
      *out<<"\nFinal Gc_prod =\n" << *Gc_prod;
  }
    
  if(Gf_prod) {
    if(out.get()&&trace) *out<<"\nGf_prod *= "<<(1.0 / (dwgt * uh_f))<<" ...\n";
    *Gf_prod *= ( 1.0 / (dwgt * uh_f) );
    if(out.get() && dump_all)
      *out<<"\nFinal Gf_prod = " << *Gf_prod << "\n";
  }

  }  // end try
  catch(...) {
    nlp->set_f( f_saved );
    if(m)  nlp->set_c( c_saved );
    if(out.get())
      *out << std::setprecision(p_saved);
    throw;
  }
  
  nlp->set_f( f_saved );
  if(m)  nlp->set_c( c_saved );
  if(out.get())
    *out << std::setprecision(p_saved);
  
  return true;
}

}  // end namespace NLPInterfacePack
