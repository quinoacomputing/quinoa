// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.  
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  A short implementation ( not all operators and 
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
// @HEADER

#include "Sacado_ConfigDefs.h"
#include "Sacado_mpl_range_c.hpp"
#include "Sacado_mpl_for_each.hpp"

// template <typename T, typename Storage> 
// template <typename S> 
// inline Sacado::ELRFad::GeneralFad<T,Storage>::GeneralFad(const Expr<S>& x) :
//   Storage(T(0.))
// {
//   int sz = x.size();

//   if (sz != this->size()) 
//     this->resize(sz);

//   if (sz) {

//     // Number of arguments
//     const int N = Expr<S>::num_args;

//     // Array to store partials
//     T partials[N];

//     // Array to store tangents
//     T dots[N];

//     // Compute partials
//     x.computePartials(T(1.), partials);
  
//     // Compute each tangent direction
//     for(int i=0; i<sz; ++i) {
//       T t = T(0.);
//       x.getTangents(i,dots);
//       for (int j=0; j<N; j++) {
// 	t += partials[j]*dots[j];
//       }
//       this->fastAccessDx(i) = t;
//     }

//   }
  
//   // Compute value
//   this->val() = x.val();
// }

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>::GeneralFad(const Expr<S>& x) :
  Storage(T(0.)),
  update_val_(x.updateValue())
{
  int sz = x.size();

  if (sz != this->size()) 
    this->resize(sz);

  if (sz) {

    if (Expr<S>::is_linear) {
      if (x.hasFastAccess())
        for(int i=0; i<sz; ++i)
          this->fastAccessDx(i) = x.fastAccessDx(i);
      else
        for(int i=0; i<sz; ++i)
          this->fastAccessDx(i) = x.dx(i);
    }
    else {

      // Number of arguments
      const int N = Expr<S>::num_args;

      if (x.hasFastAccess()) {
	// Compute partials
	FastLocalAccumOp< Expr<S> > op(x);
  
	// Compute each tangent direction
	for(op.i=0; op.i<sz; ++op.i) {
	  op.t = T(0.);

	  // Automatically unrolled loop that computes
	  // for (int j=0; j<N; j++)
	  //   op.t += op.partials[j] * x.getTangent<j>(i);
	  Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	  this->fastAccessDx(op.i) = op.t;
	}
      }
      else {
	// Compute partials
	SlowLocalAccumOp< Expr<S> > op(x);
  
	// Compute each tangent direction
	for(op.i=0; op.i<sz; ++op.i) {
	  op.t = T(0.);

	  // Automatically unrolled loop that computes
	  // for (int j=0; j<N; j++)
	  //   op.t += op.partials[j] * x.getTangent<j>(i);
	  Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	  this->fastAccessDx(op.i) = op.t;
	}
      }
 
    }

  }
  
  // Compute value
  if (update_val_)
    this->val() = x.val();
}

template <typename T, typename Storage> 
inline void 
Sacado::ELRFad::GeneralFad<T,Storage>::diff(const int ith, const int n) 
{ 
  if (this->size() != n) 
    this->resize(n);

  this->zero();
  this->fastAccessDx(ith) = T(1.);

}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator=(const T& v) 
{
  this->val() = v;

  if (this->size()) {
    this->zero();    // We need to zero out the array for future resizes
    this->resize(0);
  }

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator=(
			       const Sacado::ELRFad::GeneralFad<T,Storage>& x) 
{
  // Copy value & dx_
  Storage::operator=(x);
  update_val_ = x.update_val_;
  
  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator=(const Expr<S>& x) 
{
  int sz = x.size();

  if (sz != this->size()) 
    this->resize(sz);

  if (sz) {
 
    if (Expr<S>::is_linear) {
      if (x.hasFastAccess())
        for(int i=0; i<sz; ++i)
          this->fastAccessDx(i) = x.fastAccessDx(i);
      else
        for(int i=0; i<sz; ++i)
          this->fastAccessDx(i) = x.dx(i);
    }
    else {

      // Number of arguments
      const int N = Expr<S>::num_args;

      if (x.hasFastAccess()) {
	// Compute partials
	FastLocalAccumOp< Expr<S> > op(x);
  
	// Compute each tangent direction
	for(op.i=0; op.i<sz; ++op.i) {
	  op.t = T(0.);

	  // Automatically unrolled loop that computes
	  // for (int j=0; j<N; j++)
	  //   op.t += op.partials[j] * x.getTangent<j>(i);
	  Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	  this->fastAccessDx(op.i) = op.t;
	}
      }
      else {
	// Compute partials
	SlowLocalAccumOp< Expr<S> > op(x);
  
	// Compute each tangent direction
	for(op.i=0; op.i<sz; ++op.i) {
	  op.t = T(0.);

	  // Automatically unrolled loop that computes
	  // for (int j=0; j<N; j++)
	  //   op.t += op.partials[j] * x.getTangent<j>(i);
	  Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	  this->fastAccessDx(op.i) = op.t;
	}
      }
    }
  }
  
  // Compute value
  update_val_ = x.updateValue();
  if (update_val_)
    this->val() = x.val();
  
  return *this;
}

template <typename T, typename Storage> 
inline  Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator += (const T& v)
{
  if (update_val_)
    this->val() += v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator -= (const T& v)
{
  if (update_val_)
    this->val() -= v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator *= (const T& v)
{
  int sz = this->size();

  if (update_val_)
    this->val() *= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) *= v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator /= (const T& v)
{
  int sz = this->size();

  if (update_val_)
    this->val() /= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) /= v;

  return *this;
}

template <typename T, typename Storage> 
inline  Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::
operator += (const typename Sacado::dummy<value_type,scalar_type>::type& v)
{
  if (update_val_)
    this->val() += v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::
operator -= (const typename Sacado::dummy<value_type,scalar_type>::type& v)
{
  if (update_val_)
    this->val() -= v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::
operator *= (const typename Sacado::dummy<value_type,scalar_type>::type& v)
{
  int sz = this->size();

  if (update_val_)
    this->val() *= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) *= v;

  return *this;
}

template <typename T, typename Storage> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::
operator /= (const typename Sacado::dummy<value_type,scalar_type>::type& v)
{
  int sz = this->size();

  if (update_val_)
    this->val() /= v;
  for (int i=0; i<sz; ++i)
    this->fastAccessDx(i) /= v;

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator += (
					     const Sacado::ELRFad::Expr<S>& x)
{
  int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (Expr<S>::is_linear) {
    if (xsz) {
      if (sz) {
	if (x.hasFastAccess())
	  for (int i=0; i<sz; ++i)
	    this->fastAccessDx(i) += x.fastAccessDx(i);
	else
	  for (int i=0; i<sz; ++i)
	    this->fastAccessDx(i) += x.dx(i);
      }
      else {
	this->resize(xsz);
	if (x.hasFastAccess())
	  for (int i=0; i<xsz; ++i)
	    this->fastAccessDx(i) = x.fastAccessDx(i);
	else
	  for (int i=0; i<xsz; ++i)
	    this->fastAccessDx(i) = x.dx(i);
      }
    }
  }
  else {

    // Number of arguments
    const int N = Expr<S>::num_args;
    
    if (xsz) {

      if (sz != xsz) {
	this->resize(xsz);
	this->zero();
      }

      if (x.hasFastAccess()) {
	// Compute partials
	FastLocalAccumOp< Expr<S> > op(x);
	
	// Compute each tangent direction
	for(op.i=0; op.i<xsz; ++op.i) {
	  op.t = T(0.);
	  
	  // Automatically unrolled loop that computes
	  // for (int j=0; j<N; j++)
	  //   op.t += op.partials[j] * x.getTangent<j>(i);
	  Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	  this->fastAccessDx(op.i) += op.t;
	}
      }
      else {
	// Compute partials
	SlowLocalAccumOp< Expr<S> > op(x);
  
	// Compute each tangent direction
	for(op.i=0; op.i<xsz; ++op.i) {
	  op.t = T(0.);

	  // Automatically unrolled loop that computes
	  // for (int j=0; j<N; j++)
	  //   op.t += op.partials[j] * x.getTangent<j>(i);
	  Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	  this->fastAccessDx(op.i) += op.t;
	}
      }
      
    }

  }
  
  // Compute value
  update_val_ = x.updateValue();
  if (update_val_)
    this->val() += x.val();

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator -= (
					      const Sacado::ELRFad::Expr<S>& x)
{
  int xsz = x.size(), sz = this->size();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (Expr<S>::is_linear) {
    if (xsz) {
      if (sz) {
	if (x.hasFastAccess())
	  for(int i=0; i<sz; ++i)
	    this->fastAccessDx(i) -= x.fastAccessDx(i);
	else
	  for (int i=0; i<sz; ++i)
	    this->fastAccessDx(i) -= x.dx(i);
      }
      else {
	this->resize(xsz);
	if (x.hasFastAccess())
	  for(int i=0; i<xsz; ++i)
	    this->fastAccessDx(i) = -x.fastAccessDx(i);
	else
	  for (int i=0; i<xsz; ++i)
	    this->fastAccessDx(i) = -x.dx(i);
      }
    }
  }
  else {
    // Number of arguments
    const int N = Expr<S>::num_args;
    
    if (xsz) {
      
      if (sz != xsz) {
	this->resize(xsz);
	this->zero();
      }

      if (x.hasFastAccess()) {
	// Compute partials
	FastLocalAccumOp< Expr<S> > op(x);
	
	// Compute each tangent direction
	for(op.i=0; op.i<xsz; ++op.i) {
	  op.t = T(0.);
	  
	  // Automatically unrolled loop that computes
	  // for (int j=0; j<N; j++)
	  //   op.t += op.partials[j] * x.getTangent<j>(i);
	  Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	  this->fastAccessDx(op.i) -= op.t;
	}
      }
      else {
	// Compute partials
	SlowLocalAccumOp< Expr<S> > op(x);
  
	// Compute each tangent direction
	for(op.i=0; op.i<xsz; ++op.i) {
	  op.t = T(0.);

	  // Automatically unrolled loop that computes
	  // for (int j=0; j<N; j++)
	  //   op.t += op.partials[j] * x.getTangent<j>(i);
	  Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	  this->fastAccessDx(op.i) -= op.t;
	}	
      }
    }

  }

  update_val_ = x.updateValue();
  if (update_val_)
    this->val() -= x.val();

  return *this;
}

template <typename T, typename Storage> 
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator *= (
					     const Sacado::ELRFad::Expr<S>& x)
{
  int xsz = x.size(), sz = this->size();
  T xval = x.val();
  T v = this->val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (Expr<S>::is_linear) {
    if (xsz) {
      if (sz) {
	if (x.hasFastAccess())
	  for(int i=0; i<sz; ++i)
	    this->fastAccessDx(i) = v*x.fastAccessDx(i) + this->fastAccessDx(i)*xval;
	else
	  for (int i=0; i<sz; ++i)
	    this->fastAccessDx(i) = v*x.dx(i) + this->fastAccessDx(i)*xval;
      }
      else {
	this->resize(xsz);
	if (x.hasFastAccess())
	  for(int i=0; i<xsz; ++i)
	    this->fastAccessDx(i) = v*x.fastAccessDx(i);
	else
	  for (int i=0; i<xsz; ++i)
	    this->fastAccessDx(i) = v*x.dx(i);
      }
    }
    else {
      if (sz) {
	for (int i=0; i<sz; ++i)
	  this->fastAccessDx(i) *= xval;
      }
    }
  }
  else {

    // Number of arguments
    const int N = Expr<S>::num_args;
    
    if (xsz) {
      
      if (sz) {

	if (x.hasFastAccess()) {
	  // Compute partials
	  FastLocalAccumOp< Expr<S> > op(x);
	  
	  // Compute each tangent direction
	  for(op.i=0; op.i<xsz; ++op.i) {
	    op.t = T(0.);
	    
	    // Automatically unrolled loop that computes
	    // for (int j=0; j<N; j++)
	    //   op.t += op.partials[j] * x.getTangent<j>(i);
	    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	    this->fastAccessDx(op.i) = 
	      v * op.t + this->fastAccessDx(op.i) * xval;
	  }
	}
	else {
	  // Compute partials
	  SlowLocalAccumOp< Expr<S> > op(x);
	  
	  // Compute each tangent direction
	  for(op.i=0; op.i<xsz; ++op.i) {
	    op.t = T(0.);
	    
	    // Automatically unrolled loop that computes
	    // for (int j=0; j<N; j++)
	    //   op.t += op.partials[j] * x.getTangent<j>(i);
	    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	    this->fastAccessDx(op.i) = 
	      v * op.t + this->fastAccessDx(op.i) * xval;
	  }
	}

      }
      
      else {
	
	this->resize(xsz);
	
	if (x.hasFastAccess()) {
	  // Compute partials
	  FastLocalAccumOp< Expr<S> > op(x);
	  
	  // Compute each tangent direction
	  for(op.i=0; op.i<xsz; ++op.i) {
	    op.t = T(0.);
	    
	    // Automatically unrolled loop that computes
	    // for (int j=0; j<N; j++)
	    //   op.t += op.partials[j] * x.getTangent<j>(i);
	    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	    this->fastAccessDx(op.i) = v * op.t;
	  }
	}
	else {
	  // Compute partials
	  SlowLocalAccumOp< Expr<S> > op(x);
	  
	  // Compute each tangent direction
	  for(op.i=0; op.i<xsz; ++op.i) {
	    op.t = T(0.);
	    
	    // Automatically unrolled loop that computes
	    // for (int j=0; j<N; j++)
	    //   op.t += op.partials[j] * x.getTangent<j>(i);
	    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	    this->fastAccessDx(op.i) = v * op.t;
	  }
	}
	
      }
      
    }
    
    else {
      
      if (sz) {
	for (int i=0; i<sz; ++i)
	  this->fastAccessDx(i) *= xval;
      }
      
    }

  }

  update_val_ = x.updateValue();
  if (update_val_)
    this->val() *= xval;

  return *this;
}

template <typename T, typename Storage>
template <typename S> 
inline Sacado::ELRFad::GeneralFad<T,Storage>& 
Sacado::ELRFad::GeneralFad<T,Storage>::operator /= (
					     const Sacado::ELRFad::Expr<S>& x)
{
  int xsz = x.size(), sz = this->size();
  T xval = x.val();
  T v = this->val();

#ifdef SACADO_DEBUG
  if ((xsz != sz) && (xsz != 0) && (sz != 0))
    throw "Fad Error:  Attempt to assign with incompatible sizes";
#endif

  if (Expr<S>::is_linear) {
    if (xsz) {
      if (sz) {
	if (x.hasFastAccess())
	  for(int i=0; i<sz; ++i)
	    this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - v*x.fastAccessDx(i) )/ (xval*xval);
	else
	  for (int i=0; i<sz; ++i)
	    this->fastAccessDx(i) = ( this->fastAccessDx(i)*xval - v*x.dx(i) )/ (xval*xval);
      }
      else {
	this->resize(xsz);
	if (x.hasFastAccess())
	  for(int i=0; i<xsz; ++i)
	    this->fastAccessDx(i) = - v*x.fastAccessDx(i) / (xval*xval);
	else
	  for (int i=0; i<xsz; ++i)
	    this->fastAccessDx(i) = -v*x.dx(i) / (xval*xval);
      }
    }
    else {
      if (sz) {
	for (int i=0; i<sz; ++i)
	  this->fastAccessDx(i) /= xval;
      }
    }
  }
  else {
    // Number of arguments
    const int N = Expr<S>::num_args;
    
    if (xsz) {
      
      T xval2 = xval*xval;
      
      if (sz) {

	if (x.hasFastAccess()) {
	  // Compute partials
	  FastLocalAccumOp< Expr<S> > op(x);
	  
	  // Compute each tangent direction
	  for(op.i=0; op.i<xsz; ++op.i) {
	    op.t = T(0.);
	    
	    // Automatically unrolled loop that computes
	    // for (int j=0; j<N; j++)
	    //   op.t += op.partials[j] * x.getTangent<j>(i);
	    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	    this->fastAccessDx(op.i) = 
	      (this->fastAccessDx(op.i) * xval - v * op.t) / xval2;
	  }
	}
	else {
	  // Compute partials
	  SlowLocalAccumOp< Expr<S> > op(x);
	  
	  // Compute each tangent direction
	  for(op.i=0; op.i<xsz; ++op.i) {
	    op.t = T(0.);
	    
	    // Automatically unrolled loop that computes
	    // for (int j=0; j<N; j++)
	    //   op.t += op.partials[j] * x.getTangent<j>(i);
	    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	    this->fastAccessDx(op.i) = 
	      (this->fastAccessDx(op.i) * xval - v * op.t) / xval2;
	  }
	}
	
      }
      
      else {
	
	this->resize(xsz);
	
	if (x.hasFastAccess()) {
	  // Compute partials
	  FastLocalAccumOp< Expr<S> > op(x);
	  
	  // Compute each tangent direction
	  for(op.i=0; op.i<xsz; ++op.i) {
	    op.t = T(0.);
	    
	    // Automatically unrolled loop that computes
	    // for (int j=0; j<N; j++)
	    //   op.t += op.partials[j] * x.getTangent<j>(i);
	    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	    this->fastAccessDx(op.i) = (-v * op.t) / xval2;
	  }
	}
	else {
	  // Compute partials
	  SlowLocalAccumOp< Expr<S> > op(x);
	  
	  // Compute each tangent direction
	  for(op.i=0; op.i<xsz; ++op.i) {
	    op.t = T(0.);
	    
	    // Automatically unrolled loop that computes
	    // for (int j=0; j<N; j++)
	    //   op.t += op.partials[j] * x.getTangent<j>(i);
	    Sacado::mpl::for_each< mpl::range_c< int, 0, N > > f(op);
	  
	    this->fastAccessDx(op.i) = (-v * op.t) / xval2;
	  }
	}
	
      }
      
    }
    
    else {
      
      if (sz) {
	for (int i=0; i<sz; ++i)
	  this->fastAccessDx(i) /= xval;
      }
      
    }

  }

  update_val_ = x.updateValue();
  if (update_val_)
    this->val() /= xval;

  return *this;
}

