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
// @HEADER

#ifndef SACADO_LFAD_LOGICALSPARSE_HPP
#define SACADO_LFAD_LOGICALSPARSE_HPP

#include "Sacado_LFad_LogicalSparseTraits.hpp"
#include "Sacado_LFad_ExpressionTraits.hpp"
#include "Sacado_Fad_DynamicStorage.hpp"

namespace Sacado {

  //! Namespace for logical forward-mode AD classes
  namespace LFad {

    //! Wrapper for a generic expression template
    /*!
     * This template class serves as a wrapper for all Fad expression
     * template classes.
     */
    template <typename ExprT> class Expr {};

    /*! 
     * \brief Implementation class for computing the logical sparsity of a 
     * derivative using forward-mode AD.
     */
    template <typename ValT, typename LogT>
    class LogicalSparseImp : 
      public Fad::DynamicStorage<ValT,LogT> {

      typedef Fad::DynamicStorage<ValT,LogT> Storage;

    public:

      //! Typename of values (e.g., double)
      typedef ValT value_type;

      //! Logical type (i.e., type for derivative array components (e.g., bool)
      typedef LogT logical_type;

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor
      LogicalSparseImp() : Storage(value_type(0)) {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      LogicalSparseImp(const value_type & x) : Storage(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      LogicalSparseImp(const int sz, const value_type & x) : Storage(sz, x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      LogicalSparseImp(const int sz, const int i, const value_type & x) : 
	Storage(sz, x) { 
	this->fastAccessDx(i)=logical_type(1); 
      }

      //! Copy constructor
      LogicalSparseImp(const LogicalSparseImp& x) : 
	Storage(x) {}

      //! Copy constructor from any Expression object
      template <typename S> LogicalSparseImp(const Expr<S>& x);

      //! Destructor
      ~LogicalSparseImp() {}

      //! Set %LogicalSparseImp object as the \c ith independent variable
      /*!
       * Sets the derivative array of length \c n to the \c ith row of the
       * identity matrix and has the same affect as the 
       * Implementation(const int sz, const int i, const T & x) 
       * constructor.
       */
      void diff(const int ith, const int n);

      //! Returns whether two LFad objects have the same values
      template <typename S>
      bool isEqualTo(const Expr<S>& x) const {
	typedef IsEqual<value_type> IE;
	if (x.size() != this->size()) return false;
	bool eq = IE::eval(x.val(), this->val());
	for (int i=0; i<this->size(); i++)
	  eq = eq && IE::eval(x.dx(i), this->dx(i));
	return eq;
      }

      //@}

      /*!
       * @name Derivative accessor methods
       */
      //@{

      //! Returns true if derivative array is not empty
      bool hasFastAccess() const { return this->size()!=0;}

      //! Returns true if derivative array is empty
      bool isPassive() const { return this->size()!=0;}
      
      //! Set whether variable is constant
      void setIsConstant(bool is_const) { 
	if (is_const && this->size()!=0)
	  this->resize(0);
      }
    
      //@}

      /*!
       * @name Assignment operators
       */
      //@{

      //! Assignment operator with constant right-hand-side
      LogicalSparseImp& operator=(const value_type& val);

      //! Assignment with Expr right-hand-side
      LogicalSparseImp& 
      operator=(const LogicalSparseImp& x);

      //! Assignment operator with any expression right-hand-side
      template <typename S> 
      LogicalSparseImp& operator=(const Expr<S>& x); 

      //@}

      /*!
       * @name Unary operators
       */
      //@{

      //! Addition-assignment operator with constant right-hand-side
      LogicalSparseImp& operator += (const value_type& x);

      //! Subtraction-assignment operator with constant right-hand-side
      LogicalSparseImp& operator -= (const value_type& x);

      //! Multiplication-assignment operator with constant right-hand-side
      LogicalSparseImp& operator *= (const value_type& x);

      //! Division-assignment operator with constant right-hand-side
      LogicalSparseImp& operator /= (const value_type& x);

      //! Addition-assignment operator with Expr right-hand-side
      template <typename S> 
      LogicalSparseImp& operator += (const Expr<S>& x);

      //! Subtraction-assignment operator with Expr right-hand-side
      template <typename S> 
      LogicalSparseImp& operator -= (const Expr<S>& x);
  
      //! Multiplication-assignment operator with Expr right-hand-side
      template <typename S> 
      LogicalSparseImp& operator *= (const Expr<S>& x);

      //! Division-assignment operator with Expr right-hand-side
      template <typename S> 
      LogicalSparseImp& operator /= (const Expr<S>& x);

      //@}

    }; // class LogicalSparseImp

    //! Expression template specialization for LogicalSparse
    template <typename ValT, typename LogT> 
    class Expr< LogicalSparseImp<ValT,LogT> > : 
      public LogicalSparseImp<ValT,LogT> {

    public:

      //! Default constructor
      Expr() : 
	LogicalSparseImp<ValT,LogT>() {}

      //! Constructor with supplied value \c x
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      Expr(const ValT & x) : 
	LogicalSparseImp<ValT,LogT>(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      Expr(const int sz, const ValT & x) : 
	LogicalSparseImp<ValT,LogT>(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      Expr(const int sz, const int i, const ValT & x) : 
	LogicalSparseImp<ValT,LogT>(sz,i,x) {}

      //! Copy constructor
      Expr(const Expr& x) : 
	LogicalSparseImp<ValT,LogT>(x) {}

      //! Copy constructor from any Expression object
      template <typename S> Expr(const Expr<S>& x) :
	LogicalSparseImp<ValT,LogT>(x) {}

      //! Destructor
      ~Expr() {}

    }; // class Expr<LogicalSparseImp>

    /*! 
     * \brief User inteface class for computing the logical sparsity pattern
     * of a derivative via forward-mode AD.
     */
    template <typename ValT, typename LogT >
    class LogicalSparse : public Expr< LogicalSparseImp<ValT,LogT > > {

    public:

      //! Turn LogicalSparse into a meta-function class usable with mpl::apply
      template <typename T, typename U = LogT> 
      struct apply {
	typedef LogicalSparse<T,U> type;
      };

      /*!
       * @name Initialization methods
       */
      //@{

      //! Default constructor.
      /*!
       * Initializes value to 0 and derivative array is empty
       */
      LogicalSparse() : 
	Expr< LogicalSparseImp< ValT,LogT > >() {}

      //! Constructor with supplied value \c x of type ValueT
      /*!
       * Initializes value to \c x and derivative array is empty
       */
      LogicalSparse(const ValT& x) : 
	Expr< LogicalSparseImp< ValT,LogT > >(x) {}

      //! Constructor with size \c sz and value \c x
      /*!
       * Initializes value to \c x and derivative array 0 of length \c sz
       */
      LogicalSparse(const int sz, const ValT& x) : 
	Expr< LogicalSparseImp< ValT,LogT > >(sz,x) {}

      //! Constructor with size \c sz, index \c i, and value \c x
      /*!
       * Initializes value to \c x and derivative array of length \c sz
       * as row \c i of the identity matrix, i.e., sets derivative component
       * \c i to 1 and all other's to zero.
       */
      LogicalSparse(const int sz, const int i, const ValT & x) : 
	Expr< LogicalSparseImp< ValT,LogT > >(sz,i,x) {}

      //! Copy constructor
      LogicalSparse(const LogicalSparse& x) : 
	Expr< LogicalSparseImp< ValT,LogT > >(x) {}

      //! Copy constructor from any Expression object
      template <typename S> LogicalSparse(const Expr<S>& x) : 
	Expr< LogicalSparseImp< ValT,LogT > >(x) {}

      //@}

      //! Destructor
      ~LogicalSparse() {}

      //! Assignment operator with constant right-hand-side
      LogicalSparse& operator=(const ValT& v) {
	LogicalSparseImp< ValT,LogT >::operator=(v);
	return *this;
      }

      //! Assignment operator with LogicalSparse right-hand-side
      LogicalSparse& operator=(const LogicalSparse& x) {
	LogicalSparseImp< ValT,LogT >::operator=(static_cast<const LogicalSparseImp< ValT,LogT >&>(x));
	return *this;
      }

      //! Assignment operator with any expression right-hand-side
      template <typename S> LogicalSparse& operator=(const Expr<S>& x) 
      {
	LogicalSparseImp< ValT,LogT >::operator=(x);
	return *this;
      }
	
    }; // class LogicalSparse<ValT,LogT>

  } // namespace LFad

} // namespace Sacado

#include "Sacado_LFad_LogicalSparseImp.hpp"
#include "Sacado_LFad_LogicalSparseOps.hpp"

#endif // SACADO_LFAD_LOGICALSPARSE_HPP
