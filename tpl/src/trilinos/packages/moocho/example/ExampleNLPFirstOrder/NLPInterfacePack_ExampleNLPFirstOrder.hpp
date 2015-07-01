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

#ifndef EXAMPLE_NLP_FIRST_ORDER_INFO_H
#define EXAMPLE_NLP_FIRST_ORDER_INFO_H

#include "NLPInterfacePack_ExampleNLPObjGrad.hpp"
#include "NLPInterfacePack_NLPFirstOrder.hpp"

namespace NLPInterfacePack {

/** \brief Simple example %NLP subclass to illustrate how to implement the
 * \c NLPFirstOrder interface for a specialized \c NLP.
 *
 * The example %NLP we will use is a scalable problem where
 * the basis of the jacobian of the constraints is a diagonal
 * matrix.
 \verbatim

    min    f(x) = (1/2) * sum( x(i)^2, for i = 1..n )
    s.t.   c(x)(j) = x(j) * (x(m+j) -1) - 10 * x(m+j) = 0, for j = 1..m
          0.01 < x(i) < 20, for i = p...p+m

    where:
        m = n/2
        p = 1 if dep_bounded == true or m+1 if dep_bounded = false
 \endverbatim
 * This subclass inherits from the subclass
 * <tt>\ref NLPInterfacePack::ExampleNLPDirect "ExampleNLPDirect"</tt>
 * mostly out of lazyness but also to show how flexible these interfaces
 * can be using mutiple inheritance.
 *
 * ToDo: Finish documentation!
 */
class ExampleNLPFirstOrder
  : virtual public NLPFirstOrder
  , virtual public ExampleNLPObjGrad
{
public:

  /** \brief Constructor (see </tt>ExampleNLPDirect::ExampleNLPDirect()</tt>).
   */
  ExampleNLPFirstOrder(
    const VectorSpace::space_ptr_t&  vec_space
    ,value_type                      xo
    ,bool                            has_bounds
    ,bool                            dep_bounded
    );

  /** @name Overridden public members from NLP */
  //@{

  /** \brief . */
  void initialize(bool test_setup);
  /** \brief . */
  bool is_initialized() const;

  //@}

  /** @name Overridden public members from NLPFirstOrder */
  //@{

  /// Overridden to check the concrete type of Gc
  void set_Gc(MatrixOp* Gc);
  /** \brief . */
  const NLPFirstOrder::mat_fcty_ptr_t factory_Gc() const;
  /// Returns an ExampleBasisSystem
  const basis_sys_ptr_t basis_sys() const;

  //@}

protected:

  /** @name Overridden protected members from NLPFirstOrder */
  //@{

  /** \brief . */
  void imp_calc_Gc(const Vector& x, bool newx, const FirstOrderInfo& first_order_info) const;

  //@}

private:

  // /////////////////////////////////////////
  // Private data members

  bool                                initialized_;  // flag for if initialized has been called.
  NLPFirstOrder::mat_fcty_ptr_t       factory_Gc_;   // Factory for Gc
  NLPFirstOrder::basis_sys_ptr_t      basis_sys_;    // The basis system

  // /////////////////////////////////////////
  // Private member functions

  /** \brief . */
  void assert_is_initialized() const;

};	// end class ExampleNLPFirstOrder

// ///////////////////////////////////////////////
// Inline member functions

inline
void ExampleNLPFirstOrder::assert_is_initialized() const
{
  typedef NLPInterfacePack::NLP NLP;
  if( !is_initialized() )
    throw NLP::UnInitialized("ExampleNLPFirstOrder::assert_is_initialized() : Error, "
      "ExampleNLPFirstOrder::initialize() has not been called yet." );
}

}	// end namespace NLPInterfacePack

#endif	// EXAMPLE_NLP_FIRST_ORDER_INFO_H
