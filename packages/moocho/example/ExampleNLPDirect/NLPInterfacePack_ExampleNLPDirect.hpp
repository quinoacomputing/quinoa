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

#ifndef EXAMPLE_NLP_FIRST_ORDER_DIRECT_H
#define EXAMPLE_NLP_FIRST_ORDER_DIRECT_H

#include "NLPInterfacePack_ExampleNLPObjGrad.hpp"
#include "NLPInterfacePack_NLPDirect.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "AbstractLinAlgPack_VectorSpaceBlocked.hpp"

namespace NLPInterfacePack {

/** \brief Simple example %NLP subclass to illustrate how to implement the
 * \c NLPDirect interface for a specialized \c NLP.
 *
 * For the NLP formulation see <tt>ExampleNLPObjGrad</tt>.
 *
 * In this %NLP we will select the first <tt>n/2</tt> variables as the dependent
 * or basis variables.  Note that the variable bounds are on the dependent
 * or independent variables variables if has_bounds = true.  Otherwise
 * there are no variable bounds.  Also the starting point can be
 * varied.
 *
 * The implementation of this %NLP subclass is actually independent from the vector
 * space for the dependent and independent variables (same vector space).
 * This implementation is defined entirely based on an arbitrary
 * <tt>VectorSpace</tt> object that is passed to the constructor
 * \c ExampleNLPDirect().  This %NLP subclass uses a
 * <tt>\ref AbstractLinAlgPack::VectorSpaceBlocked "VectorSpaceBlocked"</tt>
 * object to represent the space for <tt>[ x_dep; x_indep ]</tt>
 *
 * The quantities computed by this subclass include:
 \verbatim

  py = -inv(C)*c
   
  D = -inv(C)*N
  
  where:
  
    Gc' = [ C , N ]
  
      [ x(m+1) - 1									]
      [				x(m+2) - 1						]
    C = [							.					]
      [								.				]
      [									x(m+m) - 1	]
  
  
      [ x(1) - 10										]
      [				x(2) - 10						]
    N = [							.					]
      [								.				]
      [									x(m) - 10	]
  
 \endverbatim 
 * Here \c Gc is never computed explicitly.
 *
 * To make this possible this subclass relies on some specialized RTOp operators which
 * are implemented in C (for portability).  These operator classes are declared in the header
 * file <tt>NLPInterfacePack_ExampleNLPDirect.hpp</tt> and are documented \ref explnlp2_ops_grp "here".
 */
class ExampleNLPDirect
  : virtual public NLPDirect
  , virtual public ExampleNLPObjGrad
{
public:

  /** \brief Constructor.
   *
   * @param  vec_space  [in] Smart pointer to a vector space object that will
   *                    be used to define the spaces of dependent and independent
   *                    variables.
   * @param  xo         [in] The initial starting guess for \a x.
   * @param  has_bounds [in] If \c true, then the NLP will have bounds.  If \c false
   *                    then it will not have bounds.
   * @param  dep_bouned [in] If \c true, then the bounds will be on the dependent
   *                    variables.  If \c false, then the bounds will be on the
   *                    independent variable.  This argument is ignored if
   *                    <tt>has_bounds == false</tt>.
   */
  ExampleNLPDirect(
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

  /** @name Overridden public members from NLPDirect */
  //@{

  /** \brief . */
  Range1D var_dep() const;
  /** \brief . */
  Range1D var_indep() const;
  /** \brief . */
  const mat_fcty_ptr_t factory_D() const;
  /** \brief . */
  void calc_point(
    const Vector     &x
    ,value_type      *f
    ,VectorMutable   *c
    ,bool            recalc_c
    ,VectorMutable   *Gf
    ,VectorMutable   *py
    ,VectorMutable   *rGf
    ,MatrixOp        *GcU
    ,MatrixOp        *D
    ,MatrixOp        *Uz
    ) const;
  /** \brief . */
  void calc_semi_newton_step(
    const Vector    &x
    ,VectorMutable  *c
    ,bool           recalc_c
    ,VectorMutable  *py
    ) const;

  //@}

private:

  // /////////////////////////////////////////
  // Private data members

  mat_fcty_ptr_t   factory_D_;         // Matrix space object for D

  bool             initialized_;            // flag for if initialized has been called.

  // /////////////////////////////////////////
  // Private member functions

  /** \brief . */
  void assert_is_initialized() const;

};	// end class ExampleNLPDirect

// ///////////////////////////////////////////////
// Inline member functions

inline
void ExampleNLPDirect::assert_is_initialized() const
{
  typedef NLPInterfacePack::NLP NLP;
  if( !is_initialized() )
    throw NLP::UnInitialized("ExampleNLPDirect::assert_is_initialized() : Error, "
      "ExampleNLPDirect::initialize() has not been called yet." );
}

}	// end namespace NLPInterfacePack

#endif	// EXAMPLE_NLP_FIRST_ORDER_DIRECT_H
