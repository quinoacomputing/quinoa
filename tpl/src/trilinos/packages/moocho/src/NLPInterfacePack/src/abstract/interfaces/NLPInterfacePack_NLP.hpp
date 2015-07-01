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

#ifndef NLP_H
#define NLP_H

#include <stdexcept>
#include <string>

#include "NLPInterfacePack_Types.hpp"
#include "AbstractLinAlgPack_VectorMutable.hpp"
#include "AbstractLinAlgPack_Permutation.hpp"
#include "StandardCompositionRelationshipsPack.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObject.hpp"

namespace OptionsFromStreamPack {
  class OptionsFromStream;
}

namespace NLPInterfacePack {

/** \brief %NLP interface class {abstract}.
 *
 * <b>Overview:</b>
 *
 * This class represents an abstract interface to a general nonlinear programming problem of the form
 * (in mathematical and ASCII notation):
 \f[
 \begin{array}{lcl}
 \mbox{min}  &  & f(x)                     \\
 \mbox{s.t.} &  & c(x) = 0                 \\
             &  & x^L \leq x    \leq x^U
 \end{array}
 \f]
 * where:<ul>
 * <li> \f$x, x^L, x^U \:\in\:\mathcal{X}\f$
 * <li> \f$f(x) : \:\mathcal{X} \rightarrow \Re\f$
 * <li> \f$c(x) : \:\mathcal{X} \rightarrow \mathcal{C}\f$
 * <li> \f$\mathcal{X} \:\in\:\Re\:^n\f$
 * <li> \f$\mathcal{C} \:\in\:\Re\:^m\f$
 * </ul>
 \verbatim

     min    f(x)
     s.t.   c(x) = 0
            xl <= x <= xu
  where:
            x    <: space_x
            f(x) <: space_x -> R^1
            c(x) <: space_x -> space_c 
            space_x <: R^n
            space_c <: R^n -> R^m
 \endverbatim

 * The %NLP is defined in terms of vector spaces for the unknowns \a x (\c space_x),
 * the equality constraints \a c (\c space_c)and nonlinear operator functions
 * \a f(x) and \a c(x).  In the above form, none of the variables are fixed between
 * bounds (strictly xl < xu).  It is allowed however for <tt>m == 0</tt> for the
 * elimination of general constriants.  It is also allowed for <tt>n == m</tt>
 * in which case <tt>this</tt> represents a fully determined system of nonlinear
 * equaltions.  In any case, an objective function is always
 * included in the formutation and will impact solution algorithms.
 *
 * Special types of NLPs are identified as: <ol>
 * <li> Fully general %NLP :
 *      <ul><li><tt>( xl != -Inf || xu != +Inf ) && ( m  > 0 )</tt></ul>
 * <li> General equality only constrained %NLP :
 *      <ul><li><tt>( xl == -Inf && xu == +Inf ) && ( m  > 0 )</tt></ul>
 * <li> Bound constrained %NLP :
 *      <ul><li><tt>( xl != -Inf || xu != +Inf ) && ( m == 0 )</tt></ul>
 * <li> Unconstrained %NLP :
 *      <ul><li><tt>( xl == -Inf && xu == +Inf ) && ( m == 0 )</tt></ul>
 * <li> Nonlinear Equations (NLE) :
 *      <ul><li><tt>n == m</tt></ul>
 * </ol>
 *
 * If <tt>n==m</tt> but some of the equations in <tt>c(x)</tt> are
 * dependent (but consistent) then the problem is actually an NLP
 * and not an NLE and it is possible that some of the variable bounds
 * my be active at the solution but in general they can not be.
 * An optimization algorithm may refuse to solve some of the above problems but this
 * interface allows all of these different types of mathematical programming problems
 * to be represented using this interface.
 *
 * The Lagrangian for this problem is defined by:
 \verbatim

  L = f(x) + lambda' * c(x) + nul' * ( xl - x ) + nuu' * ( x - xu )
 \endverbatim
 * The optimality conditions are given by:
 \verbatim

  del(L,x) = del(f,x) + del(c,x) * lambda + nu = 0
  c(x) = 0
  nuu(i) * ( x(i) - xu(i) ) = 0,      for i = 1...n
  nuu(i) * ( x(i) - xu(i) ) = 0,      for i = 1...n
  where:
    nu = nuu - nul
 \endverbatim
 * What is unique about this interface is that the vector objects are hidden behind
 * abstact interfaces.  Clients can create vectors from the various vector spaces
 * using the <tt>\ref AbstractLinAlgPack::VectorSpace "VectorSpace"</tt> objects returned from
 * <tt>this->space_x()</tt> (dim \c n), and <tt>this->space_c()</tt> (dim \c m).
 * In this sense, an <tt>%NLP</tt> object
 * acks as an "Abstract Factory" to create the vectors needed by an optimization
 * algorithm.  This allows optimization software to be written in a way that is
 * completly independent from the linear algebra components.
 *
 * <b>General Inequalities and Slacks</b>
 *
 * The underlying NLP may containe general inequality
 * constraints which where converted to equalities using slack
 * variables.  The actual underlying NLP may take the form
 * (in mathematical and ASCII notation):
 \f[
 \begin{array}{lcl}
 \mbox{min}  &  & \hat{f}(\hat{x})                                 \\
 \mbox{s.t.} &  & \hat{c}(\hat{x}) = 0                             \\
             &  & \hat{h}^L \leq \hat{h}(\hat{x}) \leq \hat{h}^U   \\
             &  & \hat{x}^L \leq \hat{x} \leq \hat{x}^U
 \end{array}
 \f]
 * where:<ul>
 * <li> \f$\hat{x}, \hat{x}^L, \hat{x}^U \:\in\:\hat{\mathcal{X}}\f$
 * <li> \f$\hat{f}(\hat{x}) : \:\hat{\mathcal{X}} \rightarrow \Re\f$
 * <li> \f$\hat{c}(\hat{x}) : \:\hat{\mathcal{X}} \rightarrow \hat{\mathcal{C}}\f$
 * <li> \f$h(x) : \:\hat{\mathcal{X}} \rightarrow \hat{\mathcal{H}}\f$
 * <li> \f$\hat{\mathcal{X}} \:\in\:\Re\:^{\hat{n}}\f$
 * <li> \f$\hat{\mathcal{C}} \:\in\:\Re\:^{\hat{m}}\f$
 * <li> \f$\hat{\mathcal{H}} \:\in\:\Re\:^{\hat{m^I}}\f$
 * </ul>
 \verbatim

     min    f_breve(x_breve)
     s.t.   c_breve(x_breve) = 0
            hl_breve <= h_breve(x_breve) <= hu_breve
          xl_breve <= x_breve <= xu_breve
  where:
          x_breve    <: space_x
            f_breve(x_breve) <: space_x_breve -> R^1
          c_breve(x_breve) <: space_x_breve -> space_c_breve 
          h_breve(x_breve) <: space_x_breve -> space_h_breve
            space_x_breve <: R^n_breve
            space_c_breve <: R^n -> R^m_breve
            space_h_breve <: R^n -> R^mI_breve
 \endverbatim
 *
 * ToDo: Finish!
 *
 * <b>Client Usage:</b>
 *
 * Before an %NLP object can be used, the <tt>initialize()</tt> method must be called to
 * make sure that all of the initializations needed for the NLP have been performed.
 * This method also resets counters an other information.  Before calling <tt>initialize()</tt>
 * the client can specify whether the initial point for \c x must be in bounds by calling
 * <tt>force_xinit_in_bounds(bool)</tt>.
 *
 * Smart reference counted pointers to the three vector spaces for \a x, \a c(x) and \a h(x)
 * are returned by the methods <tt>space_x()</tt>, <tt>space_c()</tt> and <tt>space_h()</tt>
 * respectively.  The vector space objects returned by these methods are ment to be more
 * than transient.  In fact, it is expected that these vector space objects should remain
 * valid for the entire run of an NLP algorithm.  Only if the underlying NLP is changed in
 * a fundamental way (i.e. \c n, or \c m changes) should the vector space objects returned
 * from these function become invalid.  In this case the client must call these methods again
 * to get updated vector space objects.
 *
 * The dimensionality of the NLP is returned by the methods <tt>n()</tt> and <tt>m()</tt>
 * but they have default implementations based on <tt>space_x()</tt> and <tt>space_c()</tt>
 * respectively.
 *
 * The number of variables \c x with finite bounds is returned by the method <tt>num_bounded_x()</tt>.
 * The methods <tt>xl()</tt> and <tt>xu()</tt> return references
 * to vector objects representing these bounds.  A lower bound is considered infinite if
 * <tt>xl().get_ele(i) == -infinite_bound()</tt> and an upper bound is considered infinite if
 * <tt>xu().get_ele(i) == +infinite_bound()</tt>.
 *
 * If <tt>ns() > 0</tt>, the methods \c hl_breve() and \c hu_breve() return references to the upper and lower
 * bounds to the general inequality constraints \a h_breve(x_breve).  While it is expected that 
 * <tt>hl_breve().get_ele(j) != -infinite_bound() || hu_breve().get_ele(j) != +infinite_bound</tt> for <tt>j = 1...ns()</tt>
 * this is not required by this interface.  On the other hand it seems silly to define general inequality constriants
 * that are not bounded but there may be some reason to include these that makes things easier for the
 * implementor of the NLP subclass.
 *
 * The initial guess for the unknowns \a x (primal variables) is returned by the method <tt>xinit()</tt>.
 * Vectors containing the initial guesses for the Lagrange multipliers can be obtained by calling the
 * method <tt>get_init_lagrange_mult()</tt>.  
 * 
 * The bread and butter of an %NLP interface is the calculation of the functions that define the objective
 * and constraints and various points \a x using the methods \c calc_f() and \c calc_c().
 * The quantities that these functions update must be set prior by calling the methods
 * \c set_f() and \c set_c() respectively.  It may seem strange not to pass these quantities
 * directly to the calculation functions but there is a good reason for this.
 * The reason is that this interface supports the efficient update of mutiple quantities at a given
 * point \a x.  For example, in many %NLPs some of the same terms are shared between the constriants
 * functions and objective function.  Therefore it is more efficient to compute \a f(x) and \a c(x)
 * simultaneously rather than computing them separately.  In order to allow for this possibility,
 * the client can set the desired quantities (i.e. \c set_f() and \c set_c()) prior to calling \c calc_f()
 * and \c calc_c().  Then, whatever quantities that have been set will be computed my any call to
 * <tt>calc_?()</tt> method.
 *
 * Once an optimization algorithm has the solution (or gives up with a suboptimal point), it should
 * report this solution to the %NLP object using the method \c report_final_solution().
 *
 * Finally, the client can get the counts for the number of function evaluations since \c initialize()
 * was called using the methods \c num_f_evals() and \c num_c_evals().
 * These counts do not include any function evaluations that may have been used internally for finite
 * difference evaluations or anything of that nature.  Also, if the client calls <tt>calc_info(x,false)</tt>
 * (where \c info = \c f or \c c) several times then the default implementation will increment the count
 * for each call even though the actual quantity may not actually be recalculated each time (i.e. if <tt>newx==false</tt>).
 * This information is not known here in this base class but the subclasses can overide this behavior if desired.
 *
 * <b>Subclass developer's notes:</b>
 *
 * The calcuation methods \c calc_f() and \c calc_c() have default implementations in this base
 * class that should meet the needs of all the subclasses.  These method implementations call the protected
 * pure virtual methods \c imp_calc_f() and \c imp_calc_c() to compute the actual quantities.
 * Subclasses must override these methods (in addition to several methods from the public interface).
 * Pointers to the quantities to be updated are passed to these methods to the subclasses in the form of
 * an aggregate \c ZeroOrderInfo object that is returned from the protected method \c zero_order_info().
 * This ensures that the only interaction between an NLP base object and its subclass objects is through
 * member functions and never through pulic or protected data members.
 *
 * <A NAME="must_override"></A>
 * The following methods must be overridden by a subclass in order to create a concrete NLP object:
 * \c force_xinit_in_bounds(bool), \c force_xinit_in_bounds(), \c is_initialized(),
 * \c space_x(), \c space_c(), \c num_bounded_x(), \c xl(), xu(), \c xinit(),]
 * \c scale_f(value_type), \c scale_f().
 *
 * <A NAME="should_override"></A>
 * The following methods should be overridden by most subclasses but do not have to be: \c initialize(),
 * \c get_init_lagrange_mult(), \c report_final_solution().
 *
 * The following methods should never have to be overridden by most subclasses except in some very
 * strange situations: \c set_f(), \c get_f(), \c f(), \c set_c(), \c get_c(), \c c(),
 * \c calc_f(), \c calc_c(), \c num_f_evals(), \c num_c_evals().
 *
 * <b>Additional notes:</b>
 *
 * The bounds on the variables can play a very critical role in many optimization algorithms.  It is desirable
 * for the functions \a f(x) and \a c(x) to be defined and
 * relatively well behaved (i.e. smooth, continous, differentiable etc.) in the region <tt>xl <= x <= xu</tt>.
 * While this is not always possible, an NLP can often be reformulated to have this properly.  For example, suppose
 * there are constraints has the form:
 \verbatim

 log(x(1) - x(5)) - 4 == 0
 x(1) - x(5) >= 1e-8
 \endverbatim
 * where \c x(1) and x(5) are unbounded.  It is clear that if <tt>x(1) < x(5)</tt> that this constraint will be
 * undefined and will return \c NaN on most computers (if you are lucky).  The constraint <tt>x(1) < x(5)</tt> is
 * very hard to enforce at every iteration in an NLP solver so that this will not happen   A better approach would be
 * to add an extra varaible (say \c x(51) for an %NLP with <tt>n == 50</tt>) and add an extra constraint:
 \verbatim

 log(x(51)) == 0
 x(51) - (x(1) - x(5)) == 0
 x(51) >= 1e-8
 \endverbatim
 * In the above expanded formulation the simple bound <tt>x(51) >= 1e-8</tt> is easy to inforce and these undefined
 * regions can be avoided.  While the property that \a f(x), \a c(x) and \a h(x) being bounded for all
 * <tt>x <: { x | xl <= x <= x}</tt> is a desireable properly, this is not required by this interface.  As a result
 * the client should be prepaired to deal with return values of \c NaN or \c Inf for \c f, \c c and \c h.
 */
class NLP : virtual public Teuchos::VerboseObject<NLP> {
public:
  
  typedef AbstractLinAlgPack::Vector         Vector;         // doxygen likes typedef?
  typedef AbstractLinAlgPack::VectorMutable  VectorMutable;  // doxygen likes typedef?
  
  /** \brief . */
  typedef Teuchos::RCP<const VectorSpace>  vec_space_ptr_t;

  /** \brief . */
  typedef Teuchos::RCP<
    const OptionsFromStreamPack::OptionsFromStream>             options_ptr_t;

  /** @name exceptions */
  //@{

  /// Thrown if any member functions are called before initialize() has been called.
  class UnInitialized : public std::logic_error
  {public: UnInitialized(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown from <tt>initialize()</tt> if some logical error occured
  class InvalidInitialization : public std::logic_error
  {public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if an incompatible object is used
  class IncompatibleType : public std::logic_error
  {public: IncompatibleType(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown some bounds do not existe
  class NoBounds : public std::logic_error
  {public: NoBounds(const std::string& what_arg) : std::logic_error(what_arg) {}};

  //@}

  /// Value for an infinite bound.
  static value_type infinite_bound();

  /** @name Constructors, Destructor */
  //@{

  /// Initialize to no reference set to calculation quanities
  NLP();
  /// Destructor that cleans all the memory it owns
  virtual ~NLP();

  //@}
  
  /** @name NLP initialization */
  //@{

  /** \brief Set if the initial point must be within the bounds.
   *
   * This method must be called before <tt>this->initialize()</tt> is called.
   *
   * Postconditions:<ul>
   * <li> <tt>this->is_initialized() == false</tt>
   * </ul>
   */
  virtual void force_xinit_in_bounds(bool force_xinit_in_bounds) = 0;
  /** \brief Returns if the initial point must be within the bounds.
    */
  virtual bool force_xinit_in_bounds() const = 0;
  /** \brief Set the options that <tt>this</tt> %NLP may be interested in.
   *
   * Note that it is allowed for the client to alter <tt>*options.get()</tt> after
   * this method is called so <tt>this</tt> had better read the options inside of
   * the <tt>this->initialize()</tt> method.
   *
   * The default implementation is to just ignore these options.
   *
   * Note that if the subclass overrides this method then it must also override
   * the <tt>get_options()</tt> method.
   */
  virtual void set_options( const options_ptr_t& options );
  /** \brief Get the <tt>OptionsFromStream</tt> object being used to extract the options from.
   *
   * The default implementation returns <tt>return.get() == NULL</tt>.
   */
  virtual const options_ptr_t& get_options() const;
  /** \brief Initialize the NLP before it is used.
   *
   * Postconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt>
   * <li> [<tt>this->force_xinit_in_bounds()==true</tt>]
   *      <tt>this->xl() <= this->xinit() <= this->xu()</tt>
   * <li> <tt>this->num_f_evals() == 0</tt>
   * <li> [<tt>this->m() > 0</tt>] <tt>this->num_c_evals() == 0</tt>
   * </ul>
   *
   * Note that subclasses must call this function to reset what needs to be
   * reset in this base object.
   */
  virtual void initialize( bool test_setup = false );
  /** \brief Return if <tt>this</tt> is initialized.
    */
  virtual bool is_initialized() const = 0;

  //@}

  /** @name Dimensionality. */
  //@{

  /** \brief Return the number of variables.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Default implementation returns <tt>this-space_x()->dim()</tt>.
   */
  virtual size_type n() const;
  /** \brief Return the number of general equality constraints.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Default implementation returns <tt>( this->space_c().get() != NULL ? this-space_c()->dim() : 0 )</tt>.
   */
  virtual size_type m() const;

  //@}

  /** @name Vector space objects */
  //@{

  /** \brief Vector space object for unknown variables x (dimension n).
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>return.get() != NULL</tt>
   * </ul>
   */
  virtual vec_space_ptr_t space_x() const = 0;
  /** \brief Vector space object for general equality constraints c(x) (dimension m).
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>this->m() > 0</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->m() == 0</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   */
  virtual vec_space_ptr_t space_c() const = 0;

  //@}

  /** @name Bounds on the unknown variables x. */
  //@{

  /** \brief Returns the number of variables in <tt>x(i)</tt> for which <tt>xl(i)> -infinite_bound()</tt>
   * or <tt>xu(i) < +infinite_bound()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   */
  virtual size_type num_bounded_x() const = 0;
  /** \brief Returns the lower bounds on the variables <tt>x</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Any bounds that are non-existant will return <tt>this->xl().get_ele(i) == -NLP::infinite_bound()</tt>.
   */
  virtual const Vector& xl() const = 0;
  /** \brief Returns a reference to the vector of upper bounds on the variables <tt>x</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Any bounds that are non-existant will return <tt>this->xu().get_ele(i) == +NLP::infinite_bound()</tt>.
   */
  virtual const Vector& xu() const = 0;

  /** \brief Set the maximum absolute value for which the variable bounds may be violated
   * by when computing function and gradient values.
   *
   * In other words the client should never never call on the NLP to compute
   * a function and gradient evaluation outside of:
   \verbatim

     xl - max_var_bounds_viol <= x <= xu + max_var_bounds_viol
   \endverbatim
   */
  virtual value_type max_var_bounds_viol() const = 0;

  //@}

  /** @name Initial guess of NLP solution */
  //@{

  /** \brief Returns a reference to the vector of the initial guess for the solution <tt>x</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>return.space().is_compatible(*this->space_x()) == true)</tt>
   * </ul>
   */
  virtual const Vector& xinit() const = 0;
  /** \brief Get the initial value of the Lagrange multipliers lambda.
   *
   * By default this function just sets them to zero.
   *
   *
   * @param lambda  [out] Pointer to lagrange multipliers for equalities.
   *                lambda == NULL is allowed in which case it will not
   *                be set.  Must have been created by <tt>this->space_c()->create_member()</tt>.
   *                Must be NULL if m() == 0.
   * @param nu      [out] Pointer to lagrange multipliers for bounds.
   *                nu == NULL is allowed in which case it will not
   *                be set.  Must have been created by <tt>this->space_x()->create_member()</tt>.
   *                Must be NULL if num_bounded_x() == 0.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   */
  virtual void get_init_lagrange_mult(
    VectorMutable*   lambda
    ,VectorMutable*  nu
    ) const;

  //@}

  /** @name Set and access storage for the objective function value f(x). */
  //@{

  /** \brief Set a pointer to an value to be updated when <tt>this->calc_f()</tt> is called.
   *
   * @param  f  [in] Pointer to objective function value.  May be \c NULL.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get_f() == f</tt>
   * </ul>
   */
  virtual void set_f(value_type* f);
  /** \brief Return pointer passed to <tt>this->set_f()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   */
  virtual value_type* get_f();
  /** \brief Returns non-<tt>const</tt> <tt>*this->get_f()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_f() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual value_type& f();
  /** \brief Returns <tt>const</tt> <tt>*this->get_f()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_f() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
   virtual const value_type& f() const;

  //@}

  /** @name Set and access storage for the residual of the general equality constriants c(x). */
  //@{

  /** \brief Set a pointer to a vector to be updated when <tt>this->calc_c()</tt> is called.
   *
   * @param  c  [in] Pointer to constraint residual vector.  May be \c NULL.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> [<tt>c != NULL</tt>] <tt>c->space().is_compatible(*this->space_c()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get_c() == c</tt>
   * </ul>
   */
  virtual void set_c(VectorMutable* c);
  /** \brief Return pointer passed to <tt>this->set_c()</tt>.
   */
  virtual VectorMutable* get_c();
  /** \brief Returns non-<tt>const</tt> <tt>*this->get_c()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_c() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual VectorMutable& c();
  /** \brief Returns <tt>const</tt> <tt>*this->get_c()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_c() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual const Vector& c() const;

  //@}

  /** @name Unset calculation quantities */
  //@{
  
  /** \brief Call to unset all storage quantities (both in this class and all subclasses).
   *
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get_f() == NULL</tt>
   * <li> <tt>this->get_c() == NULL</tt>
   * <li> <tt>this->get_c_breve() == NULL</tt>
   * <li> <tt>this->get_h_breve() == NULL</tt>
   * </ul>
   *
   * This method must be called by all subclasses that override it.
   */
  virtual void unset_quantities();

  //@}

  /** @name Calculation members */
  //@{

  /** \brief Set the scaling of the objective function.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->scale_f() == true</tt>
   * </ul>
   */
  virtual void scale_f( value_type scale_f ) = 0;
  /** \brief Return the scaling being used for the objective function.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   */
  virtual value_type scale_f() const = 0;
  /** \brief Update the value for the objective <tt>f</tt> at the point <tt>x</tt> and put it in the stored reference.
   *
   * @param  x     [in] Point at which to calculate the object function <tt>f</tt>.
   * @param  newx  [in] (default \c true) If \c false, the values in \c x are assumed to be the same as
   *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
   *               If \c true, the values in \c x are assumed to not be the same as the last call to a
   *               <tt>this->calc_*(x,newx)</tt> member.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>this->get_f() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->f()</tt> is updated to \a f(x)
   * </ul>
   *
   * The storage reference for <tt>c</tt> may also be updated at this point (if <tt>get_c() != NULL</tt>)
   * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
   * to be updated as a side effect.
   */ 
  virtual void calc_f(const Vector& x, bool newx = true) const;
  /** \brief Update the constraint residual vector for <tt>c</tt> at the point <tt>x</tt> and put it in the stored reference.
   *
   * @param  x     [in] Point at which to calculate residual to the equality constraints <tt>c</tt>.
   * @param  newx  [in] (default \c true) If \c false, the values in \c x are assumed to be the same as
   *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
   *               If \c true, the values in \c x are assumed to not be the same as the last call to a
   *               <tt>this->calc_*(x,newx)</tt> member.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>this->get_c() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->c()</tt> is updated to \a c(x)
   * </ul>
   *
   * The storage reference for <tt>f</tt> may also be updated at this point (if <tt>get_f() != NULL</tt>)
   * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
   * to be updated as a side effect.
   */ 
  virtual void calc_c(const Vector& x, bool newx = true) const;

  //@}

  /** @name Report final solution */
  //@{

  /** \brief Used by the solver to report the final solution and multipliers.
   *
   * Call this function to report the final solution of the
   * unknows x and the Lagrange multipliers for the
   * equality constriants <tt>lambda</tt> and the varaible bounds
   * <tt>nu</tt>.  If any of the Lagrange multipliers
   * are not known then you can pass <tt>NULL</tt> in for them.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * The default behavior is to just ignore this.
   */
  virtual void report_final_solution(
    const Vector&    x
    ,const Vector*   lambda
    ,const Vector*   nu
    ,bool            is_optimal
    );

  //@}

  /** @name Objective and constraint function evaluation counts. */
  //@{

  /** \brief Gives the number of object function f(x) evaluations called by the solver
   * since initialize() was called.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   */
  virtual size_type num_f_evals() const;
  /** \brief Gives the number of constraint function c(x) evaluations called by the solver
   * since initialize() was called.  Throws exception if <tt>this->m() == 0</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   */
  virtual size_type num_c_evals() const;

  //@}

  /** @name General inequalities and slack variables */
  //@{

  /** \brief Return the number of slack variables (i.e. number of general inequalities).
   *
   * Default implementation returns
   * <tt>(this->space_h_breve().get() ? this->space_h_breve()->dim() : 0)</tt>.
   */
  virtual size_type ns() const;

  /** \brief Vector space object for the original equalities <tt>c_breve(x_breve)</tt>
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>this->m() - this->ns() > 0</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->m() - this->ns() == 0</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   *
   * The default implementation returns <tt>this->space_c()</tt>.
   */
  virtual vec_space_ptr_t space_c_breve() const;

  /** \brief Vector space object for the original inequalities <tt>h_breve(x_breve)</tt>
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> [<tt>this->ns() > 0</tt>] <tt>return.get() != NULL</tt>
   * <li> [<tt>this->ns() == 0</tt>] <tt>return.get() == NULL</tt>
   * </ul>
   *
   * The default implementation returns <tt>return.get() == NULL</tt>.
   */
  virtual vec_space_ptr_t space_h_breve() const;

  /** \brief Returns a reference to the vector of lower bounds on the general inequality constraints <tt>h_breve(x_breve)</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->ns() > 0</tt> (throw <tt>std::logic_error</tt>)
   * </ul>
   *
   * Any bounds that are non-existant will return <tt>this->hl_breve().get_ele(i) == -NLP::infinite_bound()</tt>.
   *
   * The default implementation throws an exception.
   */
  virtual const Vector& hl_breve() const;

  /** \brief Returns a reference to the vector of upper bounds on the general inequality constraints <tt>h_breve(x_breve)</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Any bounds that are non-existant will return <tt>this->hu_breve().get_ele(i) == +NLP::infinite_bound()</tt>.
   *
   * The default implementation throws an exception.
   */
  virtual const Vector& hu_breve() const;

  /** \brief Set a pointer to a vector to be updated when <tt>this->calc_c_breve()</tt> is called.
   *
   * @param  c_breve  [in] Pointer to constraint residual vector.  May be \c NULL.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> [<tt>c != NULL</tt>] <tt>c->space().is_compatible(*this->space_c_breve()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get_c_breve() == c_breve</tt>
   * </ul>
   */
  virtual void set_c_breve(VectorMutable* c_breve);
  /** \brief Return pointer passed to <tt>this->set_c_breve()</tt>.
   */
  virtual VectorMutable* get_c_breve();
  /** \brief Returns non-<tt>const</tt> <tt>*this->get_c_breve()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_c() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual VectorMutable& c_breve();
  /** \brief Returns <tt>const</tt> <tt>*this->get_c_breve()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_c_breve() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual const Vector& c_breve() const;

  /** \brief Set a pointer to a vector to be updated when <tt>this->calc_h_breve()</tt> is called.
   *
   * @param  h_breve  [in] Pointer to constraint residual vector.  May be \c NULL.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> [<tt>c != NULL</tt>] <tt>c->space().is_compatible(*this->space_h_breve()) == true</tt>
   *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->get_h_breve() == h_breve</tt>
   * </ul>
   */
  virtual void set_h_breve(VectorMutable* h_breve);
  /** \brief Return pointer passed to <tt>this->set_h_breve()</tt>.
   */
  virtual VectorMutable* get_h_breve();
  /** \brief Returns non-<tt>const</tt> <tt>*this->get_h_breve()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_c() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual VectorMutable& h_breve();
  /** \brief Returns <tt>const</tt> <tt>*this->get_h_breve()</tt>.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->get_h_breve() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   */
  virtual const Vector& h_breve() const;

  /** \brief Return the permutation object for the variables.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>return.space()is_compatible(*this->space_x()) == true</tt>
   * </ul>
   *
   * The default returns <tt>return.is_identity() == true</tt>
   */
  virtual const Permutation& P_var() const;

  /** \brief Return the permutation object for the constraints.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>this->m() > 0</tt> (throw <tt>std::logic_error</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>return.space()is_compatible(*this->space_c()) == true</tt>
   * </ul>
   *
   * The default returns <tt>return.is_identity() == true</tt>
   */
  virtual const Permutation& P_equ() const;

  /** \brief Update the constraint residual vector for <tt>c_breve</tt> at the point <tt>x</tt> and put it
   * in the stored reference.
   *
   * @param  x     [in] Point at which to calculate residual to the equality constraints <tt>c_breve</tt>.
   * @param  newx  [in] (default \c true) If \c true, the values in \c x are the same as
   *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
   *               If \c false, the values in \c x are not the same as the last call to a
   *               <tt>this->calc_*(x,newx)</tt> member.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>this->get_c_breve() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->c_breve()</tt> is updated to \a c_breve(x_breve)
   * </ul>
   *
   * The storage reference for <tt>f</tt> and/or <tt>h_breve</tt> may also be updated at this point
   * (if <tt>get_f() != NULL</tt> and/or <tt>get_h_breve() != NULL</tt>) but is not guarentied to be.
   * But no other quanities from possible subclasses are allowed to be updated as a side effect.
   */ 
  virtual void calc_c_breve(const Vector& x, bool newx = true) const;

  /** \brief Update the constraint residual vector for <tt>h_breve</tt> at the point <tt>x</tt> and put it
   * in the stored reference.
   *
   * @param  x     [in] Point at which to calculate residual to the equality constraints <tt>h_breve</tt>.
   * @param  newx  [in] (default \c true) If \c true, the values in \c x are the same as
   *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
   *               If \c false, the values in \c x are not the same as the last call to a
   *               <tt>this->calc_*(x,newx)</tt> member.
   *
   * Preconditions:<ul>
   * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
   * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
   * <li> <tt>this->get_h_breve() != NULL</tt> (throw <tt>NoRefSet</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->h_breve()</tt> is updated to \a h_breve(x_breve)
   * </ul>
   *
   * The storage reference for <tt>f</tt> and/or <tt>c_breve</tt> may also be updated at this point
   * (if <tt>get_f() != NULL</tt> and/or <tt>get_c_breve() != NULL</tt>) but is not guarentied to be.
   * But no other quanities from possible subclasses are allowed to be updated as a side effect.
   */ 
  virtual void calc_h_breve(const Vector& x, bool newx = true) const;

  //@}

  /** \brief Struct for objective and constriants (pointer).
   *
   * Objects of this type are passed on to subclasses and contain pointers to
   * quantities to be updated.
   */
  struct ZeroOrderInfo {
  public:
    /** \brief . */
        ZeroOrderInfo() : f(NULL), c(NULL), h(NULL)
    {}
    /** \brief . */
    ZeroOrderInfo( value_type* f_in, VectorMutable* c_in, VectorMutable* h_in )
      : f(f_in), c(c_in), h(h_in)
    {}
    /// Pointer to objective function <tt>f</tt> (Will be NULL if not set)
    value_type*           f;
    /// Pointer to constraints residual <tt>c</tt> (Will be NULL if not set)
    VectorMutable*  c;
    /// Pointer to inequality constraints <tt>h</tt> (Will be NULL if not set)
    VectorMutable*  h;
  }; // end struct ZeroOrderInfo

  /// Return pointer to set quantities
  const ZeroOrderInfo zero_order_info() const;

  /// Return pointer to set <tt>hat</tt> quantities
  const ZeroOrderInfo zero_order_info_breve() const;

protected:

  /** @name Protected methods to be overridden by subclasses */
  //@{

  /** \brief Overridden to compute f(x) (and perhaps other quantities if set).
   *
   * Preconditions:<ul>
   * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
   * <li> <tt>zero_order_info.f != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*zero_order_info.f</tt> is updated to \a f(x).
   * </ul>
   *
   * @param x       [in]  Unknown vector (size n).
   * @param newx    [in]  True if is a new point.
   * @param zero_order_info
   *                [out] Pointers to \c f, \c c and \c h.
   *                On output, <tt>*zero_order_info.f</tt> is updated to \a f(x)
   *                If <tt>this->multi_calc() == true</tt> then
   *                any of the other quantities pointed to in \c zero_order_info may be set on
   *                output, but are not guaranteed to be.
   */
  virtual void imp_calc_f(const Vector& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;
  /** \brief Overridden to compute c(x) and perhaps f(x) and/or h(x) (if multiple calculaiton = true).
   *
   * Preconditions:<ul>
   * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
   * <li> <tt>zero_order_info.c != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*zero_order_info.c</tt> is updated to c(x).
   * </ul>
   *
   * @param x       [in]  Unknown vector (size n).
   * @param newx    [in]  True if is a new point.
   * @param zero_order_info
   *                [out] Pointers to \c f, \c c and \c h.
   *                On output, <tt>*zero_order_info.c</tt> is updated to \a c(x)
   *                If <tt>this->multi_calc() == true</tt> then
   *                any of the other quantities pointed to in \c zero_order_info may be set on
   *                output, but are not guaranteed to be.
   */
  virtual void imp_calc_c(const Vector& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;
  /** \brief Overridden to compute c_breve(x_breve) and perhaps f(x) and/or h_breve(x_breve)
   *
   * Preconditions:<ul>
   * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
   * <li> <tt>zero_order_info.c != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*zero_order_info.c</tt> is updated to c_breve(x_breve).
   * </ul>
   *
   * @param x       [in]  Unknown vector (size n).
   * @param newx    [in]  True if is a new point.
   * @param zero_order_info_breve
   *                [out] Pointers to \c f, \c c_breve and \c h_breve.
   *                On output, <tt>*zero_order_info.c</tt> is updated to \a c_breve(x_breve)
   *
   * The default implementation calls <tt>this->imp_calc_c()</tt>.
   */
  virtual void imp_calc_c_breve(const Vector& x, bool newx, const ZeroOrderInfo& zero_order_info_breve) const;
  /** \brief Overridden to compute h_breve(x_breve) and perhaps f(x) and/or c_breve(x_breve).
   *
   * Preconditions:<ul>
   * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
   * <li> <tt>zero_order_info.h != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>*zero_order_info.h</tt> is updated to <tt>h_breve(x_breve)</tt>.
   * </ul>
   *
   * @param x       [in]  Unknown vector (size n).
   * @param newx    [in]  True if is a new point.
   * @param zero_order_info_breve
   *                [out] Pointers to \c f, \c c_breve and \c h_breve.
   *                On output, <tt>*zero_order_info.h</tt> is updated to \a h_breve(x_breve)
   *
   * The default implementation throws an exception.
   */
  virtual void imp_calc_h_breve(const Vector& x, bool newx, const ZeroOrderInfo& zero_order_info_breve) const;

  //@}

  /// Assert referece has been set for a quanity
  template<class T>
  void assert_ref_set(T* p, std::string info) const {
    StandardCompositionRelationshipsPack::assert_role_name_set(p, false, info);
  }

private:

  // ////////////////////////////////////////
  // Private data members

#ifdef DOXYGEN_COMPILE
  AbstractLinAlgPack::VectorSpace *space_x;
  AbstractLinAlgPack::VectorSpace *space_c;
  AbstractLinAlgPack::VectorSpace *space_c_breve;
  AbstractLinAlgPack::VectorSpace *space_h_breve;
  Permutation                     *P_var;
  Permtuation                     *P_equ;
#else
  Teuchos::RCP<Permutation>  P_var_;
  Teuchos::RCP<Permutation>  P_equ_;
#endif
  mutable ZeroOrderInfo           first_order_info_;
  mutable ZeroOrderInfo           first_order_info_breve_;
  mutable size_type				num_f_evals_;
  mutable size_type				num_c_evals_;
  
};	// end class NLP

// /////////////////
// Inline members

inline
const NLP::ZeroOrderInfo NLP::zero_order_info() const
{
  return first_order_info_;
}

inline
const NLP::ZeroOrderInfo NLP::zero_order_info_breve() const
{
  return first_order_info_breve_;
}

}	// end namespace NLPInterfacePack 

#endif // NLP_H
