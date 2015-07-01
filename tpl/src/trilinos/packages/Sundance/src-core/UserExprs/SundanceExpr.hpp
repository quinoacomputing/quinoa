/* @HEADER@ */
// ************************************************************************
// 
//                             Sundance
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Kevin Long (kevin.long@ttu.edu)
// 

/* @HEADER@ */

#ifndef SUNDANCE_EXPR_H
#define SUNDANCE_EXPR_H

#include "SundanceDefs.hpp"
#include "SundanceExprBase.hpp"
#include "SundanceFunctionIdentifier.hpp"
#include "SundanceMap.hpp"
#include "PlayaHandle.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include <complex>


namespace Sundance
{

class ScalarExpr;
class FuncElementBase;

using namespace Sundance;
/**
 * User-level expression class. Expr is a handle to a
 * reference-counted pointer to a ExprBase subtype. As such,
 * expression copies and assignments are shallow.
 *
 * <h2> Lists </h2>
 *
 * Expressions can be grouped into lists with an arbitrary structure. Important
 * special cases are scalar, vector, and tensor expressions.
 *
 * <h3> Creating Lists </h3>
 *
 * <ul>
 * <li> Vector expressions
 * \code
 * Expr v = List(vx, vy, vz);
 * \endcode
 * <li> Heterogeneous lists
 * \code
 * // Form vector {vx, vy, vz} 
 * Expr v = List(vx, vy, vz);
 * Expr q = new TestFunction(new Lagrange(1));
 * // Form heterogeneous list {{vx, vy, vz}, q}
 * Expr state = List(v, q);
 * \endcode
 * </ul>   
 *
 * <h3> Probing Lists </h3>
 *
 * <ul>
 * <li> Finding the size of the top level of a list
 * \code
 * // Form vector {vx, vy, vz} 
 * Expr v = List(vx, vy, vz);
 * Expr q = new TestFunction(new Lagrange(1));
 * // Form heterogeneous list {{vx, vy, vz}, q}
 * Expr state = List(v, q);
 * // Check top-level size of state list {{vx, vy, vz}, q}
 * int stateSize = state.size(); // returns 2
 * \endcode
 * <li> Finding the total size of a list
 * \code
 * // Check total size of state list 
 * int totalSize = state.totalSize(); // returns 4
 * \endcode
 * </ul>
 * 
 * <h3> Manipulating Lists </h3>
 * 
 * 
 * 
 * <h2> Arithmetic Operations </h2>
 *
 * <ul>
 * <li> Addition
 * \code
 * Expr sum = x + y;
 * \endcode
 * The operands must have identical list structures.
 *
 * <li> Subtraction
 * \code
 * Expr diff = x - y;
 * \endcode
 * The operands must have identical list structures.
 * 
 * <li> Multiplication
 * \code
 * Expr product = x * y;
 * \endcode
 * The operands must have list
 * structures such that the product can be interpreted as a
 * scalar-vector product or as an inner product between vectors
 * or tensors. The multiplication operator is also used to
 * represent the application of a differential operator.
 *
 * <li> Division
 * \code
 * Expr quotient = x / y;
 * \endcode
 * The denominator must be scalar-valued.
 * </ul>
 *
 * <h2> Expression Subtypes </h2>
 * The user-level expression subtypes are listed below, along with examples of their use.
 * <ul>
 * <li> UnknownFunction - Represents an unknown function in a finite-element problem.
 * Unknown functions can be scalar-valued or vector valued.
 *
 * Example of creation of a scalar-valued unknown function:
 * \code
 * Expr u = new UnknownFunction(new Lagrange(1));
 * \endcode
 *
 * Example of creation of a vector-valued unknown function:
 * \code
 * Expr u = new UnknownFunction(List(new Lagrange(1), new Lagrange(1)));
 * \endcode
 *
 * <li> TestFunction - Represents a test function in a finite-element problem.
 * Test functions can be scalar-valued or vector valued.
 *
 * Example of creation of a scalar-valued test function:
 * \code
 * Expr v = new TestFunction(new Lagrange(1));
 * \endcode
 *
 * Example of creation of a vector-valued test function:
 * \code
 * Expr u = new TestFunction(List(new Lagrange(1), new Lagrange(1)));
 * \endcode
 * 
 * <li> Derivative - Represents a spatial derivative operator. 
 * Spatial derivatives are applied
 * using the multiplication operator. 
 * \code
 * Expr dx = new Derivative(0);
 * Expr convect = (u*dx)*u;
 * \endcode
 * Derivative expressions are scalar valued. However, vector differential operators
 * can be created using the List operator. For example,
 * \code
 * Expr dx = new Derivative(0);
 * Expr dy = new Derivative(1);
 * Expr grad = List(dx, dy);
 * \endcode
 * 
 * <li> CoordExpr - Represents a coordinate functions. 
 * \code
 * Expr x = new CoordExpr(0);
 * Expr y = new CoordExpr(1);
 * Expr r = sqrt(x*x + y*y);
 * \endcode
 * Coordinate expressions are scalar valued.
 *
 * <li> CellDiameterExpr - Represents the diameter of a cell. Cell
 * diameters are often used in stabilized methods. 
 * \code
 * Expr h = new CellDiameterExpr();
 * Expr streamlineDiffusion = eps*h*h*((u*grad)*v)*((u*grad)*T);
 * \endcode
 * Cell diameter expressions are scalar valued.
 * </ul>
 */
class Expr : public Playa::Handle<ExprBase>
{
public:
  /* boilerplate handle ctors */
  HANDLE_CTORS(Expr, ExprBase);

  /** Construct with a constant. Creates a ConstantExpr. */
  Expr(const double& c);

  /** Construct with a complex constant. Creates a ComplexExpr. */
  Expr(const std::complex<double>& c);

  /** Add two expressions. The operands must have identical list
   * structures.
   */
  Expr operator+(const Expr& other) const ;
  /** Subtract two expressions. The operands must have identical
   *  list structures. */
  Expr operator-(const Expr& other) const ;
  /** Multiply two expressions. The operands must have list
   * structures such that the product can be interpreted as a
   * scalar-vector product or as an inner product between vectors
   * or tensors. The multiplication operator is also used to
   * represent the application of a differential operator.
   */
  Expr operator*(const Expr& other) const ;
  /** Divide one expression by another. The right operand must be
      a scalar. */
  Expr operator/(const Expr& other) const ;

  /** Divide an expression by a double. */
  Expr operator/(const double& other) const 
    {return operator*(1.0/other);}

  /** Divide an expression by a complex. */
  Expr operator/(const std::complex<double>& other) const 
    {return operator*(1.0/other);}

  /** Unary minus operator */
  Expr operator-() const ;

  /** Return real part of a complex expression */
  Expr real() const ;

  /** Return imaginary part of a complex expression */
  Expr imag() const ;

  /** Return complex conjugate */
  Expr conj() const ;


  /** List element accessor */
  const Expr& operator[](int i) const ;
      
  /** Number of elements in top level of list */
  int size() const ;

  /** Total number of elements in list. */
  int totalSize() const ;

  /** Append a new element to this list */
  void append(const Expr& expr);

  /** Flatten this list */
  Expr flatten() const ;

  /** Flatten list and spectral structure */
  Expr flattenSpectral() const ;

  /** */
  std::string toString() const ;

  /** */
  XMLObject toXML() const ;

  /** */
  bool isComplex() const ;

  /** */
  bool isSpectral() const;

  /** */
  bool isTestElement() const ;

  /** */
  bool isUnknownElement() const ;

  /** */
  void setParameterValue(const double& value);
  /** */
  double getParameterValue() const ;

  /** Indicate whether the expression is independent of the given 
   * functions */
  bool isIndependentOf(const Expr& u) const ;

  
  /** 
   * Indicate whether the expression is nonlinear 
   * with respect to test functions */
  bool isLinearInTests() const ;

  /** 
   * Indicate whether every term in the expression contains test functions */
  bool everyTermHasTestFunctions() const ;

  /** 
   * Indicate whether the expression contains test functions */
  bool hasTestFunctions() const ;

  /** Indicate whether the expression is linear in the given 
   * functions */
  bool isLinearForm(const Expr& u) const ;

  /** Indicate whether the expression is quadratic in the given 
   * functions */
  bool isQuadraticForm(const Expr& u) const ;

  /** Find the unknown functions appearing in this expression */
  void getUnknowns(Set<int>& unkID, Array<Expr>& unks) const ;

  /** Find the test functions appearing in this expression */
  void getTests(Set<int>& testID, Array<Expr>& tests) const ;


#ifndef DOXYGEN_DEVELOPER_ONLY

     

  /**
   * Turn evaluation caching on
   */
  static bool& evaluationCachingOn() {static bool rtn = true; return rtn;}

  /**
   * Show parentheses around every pair of operands
   */
  static bool& showAllParens() {static bool rtn = false; return rtn;}

  /** Create a new handle for an existing ptr.  */
  static Expr handle(const RCP<ExprBase>& ptr)
    {return Expr(ptr);}

  /** */
  static Time& opTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("Expr symbolic ops"); 
      return *rtn;
    }
  /** */
  static Time& outputTimer() 
    {
      static RCP<Time> rtn 
        = TimeMonitor::getNewTimer("Expr output"); 
      return *rtn;
    }

  /** Comparison operator for ordering expr trees, DO NOT use for comparison
   * of numerical values. */
  bool lessThan(const Expr& other) const ;

  /** Comparison operator for ordering expr trees, DO NOT use for comparison
   * of numerical values. */
  bool operator<(const Expr& other) const ;

  /** Equality test operator for ordering expr trees, DO NOT use for comparison
   * of numerical values. */
  bool sameAs(const Expr& other) const ;

  /** */
  Sundance::Map<Expr, int> getSumTree() const ;



  /** */
  const ScalarExpr* scalarExpr() const ;

  /** */
  const FuncElementBase* funcElement() const ;

private:

      

  /** Add two scalar expressions */
  Expr sum(const Expr& other, int sign) const ;

  /** Multiply two scalar expressions */
  Expr multiply(const Expr& other) const ;

  /** Divide two scalar expressions */
  Expr divide(const Expr& other) const ;

  /** Try transformations of complex addition */
  bool tryAddComplex(const Expr& L, const Expr& R, int sign,
    Expr& rtn) const ;
      
  /** Try transformations of complex multiplication */
  bool tryMultiplyComplex(const Expr& L, const Expr& R, 
    Expr& rtn) const ;




#endif /* DOXYGEN_DEVELOPER_ONLY */
};

/** \relates Expr */
inline std::ostream& operator<<(std::ostream& os, const Expr& e)
{
  if (e.ptr().get()==0) {os << "Expr()"; return os;}
  return e.ptr()->toText(os, false);
}

/** \relates Expr */
inline Expr operator+(const double& a, const Expr& x)
{return Expr(a) + x;}

/** \relates Expr */
inline Expr operator-(const double& a, const Expr& x)
{return Expr(a) - x;}

/** \relates Expr */
inline Expr operator*(const double& a, const Expr& x)
{return Expr(a) * x;}

/** \relates Expr */
inline Expr operator/(const double& a, const Expr& x)
{return Expr(a) / x;}

/** \relates Expr */
inline Expr operator+(const std::complex<double>& a, const Expr& x)
{return Expr(a) + x;}

/** \relates Expr */
inline Expr operator-(const std::complex<double>& a, const Expr& x)
{return Expr(a) - x;}

/** \relates Expr */
inline Expr operator*(const std::complex<double>& a, const Expr& x)
{return Expr(a) * x;}

/** \relates Expr */
inline Expr operator/(const std::complex<double>& a, const Expr& x)
{return Expr(a) / x;}

/** \relates Expr */
inline Expr Re(const Expr& a) {return a.real();}

/** \relates Expr */
inline Expr Im(const Expr& a) {return a.imag();}

/** \relates Expr */
inline Expr conj(const Expr& a) {return a.conj();}

/** \relates Expr */
Expr Complex(const Expr& real, const Expr& imag);

/** \relates Expr */
Expr List(const Expr& a);

/** \relates Expr */
Expr List(const Expr& a, const Expr& b);

/** \relates Expr */
Expr List(const Expr& a, const Expr& b, const Expr& c);

/** \relates Expr */
Expr List(const Expr& a, const Expr& b, const Expr& c,
  const Expr& d);

/** \relates Expr */
Expr List(const Expr& a, const Expr& b, const Expr& c,
  const Expr& d, const Expr& e);

/** \relates Expr */
Expr List(const Expr& a, const Expr& b, const Expr& c,
  const Expr& d, const Expr& e, const Expr& f);

/** \relates Expr */
Expr List(const Expr& a, const Expr& b, const Expr& c,
  const Expr& d, const Expr& e, const Expr& f,
  const Expr& g);


/** \relates Expr */
Expr List(const Expr& a, const Expr& b, const Expr& c,
  const Expr& d, const Expr& e, const Expr& f,
  const Expr& g, const Expr& h);

/** \relates Expr */
Expr List(const Expr& a, const Expr& b, const Expr& c,
  const Expr& d, const Expr& e, const Expr& f,
  const Expr& g, const Expr& h, const Expr& i);

/** \relates Expr */
Expr List(const Expr& a, const Expr& b, const Expr& c,
  const Expr& d, const Expr& e, const Expr& f,
  const Expr& g, const Expr& h, const Expr& i,
  const Expr& j);


  

}

#endif
