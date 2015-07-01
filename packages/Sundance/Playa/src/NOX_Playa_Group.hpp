/* @HEADER@ */
// ************************************************************************
// 
//                 Playa: Programmable Linear Algebra
//                 Copyright 2012 Sandia Corporation
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


// $Id$ 
// $Source$ 


//   


#ifndef NOX_Playa_GROUP_H
#define NOX_Playa_GROUP_H

#include "NOX_Abstract_Group.H"	// base class
#include "Teuchos_ParameterList.hpp"	// base class

#include "NOX_Common.H"             // class data element (string)
#include "NOX_Playa_Vector.hpp"	    // class data element
#include "Teuchos_RCP.hpp"
#include "PlayaVectorType.hpp"
#include "PlayaVectorSpaceDecl.hpp"
#include "PlayaLinearSolverDecl.hpp"
#include "PlayaLinearOperatorDecl.hpp"
#include "PlayaNonlinearOperator.hpp" // nonlinear operator

// Forward declares
namespace Teuchos 
{
namespace Parameter 
{
class List;
}
}

namespace NOX {
namespace NOXPlaya {

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

class Group : public virtual NOX::Abstract::Group
{

public:

  /*! \brief Constructor.
   *
   * Construct a group given an initial condition, the nonlinear operator that 
   * describes the problem to be solved, and the linear solver
   */
  Group(const Playa::Vector<double>& initcond, 
    const Playa::NonlinearOperator<double>& nonlinOp,
    const Playa::LinearSolver<double>& solver);

  /*! \brief Constructor.
   *
   * Construct a group given an initial condition and the nonlinear operator that 
   * describes the problem to be solved.
   */
  Group(const Playa::NonlinearOperator<double>& nonlinOp,
    const Playa::LinearSolver<double>& solver);

  /*! \brief Constructor.
   *
   * Construct a group given an initial condition, the nonlinear operator that 
   * describes the problem to be solved, the linear solver, and user-specified precision.
   */
  Group(const Playa::Vector<double>& initcond, 
    const Playa::NonlinearOperator<double>& nonlinOp,
    const Playa::LinearSolver<double>& solver,
    int numdigits);

  /*! \brief Constructor.
   *
   * Construct a group given an initial condition, the nonlinear operator that 
   * describes the problem to be solved, and user-specified precision.
   */
  Group(const Playa::NonlinearOperator<double>& nonlinOp,
    const Playa::LinearSolver<double>& solver,
    int numdigits);




  /*! \brief Copy constructor
   *
   * Construct a new group given an existing group to copy from.
   */
  Group(const NOX::NOXPlaya::Group& source, NOX::CopyType type = DeepCopy);

  //! Destructor.
  ~Group();


  NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
  //! See above.
  NOX::Abstract::Group& operator=(const NOX::NOXPlaya::Group& source);

  /** @name "Compute" functions. */
  //@{

  void setX(const NOX::Abstract::Vector& y);
  //! See above
  void setX(const NOX::NOXPlaya::Vector& y);

  void computeX(const NOX::Abstract::Group& grp, 
    const NOX::Abstract::Vector& d, 
    double step);
  //! See above.
  void computeX(const NOX::NOXPlaya::Group& grp, 
    const NOX::NOXPlaya::Vector& d, 
    double step);

  NOX::Abstract::Group::ReturnType computeF();

  NOX::Abstract::Group::ReturnType computeJacobian();

  NOX::Abstract::Group::ReturnType computeGradient();

  NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);

  //@}

  /** @name Jacobian operations.
   *
   * Operations using the Jacobian matrix. These may not be defined in
   * matrix-free scenarios. */

  //@{
  
  NOX::Abstract::Group::ReturnType 
  applyJacobian(const NOX::NOXPlaya::Vector& input, 
    NOX::NOXPlaya::Vector& result) const;

  //! See above
  NOX::Abstract::Group::ReturnType 
  applyJacobian(const NOX::Abstract::Vector& input, 
    NOX::Abstract::Vector& result) const;

  NOX::Abstract::Group::ReturnType 
  applyJacobianTranspose(const NOX::NOXPlaya::Vector& input, 
    NOX::NOXPlaya::Vector& result) const;

  //! See above
  NOX::Abstract::Group::ReturnType 
  applyJacobianTranspose(const NOX::Abstract::Vector& input, 
    NOX::Abstract::Vector& result) const;

  NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(Teuchos::ParameterList& params, 
    const NOX::NOXPlaya::Vector& input, 
    NOX::NOXPlaya::Vector& result) const;

  NOX::Abstract::Group::ReturnType 
  applyJacobianInverse(Teuchos::ParameterList& params, 
    const NOX::Abstract::Vector& input, 
    NOX::Abstract::Vector& result) const;
  
  //@}

  /** @name "Is" functions
   *
   * Checks to see if various objects have been computed. Returns true
   * if the corresponding "compute" function has been called since the
   * last update to the solution vector (via instantiation or
   * computeX). */
  //@{

  bool isF() const;
  bool isJacobian() const;
  bool isGradient() const;
  bool isNewton() const;

  //@}

  /** @name "Get" functions 
   *
   * Note that these function do not check whether or not the vectors
   * are valid. Must use the "Is" functions for that purpose. */
  //@{

  const NOX::Abstract::Vector& getX() const;

  const NOX::Abstract::Vector& getF() const;
  
  double getNormF() const;

  const NOX::Abstract::Vector& getGradient() const;

  const NOX::Abstract::Vector& getNewton() const;

  //! Return RCP to solution vector.  
  virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const 
    {return rcp_dynamic_cast<const NOX::Abstract::Vector>(xVector);}

  //! Return RCP to F(x)
  virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const 
    {return rcp_dynamic_cast<const NOX::Abstract::Vector>(fVector);}



  //! Return RCP to gradient.
  virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const
    {return rcp_dynamic_cast<const NOX::Abstract::Vector>(gradientVector);}

  //! Return RCP to Newton direction.
  virtual Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const 
    {return rcp_dynamic_cast<const NOX::Abstract::Vector>(newtonVector);}
  //@}

#ifdef TRILINOS_6
  virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const;
#else
  virtual RCP<NOX::Abstract::Group> clone(NOX::CopyType type = NOX::DeepCopy) const;
#endif

  //! Print out the group
  void print() const;

protected:

  //! resets all isValid flags to false
  void resetIsValid();

protected:

  //! user-specified precision
  int precision; 
      
  /** @name Vectors */
  //@{
  //! Solution vector.
  RCP<NOX::NOXPlaya::Vector> xVector;
  //! Right-hand-side vector (function evaluation).
  RCP<NOX::NOXPlaya::Vector> fVector;
  //! Newton direction vector.
  RCP<NOX::NOXPlaya::Vector> newtonVector;
  //! Gradient vector (steepest descent vector).
  RCP<NOX::NOXPlaya::Vector> gradientVector;
  //@}

  //! Linear solver that will be used to solve J*step = resid
  mutable Playa::LinearSolver<double> solver;

  //! Linear solver that will be used to solve J*step = resid
  Playa::LinearOperator<double> jacobian;


  //! Problem interface: reference to nonlinear operator passed to group
  Playa::NonlinearOperator<double> nonlinearOp;

  /** @name IsValid flags 
   *  
   * True if the current solution is up-to-date with respect to the
   * currect xVector. */
  //@{
  bool isValidF;
  bool isValidJacobian;
  bool isValidGradient;
  bool isValidNewton;
  //@}
  
  //! Norm of F
  double normF;

};

} // namespace Playa
} // namespace NOX


#endif
