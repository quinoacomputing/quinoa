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


#include "NOX_Common.H"
#include "NOX_Playa_Group.hpp"	// class definition
#include "PlayaMPIComm.hpp"
#include "PlayaOut.hpp"
#include "PlayaTabs.hpp"

using Playa::Out;
using Playa::Tabs;
using Teuchos::RCP;

using std::ostream;
using std::cout;
using std::runtime_error;



NOX::NOXPlaya::Group::Group(const Playa::Vector<double>& initcond, 
  const Playa::NonlinearOperator<double>& nonlinOp,
  const Playa::LinearSolver<double>& linsolver) 
  :
  precision(3),  // 3 digits of accuracy is default
  xVector(rcp(new NOX::NOXPlaya::Vector(initcond, precision, DeepCopy))),
  fVector(rcp(new NOX::NOXPlaya::Vector(initcond, precision, ShapeCopy))),
  newtonVector(rcp(new NOX::NOXPlaya::Vector(initcond, precision, ShapeCopy))),
  gradientVector(rcp(new NOX::NOXPlaya::Vector(initcond, precision, ShapeCopy))),
  solver(linsolver),
  jacobian(),
  nonlinearOp(nonlinOp),
  normF(0.0)
{  
  nonlinearOp.setEvalPt(xVector->getPlayaVector());
  resetIsValid();
}

NOX::NOXPlaya::Group::Group(const Playa::NonlinearOperator<double>& nonlinOp,
  const Playa::LinearSolver<double>& linsolver) 
  :
  precision(3), // 3 digits of accuracy is default
  xVector(rcp(new NOX::NOXPlaya::Vector(nonlinOp.getInitialGuess(), precision, DeepCopy))),
  fVector(rcp(new NOX::NOXPlaya::Vector(nonlinOp.getInitialGuess(), precision, ShapeCopy))),
  newtonVector(rcp(new NOX::NOXPlaya::Vector(nonlinOp.getInitialGuess(), precision, ShapeCopy))),
  gradientVector(rcp(new NOX::NOXPlaya::Vector(nonlinOp.getInitialGuess(), precision, ShapeCopy))),
  solver(linsolver),
  jacobian(),
  nonlinearOp(nonlinOp),
  normF(0.0)
{  
  nonlinearOp.setEvalPt(xVector->getPlayaVector());
  resetIsValid();
}

NOX::NOXPlaya::Group::Group(const Playa::Vector<double>& initcond, 
  const Playa::NonlinearOperator<double>& nonlinOp,
  const Playa::LinearSolver<double>& linsolver,
  int numdigits) 
  :
  precision(numdigits),
  xVector(rcp(new NOX::NOXPlaya::Vector(initcond, precision, DeepCopy))),
  fVector(rcp(new NOX::NOXPlaya::Vector(initcond, precision, ShapeCopy))),
  newtonVector(rcp(new NOX::NOXPlaya::Vector(initcond, precision, ShapeCopy))),
  gradientVector(rcp(new NOX::NOXPlaya::Vector(initcond, precision, ShapeCopy))),
  solver(linsolver),
  jacobian(),
  nonlinearOp(nonlinOp),
  normF(0.0)
{  
  nonlinearOp.setEvalPt(xVector->getPlayaVector());
  resetIsValid();
}

NOX::NOXPlaya::Group::Group(const Playa::NonlinearOperator<double>& nonlinOp,
  const Playa::LinearSolver<double>& linsolver,
  int numdigits) 
  :
  precision(numdigits),
  xVector(rcp(new NOX::NOXPlaya::Vector(nonlinOp.getInitialGuess(), precision, DeepCopy))),
  fVector(rcp(new NOX::NOXPlaya::Vector(nonlinOp.getInitialGuess(), precision, ShapeCopy))),
  newtonVector(rcp(new NOX::NOXPlaya::Vector(nonlinOp.getInitialGuess(), precision, ShapeCopy))),
  gradientVector(rcp(new NOX::NOXPlaya::Vector(nonlinOp.getInitialGuess(), precision, ShapeCopy))),
  solver(linsolver),
  jacobian(),
  nonlinearOp(nonlinOp),
  normF(0.0)
{  
  nonlinearOp.setEvalPt(xVector->getPlayaVector());
  resetIsValid();
}


NOX::NOXPlaya::Group::Group(const NOX::NOXPlaya::Group& source, NOX::CopyType type) :
  precision(source.precision),
  xVector(rcp(new NOX::NOXPlaya::Vector(*(source.xVector), precision, type))), 
  fVector(rcp(new NOX::NOXPlaya::Vector(*(source.fVector), precision, type))),  
  newtonVector(rcp(new NOX::NOXPlaya::Vector(*(source.newtonVector), precision, type))),
  gradientVector(rcp(new NOX::NOXPlaya::Vector(*(source.gradientVector), precision, type))),
  solver(source.solver),
  jacobian(source.jacobian),
  nonlinearOp(source.nonlinearOp),
  isValidF(false),
  isValidJacobian(false),
  isValidGradient(false),
  isValidNewton(false),
  normF(0.0)
{
  switch (type) 
  {
    
    case NOX::DeepCopy:
    
      isValidF = source.isValidF;
      isValidGradient = source.isValidGradient;
      isValidNewton = source.isValidNewton;
      isValidJacobian = source.isValidJacobian;
      normF = source.normF;
      break;

    case NOX::ShapeCopy:
      resetIsValid();
      normF = 0.0;
      break;

    default:
      Out::os() << "NOX:Playa::Group - invalid CopyType for copy constructor." << std::endl;
      throw "NOX Playa Error";
  }

}

NOX::NOXPlaya::Group::~Group() 
{
}

void NOX::NOXPlaya::Group::resetIsValid() //private
{
  isValidF = false;
  isValidJacobian = false;
  isValidGradient = false;
  isValidNewton = false;
}


RCP<NOX::Abstract::Group> NOX::NOXPlaya::Group::clone(NOX::CopyType type) const 
{
  return rcp(new NOX::NOXPlaya::Group(*this, type));
}

NOX::Abstract::Group& NOX::NOXPlaya::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const NOX::NOXPlaya::Group&> (source));
}

NOX::Abstract::Group& NOX::NOXPlaya::Group::operator=(const NOX::NOXPlaya::Group& source)
{
  if (this != &source) 
  {

    // Deep Copy of the xVector
    *xVector = *(source.xVector);
    nonlinearOp = source.nonlinearOp;
    solver = source.solver;
    jacobian = source.jacobian;
    precision = source.precision;

    // Update the isValidVectors
    isValidF = source.isValidF;
    isValidGradient = source.isValidGradient;
    isValidNewton = source.isValidNewton;
    isValidJacobian = source.isValidJacobian;
    
    // Only copy vectors that are valid
    if (isValidF) 
    {
      *fVector = *(source.fVector);
      normF = source.normF;
    }

    if (isValidGradient)
      *gradientVector = *(source.gradientVector);
    
    if (isValidNewton)
      *newtonVector = *(source.newtonVector);
    
    if (isValidJacobian)
      jacobian = source.jacobian;
  }
  return *this;
}

void NOX::NOXPlaya::Group::setX(const NOX::Abstract::Vector& y) 
{
  setX(dynamic_cast<const NOX::NOXPlaya::Vector&> (y));
}

void NOX::NOXPlaya::Group::setX(const NOX::NOXPlaya::Vector& y) 
{
  resetIsValid();
  nonlinearOp.setEvalPt(y.getPlayaVector());
  *xVector = y;
}

void NOX::NOXPlaya::Group::computeX(const NOX::Abstract::Group& grp, 
  const NOX::Abstract::Vector& d, 
  double step) 
{
  // Cast to appropriate type, then call the "native" computeX
  const Group& tsfgrp = dynamic_cast<const NOX::NOXPlaya::Group&> (grp);
  const Vector& tsfd = dynamic_cast<const NOX::NOXPlaya::Vector&> (d);
  computeX(tsfgrp, tsfd, step); 
}

void NOX::NOXPlaya::Group::computeX(const Group& grp, const Vector& d, double step) 
{
  Tabs tab;
  resetIsValid();
  xVector->update(1.0, *(grp.xVector), step, d);
}

NOX::Abstract::Group::ReturnType NOX::NOXPlaya::Group::computeF() 
{
 
  if (nonlinearOp.verb() > 2)
  {
    Out::os() << "calling computeF()" << std::endl;
  }

  if (isValidF)
  {
    if (nonlinearOp.verb() > 2)
    {
      Out::os() << "reusing F" << std::endl;
    }
    return NOX::Abstract::Group::Ok;
  }
  else
  {
    if (nonlinearOp.verb() > 2)
    {
      Out::os() << "computing new F" << std::endl;
    }
    nonlinearOp.setEvalPt(xVector->getPlayaVector());
    *fVector = nonlinearOp.getFunctionValue();
    isValidF = true;
    normF = fVector->norm();
  }

  return (isValidF) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType NOX::NOXPlaya::Group::computeJacobian() 
{

  if (nonlinearOp.verb() > 2)
  {
    Out::os() << "calling computeJ()" << std::endl;
  }

  // Skip if the Jacobian is already valid
  if (isValidJacobian)
  {
    return NOX::Abstract::Group::Ok;
  }
  else
  {
    nonlinearOp.setEvalPt(xVector->getPlayaVector());
    jacobian = nonlinearOp.getJacobian();

    isValidJacobian = true;
  }
  return (isValidJacobian) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType NOX::NOXPlaya::Group::computeGradient() 
{
  if (nonlinearOp.verb() > 2)
  {
    Out::os() << "calling computeGrad()" << std::endl;
  }
  if (isValidGradient)
  {
    return NOX::Abstract::Group::Ok;
  }
  else
  {
    if (!isF()) 
    {
      Out::os() << "ERROR: NOX::NOXPlaya::Group::computeGrad() - F is out of date wrt X!" << std::endl;
      return NOX::Abstract::Group::BadDependency;
    }

    if (!isJacobian()) 
    {
      Out::os() << "ERROR: NOX::NOXPlaya::Group::computeGrad() - Jacobian is out of date wrt X!" << std::endl;
      return NOX::Abstract::Group::BadDependency;
    }
  
    // Compute Gradient = J' * F

    NOX::Abstract::Group::ReturnType status 
      = applyJacobianTranspose(*fVector,*gradientVector);
    isValidGradient = (status == NOX::Abstract::Group::Ok);

    // Return result
    return status;
  }
}

NOX::Abstract::Group::ReturnType 
NOX::NOXPlaya::Group::computeNewton(Teuchos::ParameterList& p) 
{
  if (isNewton())
  {
    return NOX::Abstract::Group::Ok;
  }
  else
  {
    if (!isF()) 
    {
      Out::os() << "ERROR: NOX::Example::Group::computeNewton() - invalid F" 
                << std::endl;
      throw "NOX Error";
    }

    if (!isJacobian()) 
    {
      Out::os() << "ERROR: NOX::Example::Group::computeNewton() - invalid Jacobian" << std::endl;
      throw "NOX Error";
    }

/*
  if (p.isParameter("Tolerance"))
  {
  double tol = p.get("Tolerance", tol);
  solver.updateTolerance(tol);
  }
*/

    NOX::Abstract::Group::ReturnType status 
      = applyJacobianInverse(p, *fVector, *newtonVector);
    isValidNewton = (status == NOX::Abstract::Group::Ok);

    
    // Scale soln by -1
    newtonVector->scale(-1.0);

    if (nonlinearOp.verb() > 0)
    {
      Out::os() << "newton step" << std::endl;
      newtonVector->getPlayaVector().print(Out::os());
    }
      
    // Return solution
    return status;
  }
}


NOX::Abstract::Group::ReturnType 
NOX::NOXPlaya::Group::applyJacobian(const Abstract::Vector& input, 
  NOX::Abstract::Vector& result) const
{
  const NOX::NOXPlaya::Vector& tsfinput = dynamic_cast<const NOX::NOXPlaya::Vector&> (input);
  NOX::NOXPlaya::Vector& tsfresult = dynamic_cast<NOX::NOXPlaya::Vector&> (result);
  return applyJacobian(tsfinput, tsfresult);
}

NOX::Abstract::Group::ReturnType 
NOX::NOXPlaya::Group::applyJacobian(const NOX::NOXPlaya::Vector& input, 
  NOX::NOXPlaya::Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
  {
    return NOX::Abstract::Group::BadDependency;
  }
  else
  {
    // Compute result = J * input
    jacobian.apply(input.getPlayaVector(),result.getPlayaVector());
    return NOX::Abstract::Group::Ok;
  }
}

NOX::Abstract::Group::ReturnType 
NOX::NOXPlaya::Group::applyJacobianTranspose(const NOX::Abstract::Vector& input, 
  NOX::Abstract::Vector& result) const
{
  const NOX::NOXPlaya::Vector& tsfinput = dynamic_cast<const NOX::NOXPlaya::Vector&> (input);
  NOX::NOXPlaya::Vector& tsfresult = dynamic_cast<NOX::NOXPlaya::Vector&> (result);
  return applyJacobianTranspose(tsfinput, tsfresult);
}

NOX::Abstract::Group::ReturnType 
NOX::NOXPlaya::Group::applyJacobianTranspose(const NOX::NOXPlaya::Vector& input, NOX::NOXPlaya::Vector& result) const
{
  // Check validity of the Jacobian
  if (!isJacobian()) 
  {
    return NOX::Abstract::Group::BadDependency;
  }
  else
  {
    // Compute result = J^T * input
    jacobian.applyTranspose(input.getPlayaVector(), result.getPlayaVector());
      
    return NOX::Abstract::Group::Ok;
  }
}

NOX::Abstract::Group::ReturnType 
NOX::NOXPlaya::Group::applyJacobianInverse(Teuchos::ParameterList& p, 
  const Abstract::Vector& input, 
  NOX::Abstract::Vector& result) const 
{
  const NOX::NOXPlaya::Vector& tsfinput = dynamic_cast<const NOX::NOXPlaya::Vector&> (input);
  NOX::NOXPlaya::Vector& tsfresult = dynamic_cast<NOX::NOXPlaya::Vector&> (result); 
  return applyJacobianInverse(p, tsfinput, tsfresult);
}

NOX::Abstract::Group::ReturnType 
NOX::NOXPlaya::Group::applyJacobianInverse(Teuchos::ParameterList& p, 
  const NOX::NOXPlaya::Vector& input, 
  NOX::NOXPlaya::Vector& result) const 
{

  if (!isJacobian()) 
  {
    Out::os() << "ERROR: NOX::NOXPlaya::Group::applyJacobianInverse() - invalid Jacobian" << std::endl;
    throw "NOX Error";

  }
/*
  if (p.isParameter("Tolerance"))
  {
  double tol = p.get("Tolerance", tol);
  solver.updateTolerance(tol);
  }
*/
  if (nonlinearOp.verb() > 4)
  {
    Out::os() << "---------------- applying J^-1 ------------------" << std::endl;
    Out::os() << "J=" << std::endl;
    jacobian.print(Out::os());
    Out::os() << "F=" << std::endl;
    input.getPlayaVector().print(Out::os());
  }

  Playa::SolverState<double> status 
    = solver.solve(jacobian, input.getPlayaVector(),
      result.getPlayaVector());
  
  if (status.finalState() != Playa::SolveConverged)
  {
    return NOX::Abstract::Group::Failed;
  }
  else
  {
    if (nonlinearOp.verb() > 2)
    {
      Out::os() << "soln=" << std::endl;
      result.getPlayaVector().print(Out::os());
    }
    return NOX::Abstract::Group::Ok;
  }

}

bool NOX::NOXPlaya::Group::isF() const 
{   
  return isValidF;
}

bool NOX::NOXPlaya::Group::isJacobian() const 
{  
  return isValidJacobian;
}

bool NOX::NOXPlaya::Group::isGradient() const 
{   
  return isValidGradient;
}

bool NOX::NOXPlaya::Group::isNewton() const 
{   
  return isValidNewton;
}

const NOX::Abstract::Vector& NOX::NOXPlaya::Group::getX() const 
{
  return *xVector;
}

const NOX::Abstract::Vector& NOX::NOXPlaya::Group::getF() const 
{  
  if (nonlinearOp.verb() > 2)
  {
    Out::os() << "calling getF()" << std::endl;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(!isF(), runtime_error, 
    "calling getF() with invalid function value");
  return *fVector;
}

double NOX::NOXPlaya::Group::getNormF() const
{
  if (nonlinearOp.verb() > 2)
  {
    Out::os() << "normF = " << normF << std::endl;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(!isF(), runtime_error, 
    "calling normF() with invalid function value");
  return normF;
}

const NOX::Abstract::Vector& NOX::NOXPlaya::Group::getGradient() const 
{ 
  return *gradientVector;
}

const NOX::Abstract::Vector& NOX::NOXPlaya::Group::getNewton() const 
{
  return *newtonVector;
}

void NOX::NOXPlaya::Group::print() const
{
  cout << "x = " << *xVector << "\n";

  if (isValidF) {
    cout << "F(x) = " << *fVector << "\n";
    cout << "|| F(x) || = " << normF << "\n";
  }
  else
    cout << "F(x) has not been computed" << "\n";
  
  cout << std::endl;
}



