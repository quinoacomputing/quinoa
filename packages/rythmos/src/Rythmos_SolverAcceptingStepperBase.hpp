//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_SOLVER_ACCEPTING_STEPPER_BASE_HPP
#define RYTHMOS_SOLVER_ACCEPTING_STEPPER_BASE_HPP


#include "Rythmos_StepperBase.hpp"
#include "Thyra_NonlinearSolverBase.hpp"


namespace Rythmos {


/** \brief Mix-in interface all implicit stepper objects that accept a
 * nonlinear solver to be used to compute the timestep.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class SolverAcceptingStepperBase : virtual public StepperBase<Scalar>
{
public:

  /** \brief . */
  virtual void setSolver(
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> > &solver
    ) = 0;

  /** \brief . */
  virtual Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >
  getNonconstSolver() = 0;

  /** \brief . */
  virtual Teuchos::RCP<const Thyra::NonlinearSolverBase<Scalar> >
  getSolver() const = 0;

};


} // namespace Rythmos


#endif // RYTHMOS_SOLVER_ACCEPTING_STEPPER_BASE_HPP
