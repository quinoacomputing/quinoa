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
#ifndef Rythmos_IMPLICIT_BDF_CONVERGENCETEST_H
#define Rythmos_IMPLICIT_BDF_CONVERGENCETEST_H

#include "Rythmos_Types.hpp"
#include "Rythmos_ConvergenceTestHelpers.hpp"
#include "Rythmos_ImplicitBDFStepper.hpp"
#include "Rythmos_TimeStepNonlinearSolver.hpp"

namespace Rythmos {

template<class Scalar>
class ImplicitBDFStepperFactory : public virtual StepperFactoryBase<Scalar>
{
  public:
    ImplicitBDFStepperFactory(RCP<ModelFactoryBase<Scalar> > modelFactory) 
    { 
      modelFactory_ = modelFactory;
      order_ = 1; 
    }
    virtual ~ImplicitBDFStepperFactory() {}
    void setOrder(int order) { order_ = order; }
    int maxOrder() { return 2; } // This can change when we ramp-up the step-sizes correctly for higher orders than 2
    RCP<StepperBase<Scalar> > getStepper() const 
    { 
      RCP<ModelEvaluator<Scalar> > model = modelFactory_->getModel();
      Thyra::ModelEvaluatorBase::InArgs<Scalar> model_ic = model->getNominalValues();
      RCP<Rythmos::TimeStepNonlinearSolver<Scalar> >
        nonlinearSolver = Rythmos::timeStepNonlinearSolver<Scalar>();
      RCP<ParameterList> nonlinearSolverPL = Teuchos::parameterList();
      nonlinearSolverPL->get("Default Tol",1.0e-9); // Set default if not set
      nonlinearSolver->setParameterList(nonlinearSolverPL);
      RCP<ImplicitBDFStepper<Scalar> > stepper = Teuchos::rcp(new ImplicitBDFStepper<Scalar>(model,nonlinearSolver));
      RCP<ParameterList> bdfPL = Teuchos::parameterList();
      Teuchos::ParameterList& stepControlPL = bdfPL->sublist("Step Control Settings");
      stepControlPL.set("minOrder",order_);
      stepControlPL.set("maxOrder",order_);
      Teuchos::ParameterList& stepControlVOPL = stepControlPL.sublist("VerboseObject");
      stepControlVOPL.set("Verbosity Level","none");
      stepper->setParameterList(bdfPL);
      stepper->setInitialCondition(model_ic);
      return stepper;
    }
  private:
    RCP<ModelFactoryBase<Scalar> > modelFactory_;
    int order_;
};
// non-member constructor
template<class Scalar>
RCP<ImplicitBDFStepperFactory<Scalar> > implicitBDFStepperFactory(
    RCP<ModelFactoryBase<Scalar> > modelFactory)
{
  RCP<ImplicitBDFStepperFactory<Scalar> > ibdfFactory = Teuchos::rcp(
      new ImplicitBDFStepperFactory<Scalar>(modelFactory)
      );
  return ibdfFactory;
}

} // namespace Rythmos 

#endif // Rythmos_IMPLICIT_BDF_CONVERGENCETEST_H

