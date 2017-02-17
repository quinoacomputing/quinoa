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

#ifndef Rythmos_ExplicitRK_STEPPER_DECL_H
#define Rythmos_ExplicitRK_STEPPER_DECL_H

#include "Rythmos_RKButcherTableauAcceptingStepperBase.hpp"
#include "Rythmos_RKButcherTableauBase.hpp"
#include "Rythmos_Types.hpp"
#include "Thyra_ModelEvaluator.hpp"

#include "Rythmos_StepControlStrategyAcceptingStepperBase.hpp"
#include "Rythmos_StepControlStrategyBase.hpp"

namespace Rythmos {

/** \brief . */
template<class Scalar>
class ExplicitRKStepper : virtual public RKButcherTableauAcceptingStepperBase<Scalar>,
  virtual public StepControlStrategyAcceptingStepperBase<Scalar>
{
  public:
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;
    
    /** \brief . */
    ExplicitRKStepper();

    /** \name Overridden from StepperBase */
    //@{
    
    /** \brief. */
    bool supportsCloning() const;

    /** \brief . */
    RCP<StepperBase<Scalar> > cloneStepperAlgorithm() const;
    
    //@}

    /** \name Overridden from RKButcherTableauAcceptingStepperBase */
    //@{

    /** \brief. */
    void setRKButcherTableau(const RCP<const RKButcherTableauBase<Scalar> > &rkbt);

    /** \brief. */
    RCP<const RKButcherTableauBase<Scalar> > getRKButcherTableau() const;

    //@}

    /** \brief . */
    Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

    /** \brief . */
    void setModel(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);

    /** \brief . */
    void setNonconstModel(const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model);

    /** \brief . */
    Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const;

    /** \brief . */
    RCP<Thyra::ModelEvaluator<Scalar> > getNonconstModel();
    
    /** \brief . */
    ~ExplicitRKStepper();

    /** \brief . */
    void setInitialCondition(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &initialCondition
      );

    /** \brief . */
    Thyra::ModelEvaluatorBase::InArgs<Scalar> getInitialCondition() const;

    /** \brief . */
    Scalar takeStep(Scalar dt, StepSizeType flag);

    /** \brief . */
    const StepStatus<Scalar> getStepStatus() const;

    /** \brief . */
    void describe(
        Teuchos::FancyOStream &out,
        const Teuchos::EVerbosityLevel verbLevel
        ) const;

    /// Redefined from InterpolationBufferBase 
    /// Add points to buffer
    void addPoints(
      const Array<Scalar>& time_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& x_vec
      ,const Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > >& xdot_vec
      );

    /// Get values from buffer
    void getPoints(
      const Array<Scalar>& time_vec
      ,Array<RCP<const VectorBase<Scalar> > >* x_vec
      ,Array<RCP<const VectorBase<Scalar> > >* xdot_vec
      ,Array<ScalarMag>* accuracy_vec) const;

    /** \brief . */
    TimeRange<Scalar> getTimeRange() const;

    /// Get interpolation nodes
    void getNodes(Array<Scalar>* time_vec) const;

    /// Remove interpolation nodes
    void removeNodes(Array<Scalar>& time_vec);

    /// Get order of interpolation
    int getOrder() const;

    /// Redefined from Teuchos::ParameterListAcceptor
    /** \brief . */
    void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
    
    /** \brief. */
    RCP<const Teuchos::ParameterList> getValidParameters() const;


  /** \name Overridden from StepControlStrategyAcceptingStepperBase */
  //@{

  /** \brief . */
  void setStepControlStrategy(
      const RCP<StepControlStrategyBase<Scalar> >& stepControlStrategy
      );

  /** \brief . */
  RCP<StepControlStrategyBase<Scalar> > 
    getNonconstStepControlStrategy();
  
  /** \brief . */
  RCP<const StepControlStrategyBase<Scalar> >
    getStepControlStrategy() const;
  
  //@}

  private:

    Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;
    Teuchos::RCP<Thyra::VectorBase<Scalar> > solution_vector_;
    Teuchos::RCP<Thyra::VectorBase<Scalar> > solution_vector_old_;
    Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > k_vector_;
    Teuchos::RCP<Thyra::VectorBase<Scalar> > ktemp_vector_;
    Teuchos::RCP<Thyra::VectorBase<Scalar> > solution_hat_vector_;

    Thyra::ModelEvaluatorBase::InArgs<Scalar> basePoint_;

    RCP<const RKButcherTableauBase<Scalar> > erkButcherTableau_;

    Scalar t_;
    Scalar t_old_;
    Scalar dt_;
    int numSteps_;
    Scalar LETvalue_;   // ck * e

    Teuchos::RCP<Teuchos::ParameterList> parameterList_;
    RCP<Thyra::VectorBase<Scalar> > ee_; // error (Sidafa) 
    EStepLETStatus stepLETStatus_; // Local Error Test Status (Sidafa)
    TimeRange<Scalar> timeRange_;

    bool isInitialized_;

    bool haveInitialCondition_;

    // Private member functions:
    void defaultInitializeAll_();
    void initialize_();

  // Sidafa 9/4/15
  int rkNewtonConvergenceStatus_;
  RCP<Rythmos::StepControlStrategyBase<Scalar> > stepControl_;
  bool isVariableStep_ = false;

  Scalar takeVariableStep_(Scalar dt, StepSizeType flag);

  Scalar takeFixedStep_(Scalar dt, StepSizeType flag);

};

// Non-member constructors
template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper();

template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model 
    );

template<class Scalar>
RCP<ExplicitRKStepper<Scalar> > explicitRKStepper(
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> >& model,
    const RCP<const RKButcherTableauBase<Scalar> >& rkbt 
    );

} // namespace Rythmos

#endif //Rythmos_ExplicitRK_STEPPER_DECL_H

