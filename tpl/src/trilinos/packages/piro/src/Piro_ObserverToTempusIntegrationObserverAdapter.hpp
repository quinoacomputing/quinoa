// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PIRO_OBSERVERTOTEMPUSINTEGRATIONOBSERVERADAPTER_HPP
#define PIRO_OBSERVERTOTEMPUSINTEGRATIONOBSERVERADAPTER_HPP

#include "Tempus_IntegratorObserver.hpp"

#include "Piro_ObserverBase.hpp"

#include "Teuchos_RCP.hpp"

namespace Piro {

template <typename Scalar>
class ObserverToTempusIntegrationObserverAdapter : public Tempus::IntegratorObserver<Scalar> {
 
 
public:


  // Constructor
  ObserverToTempusIntegrationObserverAdapter(
    const Teuchos::RCP<Tempus::SolutionHistory<Scalar> >& solutionHistory,
    const Teuchos::RCP<Tempus::TimeStepControl<Scalar> >& timeStepControl,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &wrappedObserver); 

  // Overridden from Tempus::IntegratorObserver
  
  //@{
  /// Destructor

  virtual ~ObserverToTempusIntegrationObserverAdapter(); 

  /// Observe the beginning of the time integrator.
  virtual void observeStartIntegrator();

  /// Observe the beginning of the time step loop.
  virtual void observeStartTimeStep();

  /// Observe after the next time step size is selected.
  virtual void observeNextTimeStep(Tempus::Status & integratorStatus);

  /// Observe before Stepper takes step.
  virtual void observeBeforeTakeStep();

  /// Observe after Stepper takes step.
  virtual void observeAfterTakeStep();

  /// Observe after accepting time step.
  virtual void observeAcceptedTimeStep(Tempus::Status & integratorStatus);

  /// Observe the end of the time integrator.
  virtual void observeEndIntegrator(const Tempus::Status integratorStatus);
  //@}

private:

  void observeTimeStep();
  Teuchos::RCP<Tempus::SolutionHistory<Scalar> > solutionHistory_;
  Teuchos::RCP<Tempus::TimeStepControl<Scalar> > timeStepControl_;
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  bool hasSensitivities_; 
  Teuchos::RCP<ObserverBase<Scalar> > wrappedObserver_;
};

} // namespace Piro

#include "Piro_ObserverToTempusIntegrationObserverAdapter_Def.hpp"

#endif /* PIRO_OBSERVERTOTEMPUSINTEGRATIONOBSERVERADAPTER_HPP */
