//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
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

#include "ExampleApplication1Dfem.hpp"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#  include "mpi.h"
#else
#  include "Epetra_SerialComm.h"
#endif // HAVE_MPI
#include "Epetra_CrsMatrix.h"

#ifdef EXAMPLEAPPLICATION_DEBUG
#include <iostream>
#endif

#include "Teuchos_dyn_cast.hpp"

ExampleApplication1Dfem::ExampleApplication1Dfem(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params)
{
  initialize(epetra_comm_ptr, params);
}

void ExampleApplication1Dfem::initialize(Teuchos::RCP<Epetra_Comm> &epetra_comm_ptr, Teuchos::ParameterList &params)
{
  epetra_comm_ptr_ = epetra_comm_ptr;
  numElements_ = params.get<int>( "NumElements" );
  
  // This object is derived from NOX::Epetra::Interface
  problemInterfacePtr_ = Teuchos::rcp(new TransientInterface(numElements_, *epetra_comm_ptr_, -20.0, 20.0));

  // Set the PDE nonlinear coefficient for this problem
  problemInterfacePtr_->setPDEfactor(1.0);

  // This is needed to extract the Epetra_Map for the solution
  epetra_map_ptr_ = Teuchos::rcp( new Epetra_Map( problemInterfacePtr_->getMap() ) );
//  Epetra_Vector& soln = problemInterfacePtr_->getSolution();
//  const Epetra_BlockMap& solnBlockMap = soln.Map();
//  std::cout << "typeName(solnBlockMap).name() = " << typeid(solnBlockMap) << std::endl;
//  const Epetra_Map& solnMap = Teuchos::dyn_cast<const Epetra_Map>(solnBlockMap);
//  epetra_map_ptr_ = Teuchos::rcp( new Epetra_Map( solnMap ) );

  // This is needed to extract the Epetra_CrsGraph for the Jacobian
//  Epetra_CrsMatrix &jacobian = problemInterfacePtr_->getJacobian();
//  W_graph_ = Teuchos::rcp(new Epetra_CrsGraph( jacobian.Graph() ) );
  W_graph_ = Teuchos::rcp(new Epetra_CrsGraph( problemInterfacePtr_->getGraph() ) );

}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
ExampleApplication1Dfem::get_x_map() const
{
  return epetra_map_ptr_;
}

Teuchos::RCP<const Epetra_Map>
ExampleApplication1Dfem::get_f_map() const
{
  return epetra_map_ptr_;
}

Teuchos::RCP<const Epetra_Vector>
ExampleApplication1Dfem::get_x_init() const
{
  Epetra_Vector& soln = problemInterfacePtr_->getSolution();
  Teuchos::RCP<Epetra_Vector> x_init = Teuchos::rcp(new Epetra_Vector(soln));
  return x_init;
}

Teuchos::RCP<const Epetra_Vector>
ExampleApplication1Dfem::get_x_dot_init() const
{
  Epetra_Vector& soln = problemInterfacePtr_->getSolution();
  Teuchos::RCP<Epetra_Vector> x_dot_init = Teuchos::rcp(new Epetra_Vector(soln));
  x_dot_init->PutScalar(0.0);
  return x_dot_init;
}

Teuchos::RCP<Epetra_Operator>
ExampleApplication1Dfem::create_W() const
{
  Teuchos::RCP<Epetra_Operator> W = Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));
  return W;
}

EpetraExt::ModelEvaluator::InArgs
ExampleApplication1Dfem::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_x_dot,true);
  inArgs.setSupports(IN_ARG_alpha,true);
  inArgs.setSupports(IN_ARG_beta,true);
  inArgs.setSupports(IN_ARG_t,true);
 // 2007/06/08: rabartl: We have to accept t even if we don't use it!  We have
 // to require that all transient problems define t even if the model does not
 // use it since this is the only way to completely pass an initial condition.
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
ExampleApplication1Dfem::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  return outArgs;
}

void ExampleApplication1Dfem::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  Teuchos::RCP<const Epetra_Vector> x = inArgs.get_x();
  Teuchos::RCP<const Epetra_Vector> xdot = inArgs.get_x_dot();
#ifdef EXAMPLEAPPLICATION_DEBUG
  std::cout << "ExampleApplication1Dfem::evalModel ---------------------------{" << std::endl;
  std::cout << "x = " << std::endl;
  x->Print(std::cout);
  std::cout << "xdot = " << std::endl;
  xdot->Print(std::cout);
#endif // EXAMPLEAPPLICATION_DEBUG
  Teuchos::RCP<Epetra_Vector> f;
  if( (f = outArgs.get_f()).get() ) 
  {
    NOX::Epetra::Interface::Required::FillType flag = NOX::Epetra::Interface::Required::Residual;
    problemInterfacePtr_->evaluate(flag,&*x,&*xdot,0.0,0.0,&*f,NULL);
#ifdef EXAMPLEAPPLICATION_DEBUG
    std::cout << "f = " << std::endl;
    f->Print(std::cout);
#endif // EXAMPLEAPPLICATION_DEBUG
  }
  Teuchos::RCP<Epetra_Operator> W;
  if( (W = outArgs.get_W()).get() ) 
  {
    const double alpha = inArgs.get_alpha();
    const double beta = inArgs.get_beta();
    NOX::Epetra::Interface::Required::FillType flag = NOX::Epetra::Interface::Required::Jac;
//    Epetra_CrsMatrix& jacobian = problemInterfacePtr_->getJacobian();
    Epetra_CrsMatrix& jac = Teuchos::dyn_cast<Epetra_CrsMatrix>(*W);
    problemInterfacePtr_->evaluate(flag,&*x,&*xdot,alpha,beta,NULL,&jac);
#ifdef EXAMPLEAPPLICATION_DEBUG
    std::cout << "jac = " << std::endl;
    jac.Print(std::cout);
#endif // EXAMPLEAPPLICATION_DEBUG
  }
#ifdef EXAMPLEAPPLICATION_DEBUG
  std::cout << "ExampleApplication1Dfem::evalModel ---------------------------}" << std::endl;
#endif // EXAMPLEAPPLICATION_DEBUG
}

Teuchos::RCP<Epetra_Vector> ExampleApplication1Dfem::get_exact_solution( double t ) const
{
  Epetra_Vector& x_exact = problemInterfacePtr_->getExactSoln(t);
  Teuchos::RCP<Epetra_Vector> x_exact_ptr = Teuchos::rcp(&x_exact,false);
  return(x_exact_ptr);
}
