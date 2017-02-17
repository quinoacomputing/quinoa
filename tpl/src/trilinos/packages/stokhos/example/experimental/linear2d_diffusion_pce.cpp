// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

// Stokhos Stochastic Galerkin
#include "Stokhos_Epetra.hpp"
#include "Stokhos_Sacado.hpp"

// Class implementing our problem
#include "twoD_diffusion_problem_tpetra.hpp"

// Epetra communicator
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

// solver
#include "Ifpack2_Factory.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosPCETpetraAdapter.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosPseudoBlockGmresSolMgr.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Utilities
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

// Scalar types
#include "linear2d_diffusion_scalar_types.hpp"

// Random field types
enum SG_RF { UNIFORM, LOGNORMAL };
const int num_sg_rf = 2;
const SG_RF sg_rf_values[] = { UNIFORM, LOGNORMAL };
const char *sg_rf_names[] = { "Uniform", "Log-Normal" };

// Krylov methods
enum Krylov_Method { GMRES, CG };
const int num_krylov_method = 2;
const Krylov_Method krylov_method_values[] = { GMRES, CG };
const char *krylov_method_names[] = { "GMRES", "CG" };

// Preconditioning approaches
enum SG_Prec { NONE, MEAN, STOCHASTIC };
const int num_sg_prec = 3;
const SG_Prec sg_prec_values[] = { NONE, MEAN, STOCHASTIC };
const char *sg_prec_names[] = { "None",
				"Mean-Based", 
				"Stochastic" };

// Stochastic division approaches
enum SG_Div { DIRECT, SPD_DIRECT, MEAN_DIV, QUAD };
const int num_sg_div = 4;
const SG_Div sg_div_values[] = { DIRECT, SPD_DIRECT, MEAN_DIV, QUAD };
const char *sg_div_names[] = { "Direct",
			       "SPD-Direct",
			       "Mean-Based", 
			       "Quadrature" };

int main(int argc, char *argv[]) {
  typedef double MeshScalar;
  typedef double BasisScalar;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

// Initialize MPI
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  LocalOrdinal MyPID;

  try {

    // Create a communicator for Epetra objects
    RCP<const Epetra_Comm> globalComm;
#ifdef HAVE_MPI
    globalComm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    globalComm = rcp(new Epetra_SerialComm);
#endif
    MyPID = globalComm->MyPID();

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example runs an interlaced stochastic Galerkin solvers.\n");

    int n = 32;
    CLP.setOption("num_mesh", &n, "Number of mesh points in each direction");

    bool symmetric = false;
    CLP.setOption("symmetric", "unsymmetric", &symmetric, 
		  "Symmetric discretization");

    int num_spatial_procs = -1;
    CLP.setOption("num_spatial_procs", &num_spatial_procs, "Number of spatial processors (set -1 for all available procs)");

    SG_RF randField = UNIFORM;
    CLP.setOption("rand_field", &randField, 
		   num_sg_rf, sg_rf_values, sg_rf_names,
		  "Random field type");

    double mu = 0.2;
    CLP.setOption("mean", &mu, "Mean");

    double s = 0.1;
    CLP.setOption("std_dev", &s, "Standard deviation");

    int num_KL = 2;
    CLP.setOption("num_kl", &num_KL, "Number of KL terms");

    int order = 3;
    CLP.setOption("order", &order, "Polynomial order");

    bool normalize_basis = true;
    CLP.setOption("normalize", "unnormalize", &normalize_basis, 
		  "Normalize PC basis");

    Krylov_Method solver_method = GMRES;
    CLP.setOption("solver_method", &solver_method, 
		  num_krylov_method, krylov_method_values, krylov_method_names, 
		  "Krylov solver method");

    SG_Prec prec_method = MEAN;
    CLP.setOption("prec_method", &prec_method, 
		  num_sg_prec, sg_prec_values, sg_prec_names,
		  "Preconditioner method");

    SG_Div division_method = DIRECT;
    CLP.setOption("division_method", &division_method, 
		  num_sg_div, sg_div_values, sg_div_names,
		  "Stochastic division method");

    double solver_tol = 1e-12;
    CLP.setOption("solver_tol", &solver_tol, "Outer solver tolerance");

    CLP.parse( argc, argv );

    if (MyPID == 0) {
      std::cout << "Summary of command line options:" << std::endl
		<< "\tnum_mesh           = " << n << std::endl
		<< "\tsymmetric          = " << symmetric << std::endl
		<< "\tnum_spatial_procs  = " << num_spatial_procs << std::endl
		<< "\trand_field         = " << sg_rf_names[randField] 
		<< std::endl
		<< "\tmean               = " << mu << std::endl
		<< "\tstd_dev            = " << s << std::endl
		<< "\tnum_kl             = " << num_KL << std::endl
		<< "\torder              = " << order << std::endl
		<< "\tnormalize_basis    = " << normalize_basis << std::endl
		<< "\tsolver_method      = " << krylov_method_names[solver_method] << std::endl
		<< "\tprec_method        = " << sg_prec_names[prec_method] 
		<< std::endl
		<< "\tdivision_method     = " << sg_div_names[division_method] 
		<< std::endl;
    }
    bool nonlinear_expansion = false;
    if (randField == UNIFORM)
      nonlinear_expansion = false;
    else if (randField == LOGNORMAL)
      nonlinear_expansion = true;

    {
    TEUCHOS_FUNC_TIME_MONITOR("Total PCE Calculation Time");

    // Create Stochastic Galerkin basis and expansion
    Teuchos::Array< RCP<const Stokhos::OneDOrthogPolyBasis<LocalOrdinal,BasisScalar> > > bases(num_KL); 
    for (LocalOrdinal i=0; i<num_KL; i++)
      if (randField == UNIFORM)
	bases[i] = rcp(new Stokhos::LegendreBasis<LocalOrdinal,BasisScalar>(
			 order, normalize_basis));
      else if (randField == LOGNORMAL)      
	bases[i] = rcp(new Stokhos::HermiteBasis<int,double>(
			 order, normalize_basis));
    RCP<const Stokhos::CompletePolynomialBasis<LocalOrdinal,BasisScalar> > basis = 
      rcp(new Stokhos::CompletePolynomialBasis<LocalOrdinal,BasisScalar>(bases,
		     1e-12));
    LocalOrdinal sz = basis->size();
    RCP<Stokhos::Sparse3Tensor<LocalOrdinal,BasisScalar> > Cijk = 
      basis->computeTripleProductTensor(sz);
    RCP<const Stokhos::Quadrature<int,double> > quad = 
      rcp(new Stokhos::TensorProductQuadrature<int,double>(basis));
    RCP<ParameterList> expn_params = Teuchos::rcp(new ParameterList);
    if (division_method == MEAN_DIV) {
      expn_params->set("Division Strategy", "Mean-Based");
      expn_params->set("Use Quadrature for Division", false);
    }
    else if (division_method == DIRECT) {
      expn_params->set("Division Strategy", "Dense Direct");
      expn_params->set("Use Quadrature for Division", false);
    }
    else if (division_method == SPD_DIRECT) {
      expn_params->set("Division Strategy", "SPD Dense Direct");
      expn_params->set("Use Quadrature for Division", false);
    }
    else if (division_method == QUAD) {
      expn_params->set("Use Quadrature for Division", true);
    }
    RCP<Stokhos::OrthogPolyExpansion<LocalOrdinal,BasisScalar> > expansion = 
      rcp(new Stokhos::QuadOrthogPolyExpansion<LocalOrdinal,BasisScalar>(
    	    basis, Cijk, quad, expn_params));

    if (MyPID == 0)
      std::cout << "Stochastic Galerkin expansion size = " << sz << std::endl;

    // Create stochastic parallel distribution
    ParameterList parallelParams;
    parallelParams.set("Number of Spatial Processors", num_spatial_procs);
    // parallelParams.set("Rebalance Stochastic Graph", true);
    // Teuchos::ParameterList& isorropia_params = 
    //   parallelParams.sublist("Isorropia");
    // isorropia_params.set("Balance objective", "nonzeros");
    RCP<Stokhos::ParallelData> sg_parallel_data =
      rcp(new Stokhos::ParallelData(basis, Cijk, globalComm, parallelParams));
    RCP<const EpetraExt::MultiComm> sg_comm = 
      sg_parallel_data->getMultiComm();
    RCP<const Epetra_Comm> app_comm = 
      sg_parallel_data->getSpatialComm();

    // Create Teuchos::Comm from Epetra_Comm
    RCP< Teuchos::Comm<int> > teuchos_app_comm;
#ifdef HAVE_MPI
    RCP<const Epetra_MpiComm> app_mpi_comm = 
      Teuchos::rcp_dynamic_cast<const Epetra_MpiComm>(app_comm);
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > raw_mpi_comm = 
      Teuchos::opaqueWrapper(app_mpi_comm->Comm());
    teuchos_app_comm = rcp(new Teuchos::MpiComm<int>(raw_mpi_comm));
#else
    teuchos_app_comm = rcp(new Teuchos::SerialComm<int>());
#endif

    // Create application
    typedef twoD_diffusion_problem<Scalar,MeshScalar,BasisScalar,LocalOrdinal,GlobalOrdinal,Node> problem_type;
    RCP<problem_type> model = 
      rcp(new problem_type(teuchos_app_comm, n, num_KL, s, mu, 
			   nonlinear_expansion, symmetric));

    // Create vectors and operators
    typedef problem_type::Tpetra_Vector Tpetra_Vector;
    typedef problem_type::Tpetra_CrsMatrix Tpetra_CrsMatrix;
    typedef Tpetra::MatrixMarket::Writer<Tpetra_CrsMatrix> Writer;
    RCP<Tpetra_Vector> p = Tpetra::createVector<Scalar>(model->get_p_map(0));
    RCP<Tpetra_Vector> x = Tpetra::createVector<Scalar>(model->get_x_map());
    x->putScalar(0.0);
    RCP<Tpetra_Vector> f = Tpetra::createVector<Scalar>(model->get_f_map());
    RCP<Tpetra_Vector> dx = Tpetra::createVector<Scalar>(model->get_x_map());
    RCP<Tpetra_CrsMatrix> J = model->create_W();
    RCP<Tpetra_CrsMatrix> J0;
    if (prec_method == MEAN)
      J0 = model->create_W();

    // Set PCE expansion of p
    p->putScalar(0.0);
    Teuchos::ArrayRCP<Scalar> p_view = p->get1dViewNonConst();
    for (Teuchos::ArrayRCP<Scalar>::size_type i=0; i<p_view.size(); i++)
      p_view[i].reset(expansion);
    Teuchos::Array<double> point(num_KL, 1.0);
    Teuchos::Array<double> basis_vals(sz);
    basis->evaluateBases(point, basis_vals);
    if (order > 0) {
      for (int i=0; i<num_KL; i++) {
	p_view[i].term(i,1) = 1.0 / basis_vals[i+1];
      }
    }

    // Create preconditioner
    typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
    Teuchos::RCP<Tprec> M;
    if (prec_method != NONE) {
      Teuchos::ParameterList precParams;
      std::string prec_name = "RILUK";
      precParams.set("fact: iluk level-of-fill", 1);
      precParams.set("fact: iluk level-of-overlap", 0);
      Ifpack2::Factory factory;
      if (prec_method == MEAN) 
	M = factory.create<Tpetra_CrsMatrix>(prec_name, J0);
      else if (prec_method == STOCHASTIC)
	M = factory.create<Tpetra_CrsMatrix>(prec_name, J);
      M->setParameters(precParams);
    }

    // Evaluate model
    model->computeResidual(*x, *p, *f);
    model->computeJacobian(*x, *p, *J);

    // Compute mean for mean-based preconditioner
    if (prec_method == MEAN) {
      size_t nrows = J->getNodeNumRows();
      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> values;
      J0->resumeFill();
      for (size_t i=0; i<nrows; i++) {
	J->getLocalRowView(i, indices, values);
	Teuchos::Array<Scalar> values0(values.size());
	for (LocalOrdinal j=0; j<values.size(); j++)
	  values0[j] = values[j].coeff(0);
	J0->replaceLocalValues(i, indices, values0);
      }
      J0->fillComplete(Tpetra::DoOptimizeStorage);
    }

    // compute preconditioner
    if (prec_method != NONE) {
      M->initialize();
      M->compute();
    }

    // Setup Belos solver
    RCP<ParameterList> belosParams = rcp(new ParameterList);
    belosParams->set("Num Blocks", 20);
    belosParams->set("Convergence Tolerance", solver_tol);
    belosParams->set("Maximum Iterations", 1000);
    belosParams->set("Verbosity", 33);
    belosParams->set("Output Style", 1);
    belosParams->set("Output Frequency", 10);
    typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> OP;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Belos::OperatorTraits<double,MV,OP> BOPT;
    typedef Belos::MultiVecTraits<double,MV> BMVT;
    typedef Belos::LinearProblem<double,MV,OP> BLinProb;
    RCP< BLinProb > problem = rcp(new BLinProb(J, dx, f));
    if (prec_method != NONE)
      problem->setRightPrec(M);
    problem->setProblem();
    RCP<Belos::SolverManager<double,MV,OP> > solver;
    if (solver_method == CG)
      solver = rcp(new Belos::PseudoBlockCGSolMgr<double,MV,OP>(problem,
    								belosParams));
    else if (solver_method == GMRES)
      solver = rcp(new Belos::PseudoBlockGmresSolMgr<double,MV,OP>(problem,
    								   belosParams));

    // Print initial residual norm
    std::vector<double> norm_f(1);
    BMVT::MvNorm(*f, norm_f);
    if (MyPID == 0)
      std::cout << "\nInitial residual norm = " << norm_f[0] << std::endl;

    // Solve linear system
    Belos::ReturnType ret = solver->solve();

    if (MyPID == 0) {
      if (ret == Belos::Converged)
	std::cout << "Solver converged!" << std::endl;
      else
	std::cout << "Solver failed to converge!" << std::endl;
    }

    // Update x
    x->update(-1.0, *dx, 1.0);
    Writer::writeDenseFile("stochastic_solution.mm", x);

    // Compute new residual & response function
    RCP<Tpetra_Vector> g = Tpetra::createVector<Scalar>(model->get_g_map(0));
    f->putScalar(0.0);
    model->computeResidual(*x, *p, *f);
    model->computeResponse(*x, *p, *g);

    // Print final residual norm
    BMVT::MvNorm(*f, norm_f);
    if (MyPID == 0)
      std::cout << "\nFinal residual norm = " << norm_f[0] << std::endl;

    // Print response
    std::cout << "\nResponse =      " << std::endl;
    //Writer::writeDense(std::cout, g);
    Writer::writeDenseFile("stochastic_residual.mm", f);

    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

  }
  
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std:: endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

}
