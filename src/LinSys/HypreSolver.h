//******************************************************************************
/*!
  \file      src/LinSys/HypreSolver.h
  \author    J. Bakosi
  \date      Wed 26 Aug 2015 02:28:50 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Hypre solver class
  \details   Hypre solver class.
*/
//******************************************************************************
#ifndef HypreSolver_h
#define HypreSolver_h

#include <iostream>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include <HYPRE.h>
#include <HYPRE_krylov.h>

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace tk {
namespace hypre {

//! Hypre matrix class
class HypreSolver {

  public:
    //! Create and initialize Hypre solver
    void create( ) {
      // Create Hypre solver
      HYPRE_ParCSRPCGCreate( MPI_COMM_WORLD, &m_solver );
      // Set solver parameters, see Hypre manual for more
      HYPRE_PCGSetMaxIter( m_solver, 1000 ); // max iterations
      HYPRE_PCGSetTol( m_solver, 1.0e-8 );   // conv. tolerance
      HYPRE_PCGSetTwoNorm( m_solver, 1);     // use 2-norm as stopping criteria
      HYPRE_PCGSetPrintLevel( m_solver, 1 ); // print out iteration info
      HYPRE_PCGSetLogging( m_solver, 1 );    // for run info
    }

    //! Solve the linear system
    void solve( const HypreMatrix& A,
                const HypreVector& b,
                const HypreVector& x ) const
    {
      HYPRE_ParCSRPCGSetup( m_solver, A.get(), b.get(), x.get() );
      HYPRE_ParCSRPCGSolve( m_solver, A.get(), b.get(), x.get() );
      if (CkMyPe() == 0) {
        int niter;
        HYPRE_PCGGetNumIterations( m_solver, &niter );
        double resnorm;
        HYPRE_PCGGetFinalRelativeResidualNorm( m_solver, &resnorm );
        //std::cout << "it = " << niter << ", norm = " << resnorm << std::endl;
      }
    }

    //! \brief Destructor: destroy Hypre solver
    ~HypreSolver() noexcept { HYPRE_ParCSRPCGDestroy( m_solver ); }

  private:
    HYPRE_Solver m_solver;      //!< Hypre solver
};

} // hypre::
} // tk::

#endif // HypreSolver_h
