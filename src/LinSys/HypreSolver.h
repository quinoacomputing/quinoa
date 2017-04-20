// *****************************************************************************
/*!
  \file      src/LinSys/HypreSolver.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Hypre solver class
  \details   Hypre solver class.
*/
// *****************************************************************************
#ifndef HypreSolver_h
#define HypreSolver_h

#include <iostream>

#include "NoWarning/charm.h"

#include <HYPRE.h>
#include "NoWarning/HYPRE_krylov.h"

namespace tk {
namespace hypre {

//! Hypre matrix class
class HypreSolver {

  public:
    #if defined(__clang__)
      #pragma clang diagnostic push
      #pragma clang diagnostic ignored "-Wold-style-cast"
    #endif

    //! Create and initialize Hypre solver
    void create( ) {
      // Create Hypre solver
      HYPRE_ParCSRPCGCreate( MPI_COMM_WORLD, &m_solver );
      // Set solver parameters, see Hypre manual for more
      HYPRE_PCGSetMaxIter( m_solver, 1000 ); // max iterations
      HYPRE_PCGSetTol( m_solver, 1.0e-14 );  // conv. tolerance
      HYPRE_PCGSetTwoNorm( m_solver, 1);     // use 2-norm as stopping criteria
      HYPRE_PCGSetPrintLevel( m_solver, 1 ); // print out iteration info
      HYPRE_PCGSetLogging( m_solver, 1 );    // for run info
    }

    #if defined(__clang__)
      #pragma clang diagnostic pop
    #endif

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
