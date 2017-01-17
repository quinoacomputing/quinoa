#ifndef __TACHO_EXAMPLE_SOLVER_HPP__
#define __TACHO_EXAMPLE_SOLVER_HPP__

#include "Tacho_Solver.hpp"

#ifdef HAVE_SHYLUTACHO_VTUNE
#include "ittnotify.h"
#endif

namespace Tacho {

#define TACHO_SOLVER_RUN(run, time) timer.reset(); run; const double time = timer.seconds();

  template<typename DeviceSpaceType>
  int exampleSolver(const std::string file_input,
                    const int prunecut,
                    const int max_concurrency,
                    const int nrhs,
                    const int mb,
                    const int nb,
                    const bool verbose) {
    typedef typename
      Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;
    
    const bool detail = false;
    std::cout << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, detail);
    std::cout << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, detail);

    typedef Solver<value_type,ordinal_type,size_type,DeviceSpaceType> SolverType;

#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_pause();
#endif

    int r_val = 0;
    Kokkos::Impl::Timer timer;

    SolverType tacho("Tacho::CholSolver");
    tacho.setPolicy(max_concurrency);
    tacho.setBlocksize(mb, nb);

    ///
    /// Read from matrix market
    ///
    ///     input  - file
    ///     output - AA
    ///
    struct {
      size_type *ap;
      ordinal_type *aj;
      value_type *ax;

      value_type *b;
      bool is_allocated;
    } data;

    timer.reset();
    const auto is_file_null = (file_input.compare("null") == 0);
    if (!is_file_null) {
      data.is_allocated = false;
      std::cout << "Solver:: "
                << "Matrix market input is selected"
                << std::endl;

      // matrix market example
      typename SolverType::CrsMatrixBaseHostType AA("AA");
      std::ifstream in;
      in.open(file_input);
      if (!in.good()) {
        std::cout << "Failed in open the file: " << file_input << std::endl;
        return -1;
      }
      MatrixMarket::read(AA, in);
      
      tacho.setMatrix(AA);
    } else {
      data.is_allocated = true;
      std::cout << "Solver:: "
                << "user data input is selected"
                << std::endl;

      const ordinal_type m = 12, nnz = 46;
      
      data.ap = new size_type[m+1]{
        0, 3, 7, 11, 15, 20, 23, 26, 31, 35, 39, 43, 46
      };
      data.aj = new ordinal_type[nnz]{
        0, 1, 2,
        0, 1, 3, 4,  
        0, 2, 4, 5, 
        1, 3, 6, 7, 
        1, 2, 4, 7, 8, 
        2, 5, 8, 
        3, 6, 9, 
        3, 4, 7, 9, 10, 
        4, 5, 8, 10, 
        6, 7, 9, 11, 
        7, 8, 10, 11, 
        9, 10, 11 
      };
      data.ax = new value_type[nnz]{
        10,  1,  1, 
        1, 10,  1,  2, 
        1, 10,  1,  1, 
        1, 10,  1,  1, 
        2,  1, 10,  2,  2, 
        1, 10,  1, 
        1, 10,  2, 
        1,  2, 10,  2,  1, 
        2,  1, 10,  1, 
        2,  2, 10,  2, 
        1,  1, 10,  1, 
        2, 1, 10 
      };

      typename SolverType::CrsMatrixBaseHostType AA("AA", m, m, nnz);

      // copy user data
      for (auto i=0;i<m;++i) {
        AA.RowPtrBegin(i) = data.ap[i];
        AA.RowPtrEnd(i) = data.ap[i+1];
      }

      for (auto k=0;k<nnz;++k) {
        AA.Col(k) = data.aj[k];
        AA.Value(k) = data.ax[k];
      }

      tacho.setMatrix(AA);
    }
    const double t_input = timer.seconds();    
    std::cout << "Solver:: "
              << "input generation = " << t_input << " [sec] "
              << std::endl;

    if (verbose) {
      std::cout << "- Problem A, B, X - " << std::endl;
      tacho.getA().showMe(std::cout) << std::endl;
    }

    typename SolverType::DenseMatrixBaseHostType BB("BB", tacho.getA().NumRows(), nrhs), XX("XX");
    XX.createConfTo(BB);
    
    {
      const ordinal_type m = BB.NumRows();
      if (!is_file_null) {
        srand(time(NULL));  
        for (auto rhs=0;rhs<nrhs;++rhs) 
          for (auto i=0;i<m;++i) 
            BB.Value(i, rhs) = 1.0; //((value_type)rand()/(RAND_MAX));
      } else {
        data.b = new value_type[m*nrhs];
        
        // init user data 
        srand(time(NULL));        
        for (auto rhs=0;rhs<nrhs;++rhs) 
          for (auto i=0;i<m;++i) 
            data.b[i+m*rhs] = ((value_type)rand()/(RAND_MAX));
        
        // copied to solver
        for (auto rhs=0;rhs<nrhs;++rhs)
          for (auto i=0;i<m;++i)
            BB.Value(i, rhs) = data.b[i+m*rhs];
      } 
    }

    if (verbose) {
      tacho.getB().showMe(std::cout) << std::endl;
      tacho.getX().showMe(std::cout) << std::endl;
    }

    ///
    /// Solver interface
    ///
    TACHO_SOLVER_RUN(tacho.reorder(prunecut), t_reorder);
    TACHO_SOLVER_RUN(tacho.analyze(), t_analyze);

#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_resume();
#endif
    TACHO_SOLVER_RUN(tacho.factorize(), t_factorize);
#ifdef HAVE_SHYLUTACHO_VTUNE
    __itt_pause();
#endif

    tacho.setLeftHandSide(XX);
    tacho.setRightHandSide(BB);

    TACHO_SOLVER_RUN(tacho.solve(), t_solve);

    if (verbose) {
      std::cout << "- Solution X - " << std::endl;
      tacho.getX().showMe(std::cout) << std::endl;
    }

    double norm, error;
    tacho.check(norm, error);
  
    ///
    /// Print out 
    ///
    {
      const auto prec = std::cout.precision();
      std::cout.precision(4);
      
      const auto AA = tacho.getA();
      const auto UU = tacho.getFlatU();
      const auto HU = tacho.getHierU();

      std::cout << std::scientific;
      std::cout << "Solver:: given    matrix = " << AA.NumRows() << " x " << AA.NumCols() << ", nnz = " << AA.NumNonZeros() << std::endl;
      std::cout << "Solver:: factored matrix = " << UU.NumRows() << " x " << UU.NumCols() << ", nnz = " << UU.NumNonZeros() << std::endl;
      std::cout << "Solver:: hier     matrix = " << HU.NumRows() << " x " << HU.NumCols() << ", nnz = " << HU.NumNonZeros() << std::endl;
      
      std::cout << "Solver:: "
                << "input generation = " << t_input << " [sec], "
                << std::endl
                << "Solver:: "
                << "reorder          = " << t_reorder << " [sec] "
                << std::endl
                << "Solver:: "
                << "analyze          = " << t_analyze << " [sec] "
                << std::endl
                << "Solver:: "
                << "factorize        = " << t_factorize << " [sec] "
                << std::endl
                << "Solver:: "
                << "solve            = " << t_solve << " [sec] "
                << std::endl
                << "Solver:: "
                << "norm = " << norm << ", error = " << error << ", rel error = " << (error/norm) 
                << std::endl;

      std::cout << std::endl;
      
      std::cout.unsetf(std::ios::scientific);
      std::cout.precision(prec);
    }

    if (data.is_allocated) {
      delete []data.ap;
      delete []data.aj;
      delete []data.ax;
      delete []data.b;
    }
    
    return r_val;
  }
}

#endif
