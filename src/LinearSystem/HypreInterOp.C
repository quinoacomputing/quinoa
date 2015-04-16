//******************************************************************************
/*!
  \file      src/Mesh/HypreInterOp.C
  \author    J. Bakosi
  \date      Thu 16 Apr 2015 05:35:26 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Interoperation with the Hypre library
  \details   Interoperation with the Hypre library, used for linear solvers.
*/
//******************************************************************************
#include <ExceptionMPI.h>
#include <HypreInterOp.h>

#include <HYPRE_parcsr_ls.h>

namespace tk {
namespace hypre {

void test( /*const std::map< int, std::vector< std::size_t > >& contribute*/ )
{
//   int peid;
//   MPI_Comm_rank( MPI_COMM_WORLD, &peid );
// 
//   std::cout << peid << ": ";
//   for (const auto& c : contribute) {
//      std::cout << c.first << "( ";
//      for (auto s : c.second ) std::cout << s << " ";
//      std::cout << " ) ";
//   }
//   std::cout << '\n';
// 
//   std::size_t l = std::numeric_limits< std::size_t >::max();
//   std::size_t u = 0;
//   for (const auto& c : contribute)
//     for (auto s : c.second) {
//       if (s < l) l = s;
//       if (s > u) u = s;
//     }
// 
//   std::cout << peid << ": " << l << "..." << u << '\n';

  HYPRE_IJMatrix A;
  int ilower = 0;//static_cast< int >( l );
  int iupper = 1;//static_cast< int >( u );
  HYPRE_IJMatrixCreate( MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A );
  HYPRE_IJMatrixDestroy( A );
}

} // hypre::
} // tk::
