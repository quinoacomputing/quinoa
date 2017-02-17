/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef HYBRIDFEM_NONLINEAR_HPP
#define HYBRIDFEM_NONLINEAR_HPP

#include <utility>
#include <iostream>
#include <iomanip>

#include <KokkosArray_Array.hpp>
#include <KokkosArray_MultiVector.hpp>
#include <SparseLinearSystem.hpp>
#include <SparseLinearSystemFill.hpp>
#include <FEMesh.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace HybridFEM {
namespace NonLinear {

struct PerformanceData {
  double mesh_time ;
  double graph_time ;
  double elem_time ;
  double matrix_gather_fill_time ;
  double matrix_boundary_condition_time ;
  double cg_iteration_time ;
  size_t cg_iteration_count ;
  size_t newton_iteration_count ;
  double error_max ;

  PerformanceData()
    : mesh_time(0)
    , graph_time(0)
    , elem_time(0)
    , matrix_gather_fill_time(0)
    , matrix_boundary_condition_time(0)
    , cg_iteration_time(0)
    , cg_iteration_count(0)
    , newton_iteration_count(0)
    , error_max(0)
    {}

  void best( const PerformanceData & rhs )
  {
    mesh_time = std::min( mesh_time , rhs.mesh_time );
    graph_time = std::min( graph_time , rhs.graph_time );
    elem_time = std::min( elem_time , rhs.elem_time );
    matrix_gather_fill_time = std::min( matrix_gather_fill_time , rhs.matrix_gather_fill_time );
    matrix_boundary_condition_time = std::min( matrix_boundary_condition_time , rhs.matrix_boundary_condition_time );
    cg_iteration_time = std::min( cg_iteration_time , rhs.cg_iteration_time );
    cg_iteration_count = std::min( cg_iteration_count , rhs.cg_iteration_count );
    newton_iteration_count = std::min( newton_iteration_count , rhs.newton_iteration_count );
    error_max = std::min( error_max , rhs.error_max );
  }
};

//----------------------------------------------------------------------------

template< typename ScalarType , typename ScalarCoordType , class Device >
struct ElementComputation ;

template< typename ScalarType , typename ScalarCoordType , class Device >
struct DirichletSolution ;

template< typename ScalarType , typename ScalarCoordType , class Device >
struct DirichletResidual ;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class ManufacturedSolution {
public:

  // Manufactured solution for one dimensional nonlinear PDE
  //
  //  -K T_zz + T^2 = 0 ; T(zmin) = T_zmin ; T(zmax) = T_zmax
  //
  //  Has an analytic solution of the form:
  //
  //    T(z) = ( a ( z - zmin ) + b )^(-2) where K = 1 / ( 6 a^2 )
  //
  //  Given T_0 and T_L compute K for this analytic solution.
  //
  //  Two analytic solutions:
  //
  //    Solution with singularity:
  //    , a( ( 1.0 / sqrt(T_zmax) + 1.0 / sqrt(T_zmin) ) / ( zmax - zmin ) )
  //    , b( -1.0 / sqrt(T_zmin) )
  //
  //    Solution without singularity:
  //    , a( ( 1.0 / sqrt(T_zmax) - 1.0 / sqrt(T_zmin) ) / ( zmax - zmin ) )
  //    , b( 1.0 / sqrt(T_zmin) )

  const double zmin ;
  const double zmax ;
  const double T_zmin ;
  const double T_zmax ;
  const double a ;
  const double b ;
  const double K ;

  ManufacturedSolution( const double arg_zmin ,
                        const double arg_zmax ,
                        const double arg_T_zmin ,
                        const double arg_T_zmax )
    : zmin( arg_zmin )
    , zmax( arg_zmax )
    , T_zmin( arg_T_zmin )
    , T_zmax( arg_T_zmax )
    , a( ( 1.0 / sqrt(T_zmax) - 1.0 / sqrt(T_zmin) ) / ( zmax - zmin ) )
    , b( 1.0 / sqrt(T_zmin) )
    , K( 1.0 / ( 6.0 * a * a ) )
    {}

  double operator()( const double z ) const
  {
    const double tmp = a * ( z - zmin ) + b ;
    return 1.0 / ( tmp * tmp );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename Scalar , class Device >
PerformanceData run( comm::Machine machine ,
                     const int global_count_x ,
                     const int global_count_y ,
                     const int global_count_z ,
                     const bool print_error )
{
  typedef Device                           device_type;
  typedef typename device_type::size_type  size_type ;
  typedef Scalar                           scalar_type ;

  //------------------------------------
  // The amount of nonlinearity is proportional to the ratio
  // between T(zmax) and T(zmin).  For the manufactured solution
  // 0 < T(zmin) and 0 < T(zmax)

  const ManufacturedSolution 
   exact_solution( /* zmin */ 0 ,
                   /* zmax */ global_count_z - 1 , 
                   /* T(zmin) */ 1 ,
                   /* T(zmax) */ 20 );

  //-----------------------------------
  // Convergence Criteria and perf data:

  const size_t cg_iteration_limit = 200 ;
  const double cg_tolerance = 1e-14 ;

  const size_t newton_iteration_limit = 150 ;
  const double newton_tolerance = 1e-14 ;

  size_t cg_iteration_count_total = 0 ;
  double cg_iteration_time = 0 ;

  size_t newton_iteration_count = 0 ;
  double residual_norm_init = 0 ;
  double residual_norm = 0 ;

  PerformanceData perf_data ;

  //------------------------------------
  // FEMesh types:

  enum { ElementNodeCount = 8 };

  typedef double coordinate_scalar_type ;

  typedef FEMesh< coordinate_scalar_type ,
                  ElementNodeCount , device_type > mesh_type ;

  //------------------------------------
  // Sparse linear system types:

  typedef KokkosArray::MultiVector< Scalar , Device >   vector_type ;
  typedef KokkosArray::CrsMatrix< Scalar , Device >     matrix_type ;
  typedef typename matrix_type::graph_type         matrix_graph_type ;
  typedef typename matrix_type::coefficients_type  matrix_coefficients_type ;

  typedef KokkosArray::Impl::Factory< matrix_graph_type , mesh_type > graph_factory ;

  //------------------------------------
  // Problem setup types:

  typedef ElementComputation < Scalar , Scalar , device_type >  ElementFunctor ;
  typedef DirichletSolution  < Scalar , Scalar , device_type > DirichletSolutionFunctor ;
  typedef DirichletResidual  < Scalar , Scalar , device_type > DirichletResidualFunctor ;

  //------------------------------------

  matrix_type jacobian ;
  vector_type residual ;
  vector_type delta ;
  vector_type nodal_solution ;

  typename graph_factory::element_map_type element_map ;

  //------------------------------------
  // Generate mesh and corresponding sparse matrix graph

  KokkosArray::Impl::Timer wall_clock ;

  mesh_type mesh =
    box_mesh_fixture< coordinate_scalar_type , device_type >
      ( comm::size( machine ) , comm::rank( machine ) ,
        global_count_x , global_count_y , global_count_z );

  mesh.parallel_data_map.machine = machine ;

  device_type::fence();
  perf_data.mesh_time = comm::max( machine , wall_clock.seconds() );

  const size_t element_count = mesh.elem_node_ids.dimension(0);

  //------------------------------------
  // Generate sparse matrix graph and element->graph map.

  wall_clock.reset();

  graph_factory::create( mesh , jacobian.graph , element_map );

  device_type::fence();
  perf_data.graph_time = comm::max( machine , wall_clock.seconds() );

  //------------------------------------
  // Allocate linear system coefficients and rhs:

  const size_t local_owned_length = jacobian.graph.row_map.length();
  const size_t local_total_length = mesh.node_coords.dimension(0);

  jacobian.coefficients =
    KokkosArray::create_multivector< matrix_coefficients_type >( jacobian.graph.entries.dimension(0) );

  // Nonlinear residual for owned nodes:
  residual =
    KokkosArray::create_multivector< vector_type >( local_owned_length );

  // Nonlinear solution for owned and ghosted nodes:
  nodal_solution = 
    KokkosArray::create_multivector< vector_type >( local_total_length );

  // Nonlinear solution update for owned nodes:
  delta =
    KokkosArray::create_multivector< vector_type >( local_owned_length );

  //------------------------------------
  // Allocation of arrays to fill the linear system

  typedef KokkosArray::Array< scalar_type[ElementNodeCount][ElementNodeCount] , device_type > elem_matrices_type ;
  typedef KokkosArray::Array< scalar_type[ElementNodeCount]                   , device_type >  elem_vectors_type ;

  elem_matrices_type elem_matrices ; // Jacobian matrices
  elem_vectors_type  elem_vectors ;  // Residual vectors

  if ( element_count ) {
    elem_matrices =
      KokkosArray::create_array< elem_matrices_type >( element_count );

    elem_vectors =
      KokkosArray::create_array< elem_vectors_type >( element_count );
  }

  //------------------------------------
  // For boundary condition set the correct values in the solution vector
  //   The 'zmin' face is assigned to 'T_zmin'.
  //   The 'zmax' face is assigned to 'T_zmax'.
  //   The resulting solution is one dimensional along the 'Z' axis.

  DirichletSolutionFunctor::apply( nodal_solution , mesh , 
                                   exact_solution.zmin , 
                                   exact_solution.zmax ,
                                   exact_solution.T_zmin ,
                                   exact_solution.T_zmax );
    
  for(;;) { // Nonlinear loop

#if defined( HAVE_MPI )

    { //------------------------------------
      // Import off-processor nodal solution values
      // for residual and jacobian computations

      KokkosArray::AsyncExchange< typename vector_type::value_type , Device ,
                                  KokkosArray::ParallelDataMap >
        exchange( mesh.parallel_data_map , 1 );

      KokkosArray::PackArray< vector_type >
        ::pack( exchange.buffer() ,
                mesh.parallel_data_map.count_interior ,
                mesh.parallel_data_map.count_send ,
                nodal_solution );

      exchange.setup();

      exchange.send_receive();

      KokkosArray::UnpackArray< vector_type >
        ::unpack( nodal_solution , exchange.buffer() ,
                  mesh.parallel_data_map.count_owned ,
                  mesh.parallel_data_map.count_receive );
    }

#endif

    //------------------------------------
    // Compute element matrices and vectors:

    wall_clock.reset();

    ElementFunctor::apply( mesh ,
                           elem_matrices , 
                           elem_vectors ,
                           nodal_solution ,
                           exact_solution.K );

    device_type::fence();
    perf_data.elem_time += comm::max( machine , wall_clock.seconds() );

    //------------------------------------
    // Fill linear system coefficients:

    wall_clock.reset();

    fill( 0 , jacobian.coefficients );
    fill( 0 , residual );

    GatherFill< matrix_type , mesh_type >::apply( jacobian , 
                                                  residual ,
                                                  mesh , 
                                                  element_map , 
                                                  elem_matrices , 
                                                  elem_vectors );

    device_type::fence();
    perf_data.matrix_gather_fill_time += comm::max( machine , wall_clock.seconds() );

    // Apply boundary conditions:

    wall_clock.reset();

    // Updates jacobian matrix to 1 on the diagonal, zero elsewhere,
    // and 0 in the residual due to the solution vector having the correct value
    DirichletResidualFunctor::apply( jacobian, residual, mesh ,
                                     exact_solution.zmin ,
                                     exact_solution.zmax );

    device_type::fence();
    perf_data.matrix_boundary_condition_time +=
      comm::max( machine , wall_clock.seconds() );

    //------------------------------------
    // Has the residual converged?

    residual_norm = sqrt( dot(mesh.parallel_data_map, residual) );

    if ( 0 == newton_iteration_count ) {
      residual_norm_init = residual_norm ;
    }

    if ( residual_norm / residual_norm_init < newton_tolerance ) {
      break ;
    }

    //------------------------------------
    // Solve linear sytem

    size_t cg_iteration_count = 0 ;
    double cg_residual_norm = 0 ;

    cgsolve( mesh.parallel_data_map ,
             jacobian , residual , delta ,
             cg_iteration_count ,
             cg_residual_norm ,
             cg_iteration_time ,
             cg_iteration_limit , cg_tolerance ) ;

    perf_data.cg_iteration_time += cg_iteration_time ;
    cg_iteration_count_total += cg_iteration_count ;

    // Update non-linear solution with delta...
    // delta is : - Dx = [Jacobian]^1 * Residual which is the negative update
    // LaTeX:
    // \vec {x}_{n+1} = \vec {x}_{n} - ( - \Delta \vec{x}_{n} )
    // text:
    // x[n+1] = x[n] + Dx

    waxpby( mesh.parallel_data_map,
            1.0, nodal_solution,
           -1.0, delta, nodal_solution);

    ++newton_iteration_count ;

    if ( newton_iteration_limit < newton_iteration_count ) {
      break ;
    }
  };

  if ( newton_iteration_count ) {
    perf_data.elem_time /= newton_iteration_count ;
    perf_data.matrix_gather_fill_time /= newton_iteration_count ;
    perf_data.matrix_boundary_condition_time /= newton_iteration_count ;
  }

  if ( cg_iteration_count_total ) {
    perf_data.cg_iteration_time /= cg_iteration_count_total ;
  }

  perf_data.newton_iteration_count = newton_iteration_count ;
  perf_data.cg_iteration_count = cg_iteration_count_total ;

  //------------------------------------

  {
    // For extracting the nodal solution and its coordinates:

    typename mesh_type::node_coords_type::HostMirror node_coords_host =
      KokkosArray::create_mirror( mesh.node_coords );

    typename vector_type::HostMirror nodal_solution_host =
      KokkosArray::create_mirror( nodal_solution );

    KokkosArray::deep_copy( node_coords_host , mesh.node_coords );
    KokkosArray::deep_copy( nodal_solution_host , nodal_solution );

    double tmp = 0 ;

    for ( size_t i = 0 ; i < mesh.parallel_data_map.count_owned ; ++i ) {
      const coordinate_scalar_type x = node_coords_host(i,0);
      const coordinate_scalar_type y = node_coords_host(i,1);
      const coordinate_scalar_type z = node_coords_host(i,2);

      const double Tx = exact_solution(z);
      const double Ts = nodal_solution_host(i);
      const double Te = std::abs( Tx - Ts ) / std::abs( Tx );

      tmp = std::max( tmp , Te );

      if ( print_error && 0.02 < Te ) {
        std::cout << "  node( " << x << " " << y << " " << z << " ) = "
                  << Ts << " != exact_solution " << Tx
                  << std::endl ;
      }
    }
    perf_data.error_max = comm::max( machine , tmp );
  }

  return perf_data ;
}

//----------------------------------------------------------------------------

template< typename Scalar , class Device >
void driver( const char * label ,
             comm::Machine machine , int beg , int end , int runs )
{
  if ( beg == 0 || end == 0 || runs == 0 ) return ;

  if ( comm::rank( machine ) == 0 ) {
    std::cout << std::endl ;
    std::cout << "\"KokkosArray::HybridFE::NonLinear " << label << "\"" << std::endl;
    std::cout
      << "\"Size\" ,  \"Meshing\" ,  \"Graphing\" , \"Element\" , \"Fill\" ,   \"Boundary\" ,  \"CG-Iter\" , \"CG-Iter\" , \"Newton-Iter\" , \"Max-node-error\"" 
      << std::endl
      << "\"nodes\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\" , \"millisec\" , \"total-count\" , \"total-count\" , \"ratio\""
      << std::endl ;
  }

  for(int i = beg ; i < end ; i *= 2 )
  {
    const int ix = 1 + std::max( 1 , (int) cbrt( ((double) i) / 2.0 ) );
    const int iy = 1 + ix ;
    const int iz = 2 * iy ;
    const int n  = ix * iy * iz ;

    PerformanceData perf_data , perf_best ;

    for(int j = 0; j < runs; j++){

      perf_data = run<Scalar,Device>(machine,ix,iy,iz, false );

      if( j == 0 ) {
        perf_best = perf_data ;
      }
      else {
        perf_best.best( perf_data );
      }
    }

    if ( comm::rank( machine ) == 0 ) {

      std::cout << std::setw(8) << n << " , "
                << std::setw(10) << perf_best.mesh_time * 1000 << " , "
                << std::setw(10) << perf_best.graph_time * 1000 << " , "
                << std::setw(10) << perf_best.elem_time * 1000 << " , "
                << std::setw(10) << perf_best.matrix_gather_fill_time * 1000 << " , "
                << std::setw(10) << perf_best.matrix_boundary_condition_time * 1000 << " , "
                << std::setw(10) << perf_best.cg_iteration_time * 1000 << " , "
                << std::setw(7) << perf_best.cg_iteration_count << " , "
                << std::setw(3) << perf_best.newton_iteration_count << " , "
                << std::setw(10) << perf_best.error_max
                << std::endl ;
    }
  }
}

//----------------------------------------------------------------------------

} /* namespace NonLinear */
} /* namespace HybridFEM */


#endif /* #ifndef HYBRIDFEM_IMPLICIT_HPP */

