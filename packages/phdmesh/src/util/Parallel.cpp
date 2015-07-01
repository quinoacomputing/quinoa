/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#include <util/Parallel.hpp>

/*--------------------------------------------------------------------*/
/* Parallel operations */

#if defined( HAVE_MPI )

namespace phdmesh {

unsigned parallel_machine_size( ParallelMachine m )
{
  int value = 0 ;
  if ( MPI_SUCCESS != MPI_Comm_size( m , &value ) ) { value = 0 ; }
  return value ;
}

unsigned parallel_machine_rank( ParallelMachine m )
{
  int value = 0 ;
  if ( MPI_SUCCESS != MPI_Comm_rank( m , &value ) ) { value = 0 ; }
  return value ;
}

void parallel_machine_barrier( ParallelMachine m )
{
  MPI_Barrier( m );
}

}

#else

namespace phdmesh {

unsigned parallel_machine_size( ParallelMachine ) { return 1 ; }

unsigned parallel_machine_rank( ParallelMachine ) { return 0 ; }

void parallel_machine_barrier( ParallelMachine ) {}

}

#endif

/*--------------------------------------------------------------------*/
/* Wall time */

namespace phdmesh {

double wall_dtime( double & t )
{
  const double tnew = wall_time();
  const double dt = tnew - t ; t = tnew ;
  return dt ;
}

}

#if ! defined(PHDMESH_NO_SYS_TIME)

#include <stddef.h>
#include <sys/time.h>

namespace phdmesh {

double wall_time()
{
  static bool    first = true ;
  static timeval tp_init ;

  timeval tp ;

  gettimeofday( &tp , reinterpret_cast<struct timezone *>( NULL ) );

  if ( first ) {
    tp_init.tv_usec = tp.tv_usec ;
    tp_init.tv_sec  = tp.tv_sec ;
    first = false ;
  }

  return ( (double)( tp.tv_sec  - tp_init.tv_sec  ) ) +
         ( (double)( tp.tv_usec - tp_init.tv_usec ) / 1000000.0 );
} 

}

#elif defined( HAVE_MPI )

namespace phdmesh {

double wall_time() { return MPI_Wtime(); }

}

#endif

/*--------------------------------------------------------------------*/


