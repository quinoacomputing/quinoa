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
 * @author  H. Carter Edwards  <hcedwar@sandia.gov>
 */

#ifndef util_Parallel_hpp
#define util_Parallel_hpp

// phdmesh_config.h resides in the build directory and contains the
// complete set of #define macros for build-dependent features.

#include <phdmesh_config.h>

//----------------------------------------------------------------------
// Parallel machine

#if defined( HAVE_MPI )

#include <mpi.h>

namespace phdmesh {
typedef MPI_Comm     ParallelMachine ;
typedef MPI_Datatype ParallelDatatype ;

inline ParallelMachine parallel_machine_null() { return MPI_COMM_NULL ; }

inline ParallelMachine parallel_machine_init( int * argc , char *** argv )
{
  MPI_Init( argc , argv );
  return MPI_COMM_WORLD ;
}

inline void parallel_machine_finalize()
{
  MPI_Finalize();
}

}

//----------------------------------------
// Other parallel communication machines go here
// as '#elif defined( PHDMESH_HAS_<name> )'

//----------------------------------------
// Stub for non-parallel

#else

namespace phdmesh {
typedef int ParallelMachine ;
typedef int ParallelDatatype ;

inline ParallelMachine parallel_machine_null() { return 0 ; }

inline ParallelMachine parallel_machine_init( int * , char *** )
{ return 0 ; }

inline void parallel_machine_finalize()
{}

}

#endif

//----------------------------------------------------------------------
// Common parallel machine needs.

namespace phdmesh {

double wall_time();
double wall_dtime( double & );

unsigned parallel_machine_size( ParallelMachine m );

unsigned parallel_machine_rank( ParallelMachine m );

void parallel_machine_barrier( ParallelMachine );

/** Parallel identification, sort by identifier and then processor. */

struct IdentProc {
  unsigned ident ;
  unsigned proc ;

  ~IdentProc() {}

  IdentProc() {}

  IdentProc( const IdentProc & rhs ) : ident(rhs.ident), proc(rhs.proc) {}

  IdentProc( unsigned i , unsigned p ) : ident(i), proc(p) {}

  IdentProc & operator = ( const IdentProc & rhs )
    { ident = rhs.ident ; proc = rhs.proc ; return *this ; }

  bool operator == ( const IdentProc & rhs ) const
    { return ident == rhs.ident && proc == rhs.proc ; }

  bool operator != ( const IdentProc & rhs ) const
    { return ident != rhs.ident || proc != rhs.proc ; }

  bool operator < ( const IdentProc & rhs ) const
    { return ident != rhs.ident ? ident < rhs.ident : proc < rhs.proc ; }

  bool operator > ( const IdentProc & rhs ) const
    { return rhs.operator<( *this ); }

  bool operator <= ( const IdentProc & rhs ) const
    { return ! rhs.operator<( *this ); }

  bool operator >= ( const IdentProc & rhs ) const
    { return ! this->operator<( rhs ); }
};

}

//----------------------------------------------------------------------

#endif

