/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_parallel_Parallel_hpp
#define stk_util_parallel_Parallel_hpp

// stk_config.h resides in the build directory and contains the
// complete set of #define macros for build-dependent features.

#include <stk_util/stk_config.h>

//----------------------------------------------------------------------
// Parallel machine

#if defined( STK_HAS_MPI )

#include <mpi.h>

namespace stk_classic {

/** \addtogroup parallel_module
 * @{
 */

/// \todo REFACTOR: Figure out a better way to typedef for non-MPI builds

typedef MPI_Comm     ParallelMachine ;

/// \todo REFACTOR: Figure out a better way to typedef for non-MPI builds

typedef MPI_Datatype ParallelDatatype ;

/**
 * @brief <b>parallel_machine_null</b> returns MPI_COMM_NULL if MPI is enabled.
 *
 * @return			a <b>ParallelMachine</b> ...
 */
inline ParallelMachine parallel_machine_null() { return MPI_COMM_NULL ; }

/**
 * @brief <b>parallel_machine_init</b> calls MPI_Init.
 *
 * @param argc
 *
 * @param argv
 *
 * @return <b>ParallelMachine</b> (MPI_COMM_WORLD)
 */
inline ParallelMachine parallel_machine_init( int * argc , char *** argv )
{
  MPI_Init( argc , argv );
  return MPI_COMM_WORLD ;
}

/**
 * @brief <b>parallel_machine_finalize</b> calls MPI_Finalize.
 *
 */
inline void parallel_machine_finalize()
{
  MPI_Finalize();
}

/** \} */

}

//----------------------------------------
// Other parallel communication machines go here
// as '#elif defined( STK_HAS_<name> )'

//----------------------------------------
// Stub for non-parallel

#else

// Some needed stubs
#define MPI_Comm int
#define MPI_COMM_WORLD 0
#define MPI_COMM_SELF 0
#define MPI_Barrier( a ) (void)a

namespace stk_classic {

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

namespace stk_classic {

/**
 * @brief Member function <b>parallel_machine_size</b> ...
 *
 * @param m			a <b>ParallelMachine</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
unsigned parallel_machine_size( ParallelMachine parallel_machine );

/**
 * @brief Member function <b>parallel_machine_rank</b> ...
 *
 * @param m			a <b>ParallelMachine</b> ...
 *
 * @return			an <b>unsigned int</b> ...
 */
unsigned parallel_machine_rank( ParallelMachine parallel_machine );

/**
 * @brief Member function <b>parallel_machine_barrier</b> ...
 *
 */
void parallel_machine_barrier( ParallelMachine parallel_machine);
}

//----------------------------------------------------------------------

#endif

