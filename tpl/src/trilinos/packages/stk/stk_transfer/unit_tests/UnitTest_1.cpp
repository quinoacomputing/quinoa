/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

namespace stk_transfer_unit_tests {

//The following ::testUnit() function is where the actual unit-test is:
//(This one doesn't do much, it's mainly a template to demonstrate how
// to create a unit-test.)
//
//To create another unit-test, copy this file, change all occurrences of
//'UnitTest_1' to something else, and write your test code in its body.

void UnitTest_1( MPI_Comm comm )
{
  int mpi_rank = 0;
  int mpi_size = 1;
  
#ifdef STK_HAS_MPI
  STKUNIT_ASSERT_EQUAL(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank), MPI_SUCCESS);
  STKUNIT_ASSERT_EQUAL(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size), MPI_SUCCESS);
#endif
  
  STKUNIT_ASSERT(mpi_rank < mpi_size);
}

} // namespace stk_transfer_unit_tests

STKUNIT_UNIT_TEST(UnitTestingOfTransfer, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_transfer_unit_tests::UnitTest_1 ( MPI_COMM_WORLD );
}

