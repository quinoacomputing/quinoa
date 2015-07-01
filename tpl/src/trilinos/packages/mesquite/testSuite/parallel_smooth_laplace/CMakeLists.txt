# the "2" is the number of MPI procs
MSQ_ADD_TEST_PARALLEL( par_hex_smooth_laplace par_hex_smooth_laplace_main 2 par_hex_smooth_laplace.cpp )
MSQ_ADD_TEST_PARALLEL( par_hex_smooth_laplace_1 par_hex_smooth_laplace_main_1 1 par_hex_smooth_laplace.cpp )
