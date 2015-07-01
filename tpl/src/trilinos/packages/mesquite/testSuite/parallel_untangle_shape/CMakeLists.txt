# the "2" is the number of MPI procs
MSQ_ADD_TEST_PARALLEL( par_hex_untangle_shape par_hex_untangle_shape_main 2 par_hex_untangle_shape.cpp)
MSQ_ADD_TEST_PARALLEL( par_hex_untangle_shape_1 par_hex_untangle_shape_main_1 1 par_hex_untangle_shape.cpp)
