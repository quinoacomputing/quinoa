# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(asc2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube_coarse.asc_mesh
                    ARGS -i cube_coarse.asc_mesh -o cube_coarse.exo -v
                    BIN_BASELINE cube_coarse.exo.std
                    BIN_RESULT cube_coarse.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(netgen2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.mesh
                    ARGS -i box_24.mesh -o box_24.exo -v
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(gmshtxt2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.txt.msh
                    ARGS -i box_24.txt.msh -o box_24.exo -v
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(gmshbin2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.msh
                    ARGS -i box_24.msh -o box_24.exo -v
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(multiblockexo2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_5blocks.exo
                    ARGS -i shear_5blocks.exo -o shear.exo -v
                    BIN_BASELINE shear.exo.std
                    BIN_RESULT shear.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

# With reordering

add_regression_test(asc2reordexo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube_coarse.asc_mesh
                    ARGS -i cube_coarse.asc_mesh -o cube_coarse_reord.exo -v -r
                    BIN_BASELINE cube_coarse_reord.exo.std
                    BIN_RESULT cube_coarse_reord.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)
