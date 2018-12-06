# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(asc2netgen ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES cube_coarse.asc_mesh
                    ARGS -i cube_coarse.asc_mesh -o cube_coarse.mesh -v
                    TEXT_BASELINE cube_coarse.mesh.std
                    TEXT_RESULT cube_coarse.mesh
                    TEXT_DIFF_PROG_CONF netgen.ndiff.cfg)

add_regression_test(exo2netgen ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.exo
                    ARGS -i box_24.exo -o box_24.mesh -v
                    TEXT_BASELINE box_24.mesh.std
                    TEXT_RESULT box_24.mesh
                    TEXT_DIFF_PROG_CONF netgen.ndiff.cfg)

add_regression_test(gmshtxt2netgen ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.txt.msh
                    ARGS -i box_24.txt.msh -o box_24.mesh -v
                    TEXT_BASELINE box_24.mesh.std
                    TEXT_RESULT box_24.mesh
                    TEXT_DIFF_PROG_CONF netgen.ndiff.cfg)

add_regression_test(gmshbin2netgen ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.msh
                    ARGS -i box_24.msh -o box_24.mesh -v
                    TEXT_BASELINE box_24.mesh.std
                    TEXT_RESULT box_24.mesh
                    TEXT_DIFF_PROG_CONF netgen.ndiff.cfg)
