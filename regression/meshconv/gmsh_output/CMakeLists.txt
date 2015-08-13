# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(exo2gmsh ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.exo gmsh-check_box_24.txt.std
                               gmsh.ndiff.cfg
                    ARGS -i box_24.exo -o box_24.msh -v
                    POSTPROCESS_PROG ${GMSH}
                    POSTPROCESS_PROG_ARGS -check box_24.msh
                    POSTPROCESS_PROG_OUTPUT gmsh-check.txt
                    TEXT_BASELINE gmsh-check_box_24.txt.std
                    TEXT_RESULT gmsh-check.txt
                    TEXT_DIFF_PROG_CONF gmsh.ndiff.cfg)

add_regression_test(gmshtxt2gmsh ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.txt.msh gmsh-check_box_24.txt.std
                               gmsh.ndiff.cfg
                    ARGS -i box_24.txt.msh -o box_24.msh -v
                    POSTPROCESS_PROG ${GMSH}
                    POSTPROCESS_PROG_ARGS -check box_24.msh
                    POSTPROCESS_PROG_OUTPUT gmsh-check.txt
                    TEXT_BASELINE gmsh-check_box_24.txt.std
                    TEXT_RESULT gmsh-check.txt
                    TEXT_DIFF_PROG_CONF gmsh.ndiff.cfg)

add_regression_test(netgen2gmsh ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES torus.mesh gmsh-check_torus.txt.std
                               gmsh.ndiff.cfg
                    ARGS -i torus.mesh -o torus.msh -v
                    POSTPROCESS_PROG ${GMSH}
                    POSTPROCESS_PROG_ARGS -check torus.msh
                    POSTPROCESS_PROG_OUTPUT gmsh-check.txt
                    TEXT_BASELINE gmsh-check_torus.txt.std
                    TEXT_RESULT gmsh-check.txt
                    TEXT_DIFF_PROG_CONF gmsh.ndiff.cfg)
