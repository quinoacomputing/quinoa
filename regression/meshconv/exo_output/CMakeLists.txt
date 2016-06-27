# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(netgen2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.mesh box_24.exo.std exodiff.cfg
                    ARGS -i box_24.mesh -o box_24.exo -v
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(gmshtxt2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.txt.msh box_24.exo.std exodiff.cfg
                    ARGS -i box_24.txt.msh -o box_24.exo -v
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(gmshbin2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES box_24.msh box_24.exo.std exodiff.cfg
                    ARGS -i box_24.msh -o box_24.exo -v
                    BIN_BASELINE box_24.exo.std
                    BIN_RESULT box_24.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)

add_regression_test(multiblockexo2exo ${MESHCONV_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES shear_5blocks.exo shear.exo.std exodiff.cfg
                    ARGS -i shear_5blocks.exo -o shear.exo -v
                    BIN_BASELINE shear.exo.std
                    BIN_RESULT shear.exo
                    BIN_DIFF_PROG_CONF exodiff.cfg)
