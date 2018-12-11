# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(Velocity_SLM_HomogeneousShear ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES slm_homogeneous_shear.q
                    ARGS -c slm_homogeneous_shear.q -v
                    TEXT_BASELINE slm_homogeneous_shear.stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF velocity.ndiff.cfg
                    LABELS verification)

add_regression_test(Velocity_SLM_HomogeneousShear ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES slm_homogeneous_shear.q
                    ARGS -c slm_homogeneous_shear.q -v
                    TEXT_BASELINE slm_homogeneous_shear.stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF velocity.ndiff.cfg
                    LABELS verification)

add_regression_test(Velocity_SLM_HomogeneousShear ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slm_homogeneous_shear.q
                    ARGS -c slm_homogeneous_shear.q -v
                    TEXT_BASELINE slm_homogeneous_shear.stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF velocity.ndiff.cfg
                    LABELS verification)

add_regression_test(Velocity_SLM_HomogeneousShear_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES slm_homogeneous_shear.q
                    ARGS -c slm_homogeneous_shear.q -v -u 0.9
                    TEXT_BASELINE slm_homogeneous_shear.stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF velocity.ndiff.cfg
                    LABELS verification)

add_regression_test(Velocity_GLM_HomogeneousShear ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES glm_homogeneous_shear.q
                    ARGS -c glm_homogeneous_shear.q -v
                    TEXT_BASELINE glm_homogeneous_shear.stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF velocity.ndiff.cfg
                    LABELS verification)

add_regression_test(Velocity_GLM_HomogeneousShear ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES glm_homogeneous_shear.q
                    ARGS -c glm_homogeneous_shear.q -v
                    TEXT_BASELINE glm_homogeneous_shear.stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF velocity.ndiff.cfg
                    LABELS verification)

add_regression_test(Velocity_GLM_HomogeneousShear_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES glm_homogeneous_shear.q
                    ARGS -c glm_homogeneous_shear.q -v -u 0.9
                    TEXT_BASELINE glm_homogeneous_shear.stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF velocity.ndiff.cfg
                    LABELS verification)
