# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Velocity coefficients: constant imposed mean shear

add_regression_test(Velocity_SLM_HomogeneousShear ${WALKER_EXECUTABLE}
                    NUMPES 1
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

# Velocity coefficients: force stationary PDF

add_regression_test(Velocity_stationary ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES stationary.q
                    ARGS -c stationary.q -v
                    TEXT_BASELINE stationary.stat.txt.std
                                  stationary_pdf_u1.txt.std
                                  stationary_pdf_u1.txt.std
                                  stationary_pdf_u1.txt.std
                    TEXT_RESULT stat.txt
                                pdf_u1.txt
                                pdf_u2.txt
                                pdf_u3.txt
                    TEXT_DIFF_PROG_CONF stationary_velocity.ndiff.cfg
                                        stationary_velocity_pdf.ndiff.cfg
                                        stationary_velocity_pdf.ndiff.cfg
                                        stationary_velocity_pdf.ndiff.cfg)
