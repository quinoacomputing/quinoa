
  add_regression_test(compflow_euler_m2m ${INCITER_EXECUTABLE}
                      NUMPES 3
                      INPUTFILES m2m.q unitcube_1k.exo sphere_r0375_2345.exo
                      ARGS -c m2m.q -v
                      BIN_BASELINE cube.0.e-s.0.3.0
                                   cube.0.e-s.0.3.1
                                   cube.0.e-s.0.3.2
                                   sphere.1.e-s.0.3.0
                                   sphere.1.e-s.0.3.1
                                   sphere.1.e-s.0.3.2
                      BIN_RESULT out.0.e-s.0.3.0
                                 out.0.e-s.0.3.1
                                 out.0.e-s.0.3.2
                                 out.1.e-s.0.3.0
                                 out.1.e-s.0.3.1
                                 out.1.e-s.0.3.2
                      BIN_DIFF_PROG_ARGS -m
                      BIN_DIFF_PROG_CONF m2m_cube.cfg
                                         m2m_cube.cfg
                                         m2m_cube.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                      LABELS "alecg;m2m")
