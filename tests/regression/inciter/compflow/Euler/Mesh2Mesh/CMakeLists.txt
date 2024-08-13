# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Parallel

  add_regression_test(compflow_euler_icoverset ${INCITER_EXECUTABLE}
                      NUMPES 3
                      INPUTFILES ic_overset.lua freestream_BGmesh_15k.exo sphere_OSmesh_12k.exo
                      ARGS -c ic_overset.lua -v
                      BIN_BASELINE background_pe3.std.exo.0
                                   background_pe3.std.exo.1
                                   background_pe3.std.exo.2
                                   overset_pe3.std.exo.0
                                   overset_pe3.std.exo.1
                                   overset_pe3.std.exo.2
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
                      LABELS "oversetfe;m2m")

  add_regression_test(compflow_euler_movingmesh_overset ${INCITER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES movingmesh_overset.lua freestream_BGmesh_15k.exo sphere_OSmesh_12k.exo
                      ARGS -c movingmesh_overset.lua -v
                      BIN_BASELINE statbackground_pe4.std.exo.0
                                   statbackground_pe4.std.exo.1
                                   statbackground_pe4.std.exo.2
                                   statbackground_pe4.std.exo.3
                                   mvngoverset0_pe4.std.exo.0
                                   mvngoverset0_pe4.std.exo.1
                                   mvngoverset0_pe4.std.exo.2
                                   mvngoverset0_pe4.std.exo.3
                                   mvngoverset13_pe4.std.exo.0
                                   mvngoverset13_pe4.std.exo.1
                                   mvngoverset13_pe4.std.exo.2
                                   mvngoverset13_pe4.std.exo.3
                                   mvngoverset38_pe4.std.exo.0
                                   mvngoverset38_pe4.std.exo.1
                                   mvngoverset38_pe4.std.exo.2
                                   mvngoverset38_pe4.std.exo.3
                      BIN_RESULT out.0.e-s.0.4.0
                                 out.0.e-s.0.4.1
                                 out.0.e-s.0.4.2
                                 out.0.e-s.0.4.3
                                 out.1.e-s.0.4.0
                                 out.1.e-s.0.4.1
                                 out.1.e-s.0.4.2
                                 out.1.e-s.0.4.3
                                 out.1.e-s.13.4.0
                                 out.1.e-s.13.4.1
                                 out.1.e-s.13.4.2
                                 out.1.e-s.13.4.3
                                 out.1.e-s.38.4.0
                                 out.1.e-s.38.4.1
                                 out.1.e-s.38.4.2
                                 out.1.e-s.38.4.3
                      BIN_DIFF_PROG_ARGS -m
                      BIN_DIFF_PROG_CONF m2m_cube.cfg
                                         m2m_cube.cfg
                                         m2m_cube.cfg
                                         m2m_cube.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                                         m2m_sphere.cfg
                      TEXT_BASELINE mvngoverset_hist.p1.std
                      TEXT_RESULT out.1.hist.p1
                      TEXT_DIFF_PROG_CONF mvngoverset_hist.ndiff.cfg
                      LABELS "oversetfe;m2m")
