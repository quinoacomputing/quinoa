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
                      INPUTFILES bouncing_sphere_overset.lua halfcube_BGmesh_266k.exo sphere_inhalfcube_OSmesh_11k.exo
                      ARGS -c bouncing_sphere_overset.lua -v
                      BIN_BASELINE statbackground_pe4.std.exo.0
                                   statbackground_pe4.std.exo.1
                                   statbackground_pe4.std.exo.2
                                   statbackground_pe4.std.exo.3
                                   mvngoverset0_pe4.std.exo.0
                                   mvngoverset0_pe4.std.exo.1
                                   mvngoverset0_pe4.std.exo.2
                                   mvngoverset0_pe4.std.exo.3
                                   mvngoverset2_pe4.std.exo.0
                                   mvngoverset2_pe4.std.exo.1
                                   mvngoverset2_pe4.std.exo.2
                                   mvngoverset2_pe4.std.exo.3
                      BIN_RESULT out.0.e-s.0.4.0
                                 out.0.e-s.0.4.1
                                 out.0.e-s.0.4.2
                                 out.0.e-s.0.4.3
                                 out.1.e-s.0.4.0
                                 out.1.e-s.0.4.1
                                 out.1.e-s.0.4.2
                                 out.1.e-s.0.4.3
                                 out.1.e-s.2.4.0
                                 out.1.e-s.2.4.1
                                 out.1.e-s.2.4.2
                                 out.1.e-s.2.4.3
                      BIN_DIFF_PROG_ARGS -m
                      BIN_DIFF_PROG_CONF m2m_moving.cfg
                      LABELS "oversetfe;m2m")
