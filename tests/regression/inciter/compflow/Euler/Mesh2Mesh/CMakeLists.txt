
  ## Test fails with new solution transfer logic in Discretization. Overset is
  ## not intended for use with alecg anyway.
  #add_regression_test(compflow_euler_m2m ${INCITER_EXECUTABLE}
  #                    NUMPES 3
  #                    INPUTFILES m2m.q unitcube_1k.exo sphere_r0375_2345.exo
  #                    ARGS -c m2m.q -v
  #                    BIN_BASELINE cube.0.e-s.0.3.0
  #                                 cube.0.e-s.0.3.1
  #                                 cube.0.e-s.0.3.2
  #                                 sphere.1.e-s.0.3.0
  #                                 sphere.1.e-s.0.3.1
  #                                 sphere.1.e-s.0.3.2
  #                    BIN_RESULT out.0.e-s.0.3.0
  #                               out.0.e-s.0.3.1
  #                               out.0.e-s.0.3.2
  #                               out.1.e-s.0.3.0
  #                               out.1.e-s.0.3.1
  #                               out.1.e-s.0.3.2
  #                    BIN_DIFF_PROG_ARGS -m
  #                    BIN_DIFF_PROG_CONF m2m_cube.cfg
  #                                       m2m_cube.cfg
  #                                       m2m_cube.cfg
  #                                       m2m_sphere.cfg
  #                                       m2m_sphere.cfg
  #                                       m2m_sphere.cfg
  #                    LABELS "alecg;m2m")

  add_regression_test(compflow_euler_icoverset ${INCITER_EXECUTABLE}
                      NUMPES 3
                      INPUTFILES ic_overset.q freestream_BGmesh_15k.exo sphere_OSmesh_12k.exo
                      ARGS -c ic_overset.q -v
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

  add_regression_test(compflow_euler_icoverset_intergrid ${INCITER_EXECUTABLE}
                      NUMPES 3
                      INPUTFILES ic_overset_intergrid.q freestream_BGmesh_15k.exo sphere_OSmesh_12k.exo
                      ARGS -c ic_overset_intergrid.q -v
                      BIN_BASELINE background_intergrid_pe3.std.exo.0
                                   background_intergrid_pe3.std.exo.1
                                   background_intergrid_pe3.std.exo.2
                                   overset_intergrid_pe3.std.exo.0
                                   overset_intergrid_pe3.std.exo.1
                                   overset_intergrid_pe3.std.exo.2
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
                      INPUTFILES movingmesh_overset.q freestream_BGmesh_15k.exo sphere_OSmesh_12k.exo
                      ARGS -c movingmesh_overset.q -v
                      BIN_BASELINE statbackground_pe4.std.exo.0
                                   statbackground_pe4.std.exo.1
                                   statbackground_pe4.std.exo.2
                                   statbackground_pe4.std.exo.3
                                   mvngoverset0_pe4.std.exo.0
                                   mvngoverset0_pe4.std.exo.1
                                   mvngoverset0_pe4.std.exo.2
                                   mvngoverset0_pe4.std.exo.3
                                   mvngoverset39_pe4.std.exo.0
                                   mvngoverset39_pe4.std.exo.1
                                   mvngoverset39_pe4.std.exo.2
                                   mvngoverset39_pe4.std.exo.3
                                   mvngoverset114_pe4.std.exo.0
                                   mvngoverset114_pe4.std.exo.1
                                   mvngoverset114_pe4.std.exo.2
                                   mvngoverset114_pe4.std.exo.3
                      BIN_RESULT out.0.e-s.0.4.0
                                 out.0.e-s.0.4.1
                                 out.0.e-s.0.4.2
                                 out.0.e-s.0.4.3
                                 out.1.e-s.0.4.0
                                 out.1.e-s.0.4.1
                                 out.1.e-s.0.4.2
                                 out.1.e-s.0.4.3
                                 out.1.e-s.39.4.0
                                 out.1.e-s.39.4.1
                                 out.1.e-s.39.4.2
                                 out.1.e-s.39.4.3
                                 out.1.e-s.114.4.0
                                 out.1.e-s.114.4.1
                                 out.1.e-s.114.4.2
                                 out.1.e-s.114.4.3
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
