# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

if (HAS_MKL)
  add_regression_test(OrnsteinUhlenbeck_mkl ${WALKER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES ou_mkl.q stat.txt.std ou.ndiff.cfg
                      ARGS -c ou_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

  add_regression_test(OrnsteinUhlenbeck_mkl ${WALKER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES ou_mkl.q stat.txt.std ou.ndiff.cfg
                      ARGS -c ou_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

  add_regression_test(OrnsteinUhlenbeck_mkl ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES ou_mkl.q stat.txt.std ou.ndiff.cfg
                      ARGS -c ou_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

  add_regression_test(OrnsteinUhlenbeck_mkl_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES ou_mkl.q stat.txt.std ou.ndiff.cfg
                      ARGS -c ou_mkl.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

  add_regression_test(OrnsteinUhlenbeckPDF_mkl ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES ou_pdf_mkl.q pdf_f1.txt.std ou_pdf.ndiff.cfg
                      ARGS -c ou_pdf_mkl.q -v
                      LABELS verification stringent
                      TEXT_BASELINE pdf_f1.txt.std
                      TEXT_RESULT pdf_f1.txt
                      TEXT_DIFF_PROG_CONF ou_pdf.ndiff.cfg)
endif()

if (HAS_RGNSSE2)
  add_regression_test(OrnsteinUhlenbeck_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES ou_rngsse.q stat.txt.std ou.ndiff.cfg
                      ARGS -c ou_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

  add_regression_test(OrnsteinUhlenbeck_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES ou_rngsse.q stat.txt.std ou.ndiff.cfg
                      ARGS -c ou_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

  add_regression_test(OrnsteinUhlenbeck_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES ou_rngsse.q stat.txt.std ou.ndiff.cfg
                      ARGS -c ou_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

  add_regression_test(OrnsteinUhlenbeck_rngsse_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES ou_rngsse.q stat.txt.std ou.ndiff.cfg
                      ARGS -c ou_rngsse.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

  add_regression_test(OrnsteinUhlenbeckPDF_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES ou_pdf_rngsse.q pdf_f1.txt.std ou_pdf.ndiff.cfg
                      ARGS -c ou_pdf_rngsse.q -v
                      LABELS verification stringent
                      TEXT_BASELINE pdf_f1.txt.std
                      TEXT_RESULT pdf_f1.txt
                      TEXT_DIFF_PROG_CONF ou_pdf.ndiff.cfg)
endif()

add_regression_test(OrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeck_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES ou.q stat.txt.std ou.ndiff.cfg
                    ARGS -c ou.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF ou.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeckPDF ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES ou_pdf.q pdf_f1.txt.std ou_pdf.ndiff.cfg
                    ARGS -c ou_pdf.q -v
                    LABELS verification stringent
                    TEXT_BASELINE pdf_f1.txt.std
                    TEXT_RESULT pdf_f1.txt
                    TEXT_DIFF_PROG_CONF ou_pdf.ndiff.cfg)

add_regression_test(OrnsteinUhlenbeckPDF_exo ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES ou_pdf_exo.q pdf_f2.exo.std pdf_f3o.exo.std
                               pdf_f3c.exo.std exodiff_f2.cfg exodiff_f3o.cfg
                               exodiff_f3c.cfg
                    ARGS -c ou_pdf_exo.q -v
                    LABELS verification stringent
                    BIN_BASELINE pdf_f2.exo.std pdf_f3o.exo.std pdf_f3c.exo
                    BIN_RESULT pdf_f2.exo pdf_f3o.exo pdf_f3c.exo
                    BIN_DIFF_PROG_CONF exodiff_f2.cfg exodiff_f3o.cfg
                                       exodiff_f3c.cfg)
