# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

if (HAS_MKL)
  add_regression_test(SkewNormal_mkl ${WALKER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES skew_mkl.q
                      ARGS -c skew_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

  add_regression_test(SkewNormal_mkl ${WALKER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES skew_mkl.q
                      ARGS -c skew_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

  add_regression_test(SkewNormal_mkl ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES skew_mkl.q
                      ARGS -c skew_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

  add_regression_test(SkewNormal_mkl_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES skew_mkl.q
                      ARGS -c skew_mkl.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

  add_regression_test(SkewNormalPDF_mkl ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES skew_pdf_mkl.q
                      ARGS -c skew_pdf_mkl.q -v
                      LABELS verification stringent
                      TEXT_BASELINE pdf_p1.txt.std pdf_p2.txt.std
                      TEXT_RESULT pdf_p1.txt pdf_p2.txt
                      TEXT_DIFF_PROG_CONF skew_pdf.ndiff.cfg)
endif()

if (HAS_RNGSSE2)
  add_regression_test(SkewNormal_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES skew_rngsse.q
                      ARGS -c skew_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

  add_regression_test(SkewNormal_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES skew_rngsse.q
                      ARGS -c skew_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

  add_regression_test(SkewNormal_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES skew_rngsse.q
                      ARGS -c skew_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

  add_regression_test(SkewNormal_rngsse_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES skew_rngsse.q
                      ARGS -c skew_rngsse.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

  add_regression_test(SkewNormalPDF_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES skew_pdf_rngsse.q
                      ARGS -c skew_pdf_rngsse.q -v
                      LABELS verification stringent
                      TEXT_BASELINE pdf_p1.txt.std pdf_p2.txt.std
                      TEXT_RESULT pdf_p1.txt pdf_p2.txt
                      TEXT_DIFF_PROG_CONF skew_pdf.ndiff.cfg)
endif()

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES skew.q
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES skew.q
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES skew.q
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES skew.q
                    ARGS -c skew.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormalPDF ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES skew_pdf.q
                    ARGS -c skew_pdf.q -v
                    LABELS verification stringent
                    TEXT_BASELINE pdf_p1.txt.std pdf_p2.txt.std
                    TEXT_RESULT pdf_p1.txt pdf_p2.txt
                    TEXT_DIFF_PROG_CONF skew_pdf.ndiff.cfg)
