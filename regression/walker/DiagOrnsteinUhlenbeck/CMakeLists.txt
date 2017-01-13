# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

if (HAS_MKL)
  add_regression_test(DiagOrnsteinUhlenbeck_mkl ${WALKER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES diagou_mkl.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeck_mkl ${WALKER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES diagou_mkl.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeck_mkl ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES diagou_mkl.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_mkl.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeck_mkl_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES diagou_mkl.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_mkl.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeck_mkl_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES diagou_mkl.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_mkl.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeckPDF_mkl ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES diagou_pdf_mkl.q pdf_f2.txt.std
                                 diagou_pdf.ndiff.cfg
                      ARGS -c diagou_pdf_mkl.q -v
                      LABELS verification
                      TEXT_BASELINE pdf_f2.txt.std
                      TEXT_RESULT pdf_f2.txt
                      TEXT_DIFF_PROG_CONF diagou_pdf.ndiff.cfg)
endif()

if (HAS_RNGSSE2)
  add_regression_test(DiagOrnsteinUhlenbeck_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 1
                      INPUTFILES diagou_rngsse.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeck_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 2
                      INPUTFILES diagou_rngsse.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeck_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES diagou_rngsse.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_rngsse.q -v
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeck_rngsse_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 4
                      INPUTFILES diagou_rngsse.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_rngsse.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeck_rngsse_u0.9 ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES diagou_rngsse.q stat.txt.std diagou.ndiff.cfg
                      ARGS -c diagou_rngsse.q -v -u 0.9
                      TEXT_BASELINE stat.txt.std
                      TEXT_RESULT stat.txt
                      TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

  add_regression_test(DiagOrnsteinUhlenbeckPDF_rngsse ${WALKER_EXECUTABLE}
                      NUMPES 8
                      INPUTFILES diagou_pdf_rngsse.q pdf_f2.txt.std
                                 diagou_pdf.ndiff.cfg
                      ARGS -c diagou_pdf_rngsse.q -v
                      LABELS verification
                      TEXT_BASELINE pdf_f2.txt.std
                      TEXT_RESULT pdf_f2.txt
                      TEXT_DIFF_PROG_CONF diagou_pdf.ndiff.cfg)
endif()

add_regression_test(DiagOrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeck ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeck_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeck_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES diagou.q stat.txt.std diagou.ndiff.cfg
                    ARGS -c diagou.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF diagou.ndiff.cfg)

add_regression_test(DiagOrnsteinUhlenbeckPDF ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES diagou_pdf.q pdf_f2.txt.std
                               diagou_pdf.ndiff.cfg
                    ARGS -c diagou_pdf.q -v
                    LABELS verification
                    TEXT_BASELINE pdf_f2.txt.std
                    TEXT_RESULT pdf_f2.txt
                    TEXT_DIFF_PROG_CONF diagou_pdf.ndiff.cfg)
