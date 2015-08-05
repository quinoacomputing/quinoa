# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 1
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 2
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormal_u0.9 ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES skew.q stat.txt.std skew.ndiff.cfg
                    ARGS -c skew.q -v -u 0.9
                    TEXT_BASELINE stat.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF skew.ndiff.cfg)

add_regression_test(SkewNormalPDF ${WALKER_EXECUTABLE}
                    NUMPES 8
                    INPUTFILES skew_pdf.q pdf_p1.txt.std pdf_p2.txt.std
                               skew_pdf.ndiff.cfg
                    ARGS -c skew_pdf.q -v
                    LABELS verification
                    TEXT_BASELINE pdf_p1.txt.std pdf_p2.txt.std
                    TEXT_RESULT pdf_p1.txt pdf_p2.txt
                    TEXT_DIFF_PROG_CONF skew_pdf.ndiff.cfg)
