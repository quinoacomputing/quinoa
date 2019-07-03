# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

add_regression_test(MixDirichlet_constpdf ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_constpdf.q
                    ARGS -c mixdir_constpdf.q -v
                    POSTPROCESS_PROG sort
                    POSTPROCESS_PROG_ARGS -g pdf_p1.txt COMMAND grep -v \#
                    POSTPROCESS_PROG_OUTPUT pdf_p1_sorted.txt
                    TEXT_BASELINE p1_analytic.txt
                    TEXT_RESULT pdf_p1_sorted.txt
                    TEXT_DIFF_PROG_CONF mixdir.ndiff.cfg
                    LABELS verification)

add_regression_test(MixDirichlet_constpdf_ijsa ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_constpdf_ijsa.q
                    ARGS -c mixdir_constpdf_ijsa.q -v
                    POSTPROCESS_PROG sort
                    POSTPROCESS_PROG_ARGS -g pdf_p1.txt COMMAND grep -v \#
                    POSTPROCESS_PROG_OUTPUT pdf_p1_sorted.txt
                    TEXT_BASELINE p1_analytic_ijsa.txt
                    TEXT_RESULT pdf_p1_sorted.txt
                    TEXT_DIFF_PROG_CONF mixdir.ndiff.cfg
                    LABELS verification)
