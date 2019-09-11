# See cmake/add_regression_test.cmake for documentation on the arguments to
# add_regression_test().

# Prescribe constant SDE coefficients, The PDF must not change in time. Only
# testing a single marginal. The PDF should not change because of constant
# coefficients for S_alpha, b_alpha, and kappa_alpha.

add_regression_test(MixDirichlet_constpdf_nomix ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_constpdf_nomix.q
                    ARGS -c mixdir_constpdf_nomix.q -v
                    POSTPROCESS_PROG sort
                    POSTPROCESS_PROG_ARGS -g pdf_p1.txt COMMAND grep -v \#
                    POSTPROCESS_PROG_OUTPUT pdf_p1_sorted.txt
                    TEXT_BASELINE p1_analytic_nomix.txt
                    TEXT_RESULT pdf_p1_sorted.txt
                    TEXT_DIFF_PROG_CONF mixdir.pdf.ndiff.cfg
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
                    TEXT_DIFF_PROG_CONF mixdir.pdf.ndiff.cfg
                    LABELS verification)

# Constrain S_alpha to force d<rho>/dt = 0. The PDF must not change in time
# because b_alpha and kappa_alpha are consistent with the omega_alpha of the
# initial Dirichlet PDF, and (interestingly), the specification for S_alpha to
# force d<rho>/dt = 0, yields the S_alpha that would be required to leave the
# PDF unchanged. Using heavy-fluid normalization.

add_regression_test(MixDirichlet_homogeneous_nomix_heavy ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_homogeneous_nomix_heavy.q
                    ARGS -c mixdir_homogeneous_nomix_heavy.q -v
                    POSTPROCESS_PROG sort
                    POSTPROCESS_PROG_ARGS -g pdf_p1.txt COMMAND grep -v \#
                    POSTPROCESS_PROG_OUTPUT pdf_p1_sorted.txt
                    TEXT_BASELINE p1_analytic_nomix.txt
                    TEXT_RESULT pdf_p1_sorted.txt
                    TEXT_DIFF_PROG_CONF mixdir.pdf.ndiff.cfg
                    LABELS verification)

add_regression_test(MixDirichlet_homogeneous_ijsa_heavy ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_homogeneous_ijsa_heavy.q
                    ARGS -c mixdir_homogeneous_ijsa_heavy.q -v
                    POSTPROCESS_PROG sort
                    POSTPROCESS_PROG_ARGS -g pdf_p1.txt COMMAND grep -v \#
                    POSTPROCESS_PROG_OUTPUT pdf_p1_sorted.txt
                    TEXT_BASELINE p1_analytic_ijsa.txt
                    TEXT_RESULT pdf_p1_sorted.txt
                    TEXT_DIFF_PROG_CONF mixdir.pdf.ndiff.cfg
                    LABELS verification)

# Constrain S_alpha to force d<rho>/dt = 0. The PDF must not change in time
# because b_alpha and kappa_alpha are consistent with the omega_alpha of the
# initial Dirichlet PDF, and interestingly, the specification for S_alpha to
# force d<rho>/dt = 0, yields the S_alpha that would be required to leave the
# PDF unchanged. Using light-fluid normalization.

# Compared to the same tests using heavy-fluid normalization, light-fluid
# normalization yields a different <rho> but still constant in time. This is,
# however, not tested by the following tests, only one of the marginal PDFs,
# thus comparing to the same baseline as the above tests.

add_regression_test(MixDirichlet_homogeneous_nomix_light ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_homogeneous_nomix_light.q
                    ARGS -c mixdir_homogeneous_nomix_light.q -v
                    POSTPROCESS_PROG sort
                    POSTPROCESS_PROG_ARGS -g pdf_p1.txt COMMAND grep -v \#
                    POSTPROCESS_PROG_OUTPUT pdf_p1_sorted.txt
                    TEXT_BASELINE p1_analytic_nomix.txt
                    TEXT_RESULT pdf_p1_sorted.txt
                    TEXT_DIFF_PROG_CONF mixdir.pdf.ndiff.cfg
                    LABELS verification)

add_regression_test(MixDirichlet_homogeneous_ijsa_light ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_homogeneous_ijsa_light.q
                    ARGS -c mixdir_homogeneous_ijsa_light.q -v
                    POSTPROCESS_PROG sort
                    POSTPROCESS_PROG_ARGS -g pdf_p1.txt COMMAND grep -v \#
                    POSTPROCESS_PROG_OUTPUT pdf_p1_sorted.txt
                    TEXT_BASELINE p1_analytic_ijsa.txt
                    TEXT_RESULT pdf_p1_sorted.txt
                    TEXT_DIFF_PROG_CONF mixdir.pdf.ndiff.cfg
                    LABELS verification)

# Constrain S_alpha to force d<rho>/dt = 0. The PDF now changes in time because
# we start from an initial Dirichlet PDF, given by omega_alpha that does not
# correspond to those specified by b_alpha, and kappa_alpha. Still d<rho/dt = 0
# must still be true as S_alpha is constrained to force this. The means are
# diffed, including <rho>. Using heavy-, and light-fluid normalization. Note
# that the value of <rho> is (correctly) different between heavy-, vs
# light-fluid normalization, but still must be constant in time.

add_regression_test(MixDirichlet_homogeneous_evolving_pdf_heavy
                    ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_homogeneous_evolving_pdf_heavy.q
                    ARGS -c mixdir_homogeneous_evolving_pdf_heavy.q -v
                    TEXT_BASELINE stat_evolving_pdf_heavy.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF mixdir.stat.ndiff.cfg)

add_regression_test(MixDirichlet_homogeneous_evolving_pdf_light
                    ${WALKER_EXECUTABLE}
                    NUMPES 4
                    INPUTFILES mixdir_homogeneous_evolving_pdf_light.q
                    ARGS -c mixdir_homogeneous_evolving_pdf_light.q -v
                    TEXT_BASELINE stat_evolving_pdf_light.txt.std
                    TEXT_RESULT stat.txt
                    TEXT_DIFF_PROG_CONF mixdir.stat.ndiff.cfg)
