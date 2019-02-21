# Naming convention for AMR regression tests
#
# Example: amr_t0ref_uu_trans_reord_dg_u0.4
#
# Legend:
# * amr - exercising adaptive mesh refinement (AMR)
# * t0ref - AMR before t=0 (initial AMR),
#          other possibilities instead of 't0ref':
#          - dtref - AMR durging time stepping
# * uu - uniform refinement (1:8 for all tets), two uniform steps
#        other possibilities instead of 'uu' are arbitrary combinations of
#        - i - initial conditions based adaptive refinement
#        - c - coordinate based refinement
#        - l - refine a list of tagged edges
# * trans - type of physics, e.g., compflow, transport
#           other possibilities instead of 'trans':
#           - compflow - compressible flow
# * reord - perform PE-locality mesh node reordering during setup
# * dg - discontinuous Galerkin discretization scheme, 1st order, DG = DG(P0),
#        other possibilities instead of 'dg':
#        - dgp1 - DG(P1), dgp2 - 2nd-order DG, etc.
#        - diagcg - lumped-mass (diagonal) continuous Galerkin
# * u0.4 - Non-zero virtualization u = 0.4
