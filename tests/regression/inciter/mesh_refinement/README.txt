# vim: filetype=cmake:
#
# Rough naming convention for AMR regression tests
#
# Example: amr_t0ref_uu_trans_reord_dg_u0.8_migr
#
# Legend:
# * amr - exercising adaptive mesh refinement (AMR)
# * t0ref/dtref - AMR at t<0 (t0ref), AMR at t>0 (dtref)
# * uu - 2 uniform refinement, 1:8 for all tets,
#        other examples: u, ui, ii, iu, iii, cu, cli, eci, ..., legend:
#        - u: uniform refinement (1->8) (keyword: uniform)
#        - i: initial conditions based adaptive refinement (keyword: ic)
#        - c: coordinate based refinement (keyword: coords)
#        - e: refine a list of tagged edges (keyword: edgelist)
#        - d : uniform de-refinement (keyword: uniform_derefine)
# * trans - type of physics, e.g., compflow, transport
# * reord - perform PE-locality mesh node reordering during setup
# * dg - discontinuous Galerkin discretization, e.g., dgp0, dpg2, ...
# * u0.8 - Non-zero virtualization u = 0.8
# * migr - Use migration
