inciter = {

  title = "Initial uniform mesh refinement",

  nstep = 1,     -- Max number of time steps
  cfl = 0.2,     -- CFL coefficient
  ttyi = 1,      -- TTY output interval
  scheme = "alecg",

  pelocal_reorder = true,
  partitioning = "mj",

  amr = {
    t0ref = true,
    initial = {"initial_conditions", "initial_conditions"},
    error = "jump",
    tol_refine = 0.8
  },

  transport = {
    physics = "advection",
    problem = "slot_cyl"
  },

  field_output = {
    interval = 1
  }

}
