inciter = {

  title = "Initial uniform mesh refinement",

  nstep = 10,     -- Max number of time steps
  cfl = 0.8,     -- CFL coefficient
  ttyi = 1,      -- TTY output interval
  scheme = "alecg",

  pelocal_reorder = true,
  partitioning = "mj",

  amr = {
    t0ref = true,
    initial = {"uniform"}
  },

  transport = {
    physics = "advection",
    problem = "slot_cyl"
  },

  field_output = {
    interval = 2
  }

}
