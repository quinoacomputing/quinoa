inciter = {

  title = "Initial uniform mesh refinement",

  nstep = 1,     -- Max number of time steps
  cfl = 0.8,     -- CFL coefficient
  ttyi = 1,      -- TTY output interval
  scheme = "alecg",

  partitioning = "mj",

  amr = {
    t0ref = true,
    initial = {"coords", "coords"},
    coords = {
      xminus = 0.5
    }
  },

  transport = {
    physics = "advection",
    problem = "slot_cyl"
  },

  field_output = {
    interval = 1
  }

}
