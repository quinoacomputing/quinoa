inciter = {

  title = "Zalesak's slotted cylinder",

  nstep = 5,     -- Max number of time steps
  dt = 5.0e-4,   -- Time step size
  ttyi = 1,      -- TTY output interval
  scheme = "dgp1",

  transport = {
    physics = "advection",
    problem = "slot_cyl"
  },

  bc = {
    {
      dirichlet = {1}
    }
  },

  field_output = {
    interval = 1
  }

}
