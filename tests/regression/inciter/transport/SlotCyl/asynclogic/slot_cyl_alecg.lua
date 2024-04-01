inciter = {

  title = "Zalesak's slotted cylinder",

  nstep = 10,    -- Max number of time steps
  dt = 0.001,    -- Time step size
  ttyi = 1,      -- TTY output interval
  scheme = "alecg",

  transport = {
    problem = "slot_cyl"
  },

  physics = "advection",

  field_output = {
    interval = 5
  }

}
