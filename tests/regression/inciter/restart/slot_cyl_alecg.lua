inciter = {

  title = "Zalesak's slotted cylinder",

  nstep = 5,     -- Max number of time steps
  cfl = 0.8,     -- CFL coefficient
  ttyi = 1 ,     -- TTY output interval

  scheme = "alecg",

  transport = {
    physics = "advection",
    problem = "slot_cyl"
  },

  field_output = {
    interval = 1,
    filetype = "exodusii",
    nodevar = {"analytic", "C1"},
    nodealias = {"", "c0_numerical" }
  }

}
