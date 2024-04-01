inciter = {

  title = "Zalesak's slotted cylinder - restarted",

  nstep = 10,     -- Max number of time steps
  cfl = 0.8,     -- CFL coefficient
  ttyi = 1 ,     -- TTY output interval

  scheme = "alecg",

  transport = {
    problem = "slot_cyl"
  },

  physics = "advection",

  field_output = {
    interval = 1,
    filetype = "exodusii",
    nodevar = {"analytic", "C1"},
    nodealias = {"", "c0_numerical" }
  }

}
