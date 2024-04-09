inciter = {

  title = "Test GCL with ALECG",

  nstep = 10,    -- Max number of time steps
  ttyi = 1,      -- TTY output interval
  cfl = 0.5,
  scheme = "alecg",

  partitioning = "mj",

  ale = {
    dvcfl = 1.0,
    mesh_velocity = "sine"
  },

  compflow = {
    physics = "euler",
    problem = "user_defined"
  },

  mesh = {
    { filename = "rectangle_01_1.5k.exo" }
  },

  ic = {
    density = 1.0,
    velocity = { 0.0, 0.0, 0.0 },
    pressure = 1.0
  },

  material = {
    {
      gamma = { 1.4 }
    }
  },

  bc = {
    {
      symmetry = { 2, 4, 5, 6 }
    }
  },

  field_output = {
    nodevar = {
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
    },
    interval = 10
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  }

}
