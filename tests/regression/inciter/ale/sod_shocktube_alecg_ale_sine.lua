inciter = {

  title = "Sod shock-tube, ALECG, sine mesh motion",

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
    problem = "sod_shocktube"
  },

  mesh = {
    { filename = "rectangle_01_1.5k.exo" }
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
    interval = 10,
    nodevar = {
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  }

}
