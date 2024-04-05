inciter = {

  title = "Sedov blast wave",

  nstep = 10,   -- Max number of time steps
  cfl = 0.3,
  ttyi = 5,      -- TTY output interval
  scheme = "pdg",
  limiter = "superbeep1",

  compflow = {
    physics = "euler",
    problem = "sedov_blastwave"
  },

  material = {
    {
      gamma = { 1.4 }
    }
  },

  bc = {
    {
      symmetry = { 1, 2 },
      extrapolate = { 3 }
    }
  },

  pref = {
    ndofmax = 4
  },

  diagnostics = {
    interval = 5,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 5,
    elemvar = {
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
   },
   elemalias = {
      "density_numerical",
      "x-velocity_numerical",
      "y-velocity_numerical",
      "z-velocity_numerical",
      "specific_total_energy_numerical",
      "pressure_numerical"
   }
  }

}
