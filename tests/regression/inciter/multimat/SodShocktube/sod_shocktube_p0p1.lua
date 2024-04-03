inciter = {

  title = "Sod shock-tube",

  nstep = 25,   -- Max number of time steps
  dt = 1.0e-3,  -- Time step size
  ttyi = 10,    -- TTY output interval
  scheme = "p0p1",
  limiter = "vertexbasedp1",
  shock_detector_coeff = 0.0,
  limsol_projection = false,

  partitioning = "mj",

  multimat = {
    physics = "euler",
    problem = "sod_shocktube",
    prelax = 0,
    nmat = 2
  },

  material = {
    {
      id = { 1, 2 },
      gamma = { 1.4, 1.4 }  -- ratio of specific heats
    }
  },

  bc = {
    {
      extrapolate = { 1, 3 },
      symmetry = { 2, 4, 5, 6 }
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 25,
    elemvar = {
      "F1",
      "F2",
      "density",
      "pressure",
      "specific_total_energy",
      "x-velocity",
      "y-velocity",
      "z-velocity"
    },
    elemalias = {
      "volfrac1_numerical",
      "volfrac2_numerical",
      "density_numerical",
      "pressure_numerical",
      "total_energy_density_numerical",
      "x-velocity_numerical",
      "y-velocity_numerical",
      "z-velocity_numerical"
    }
  }

}
