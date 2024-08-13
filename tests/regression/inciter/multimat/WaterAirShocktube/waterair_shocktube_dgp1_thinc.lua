inciter = {

  title = "Water-air shock-tube",

  nstep = 25,
  cfl = 0.7,
  ttyi = 10,    -- TTY output interval
  scheme = "dgp1",
  limiter = "vertexbasedp1",
  shock_detector_coeff = 1.0,
  limsol_projection = false,
  lowspeed_kp = 1.0,

  partitioning = "mj",

  multimat = {
    physics = "euler",
    problem = "waterair_shocktube",
    prelax = 1,
    prelax_timescale = 0.25,
    intsharp = 1,
    intsharp_param = 2.5,
    nmat = 2
  },

  material = {
    {
      id = { 1, 2 },
      gamma = { 4.4, 1.4 },  -- ratio of specific heats
      cv = { 951.36, 717.5 },  -- specific heat at const volume
      pstiff = { 6.0e8, 0.0 }  -- sg-eos stiffness parameter
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
