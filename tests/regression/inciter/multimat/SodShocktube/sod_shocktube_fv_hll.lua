inciter = {

  title = "Sod shock-tube",

  nstep = 25,
  cfl = 0.5,
  ttyi = 5,  -- TTY output interval
  scheme = "fv",
  limiter = "vertexbasedp1",
  flux = "hll",

  partitioning = "mj",

  multimat = {
    physics = "euler",
    problem = "sod_shocktube",
    prelax = 0,
    nmat = 2
  },

  material = {
    {
      eos = "stiffenedgas",
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
      "material_indicator",
      "density",
      "pressure",
      "x-velocity"
    },
    elemalias = {
      "material_indicator_numerical",
      "density_numerical",
      "pressure_numerical",
      "x-velocity_numerical"
    }
  }

}
