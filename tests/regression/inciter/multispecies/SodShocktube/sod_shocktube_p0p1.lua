inciter = {

  title = "Multi-species Sod shock tube problem",

  nstep = 25,
  dt = 1.0e-4,
  ttyi = 25,  -- TTY output interval
  scheme = "p0p1",
  limiter = "vertexbasedp1",

  partitioning = "mj",

  multispecies = {
    nspec = 2,
  },

  material = {
    {
      eos    = "thermallyperfectgas",
    }
  },

  species = {
    {
      id       = { 1, 2 },
      R        = { 2.870025066673538e+02, 2.870025066673538e+02 },  -- specific gas constant
      cp_coeff = { { {0, 0, 3.5, 0, 0, 0, 0, 0}, {0, 0, 3.5, 0, 0, 0, 0, 0}, {0, 0, 3.5, 0, 0, 0, 0, 0} },
                 { {0, 0, 3.5, 0, 0, 0, 0, 0}, {0, 0, 3.5, 0, 0, 0, 0, 0}, {0, 0, 3.5, 0, 0, 0, 0, 0} } }, -- Coefficients for cP / R (NASA GLENN polynomial)
      t_range  = { {1e-8, 1000, 6000, 20000}, {1e-8, 1000, 6000, 20000} }, -- Temperature range (K) over which the polynomials are valid
      dH_ref   = { 0, 0 } -- Reference enthalpy, h(t = 298.15 K) - h(t = 0 K)
    }
  },

  ic = {
    -- background (right-side conditions)
    pressure = 0.1,
    temperature = 2.787456446e-3,
    velocity = { 0.0, 0.0, 0.0 },
    mass_fractions = { 0.0, 1.0 },

    -- left-side conditions
    box = {
      {
        mass_fractions = { 1.0, 0.0 },
        xmin = -1e-10, xmax = 0.5,
        ymin = -1.0, ymax = 1.0,
        zmin = -1.0, zmax = 1.0,
        pressure = 1.0,
        temperature = 0.0034843206,
        velocity = { 0.0, 0.0, 0.0 }
      }
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
      "D1", "D2",
      "density",
      "specific_total_energy",
      "x-velocity",
      "y-velocity",
      "z-velocity"
    },
  }

}
