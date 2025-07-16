-- This is a comment
-- Keywords are case-sensitive

inciter = {

  title = "Fluid-fluid interface advection",

  nstep = 125,
  cfl = 0.8,
  ttyi = 25,  -- TTY output interval
  scheme = "dg",
  flux = "hllc",

  partitioning = "mj",

  multimat = {
    physics = "euler",
    problem = "user_defined",
    prelax = 0,
    min_volumefrac = 1e-8,
    nmat = 2
  },

  material = {
    {
      eos = "stiffenedgas",
      id     = { 1,      2     },
      gamma  = { 4.4,    1.4   },  -- ratio of specific heats
      cv     = { 951.36, 717.5 },  -- specific heat at const volume
      pstiff = { 6.0e8,  0.0   }  -- sg-eos stiffness parameter
    }
  },

  ic = {
    -- background (right-side conditions)
    materialid = 2,
    pressure = 1.0e5,
    temperature = 300.0,
    velocity = { 1000.0, 0.0, 0.0 },

    -- left-side conditions
    box = {
      {
        materialid = 1,
        xmin = -1e-10, xmax = 0.5,
        ymin = -1.0, ymax = 1.0,
        zmin = -1.0, zmax = 1.0,
        pressure = 1.0e5,
        temperature = 300.0,
        velocity = { 1000.0, 0.0, 0.0 }
      }
    }
  },

  bc = {
    {
      symmetry = { 2, 4, 5, 6 },
      extrapolate = { 3 },
      farfield = { 1 },
      materialid = 1,
      pressure = 1.0e5,
      temperature = 300.0,
      velocity = { 1000.0, 0.0, 0.0 }
    }
  },

  diagnostics = {
    interval = 25,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 25,
    elemvar = {
      "F1",
      "density",
      "pressure",
      "specific_total_energy",
      "x-velocity",
      "y-velocity",
      "z-velocity"
    }
  }

}
