inciter = {

  title = "Solid-fluid interface advection from Favrie and Gavrilyuk JCP 2009",

  nstep = 10,
  cfl = 0.8,
  ttyi = 1,  -- TTY output interval
  scheme = "p0p1",
  limiter = "vertexbasedp1",
  flux = "laxfriedrichs",

  partitioning = "mj",

  multimat = {
    physics = "euler",
    problem = "user_defined",
    prelax = 0,
    intsharp = 0,
    nmat = 2
  },

  material = {
    -- Copper
    {
      eos = "smallshearsolid",
      id = { 1 },
      gamma = { 4.22 },  -- ratio of specific heats
      cv = { 3978.0 },  -- specific heat at const volume
      pstiff = { 342.0e8 },  -- sg-eos stiffness parameter
      mu = { 9.2e10 }  -- shear modulus
    },
    -- Air
    {
      id = { 2 },
      gamma = { 1.4 },
      cv = { 717.5 }
    }
  },

  ic = {
    -- background (right-side conditions)
    materialid = 2,
    pressure = 1.0e5,
    temperature = 300.0,
    velocity = { 100.0, 0.0, 0.0 },

    -- left-side conditions
    box = {
      {
        materialid = 1,
        xmin = -1e-10, xmax = 0.4,
        ymin = -1.0, ymax = 1.0,
        zmin = -1.0, zmax = 1.0,
        pressure = 1.0e5,
        temperature = 300.0,
        velocity = { 100.0, 0.0, 0.0 }
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
    interval = 5,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 5,
    elemvar = {
      "F1",
      "density",
      "pressure",
      "specific_total_energy",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "g_tensor"
    },
    elemalias = {
      "volfrac1",
      "density",
      "pressure",
      "total_energy_density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      ""
    }
  }

}
