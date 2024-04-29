inciter = {

  title = "Solid-solid (Cu-Cu) shock tube",

  term = 3e-5,
  nstep = 10,
  cfl = 0.8,
  ttyi = 5,  -- TTY output interval
  scheme = "p0p1",
  limiter = "vertexbasedp1",
  flux = "laxfriedrichs",

  partitioning = "mj",

  multimat = {
    physics = "euler",
    problem = "user_defined",
    prelax = 0,
    intsharp = 1,
    nmat = 2
  },

  material = {
    -- Copper
    {
      eos = "smallshearsolid",
      id = { 1, 2 },
      gamma = { 4.22, 4.22 }, -- ratio of specific heats
      cv = { 3978.0, 3978.0 }, -- specific heat at const volume
      pstiff = { 342.0e8, 342.0e8 },  -- sg-eos stiffness parameter
      mu = { 9.2e10, 9.2e10 },  -- shear modulus
      rho0 = { 8900.0, 8900.0 }  -- initial density
    }
  },

  ic = {
    -- background (right-side conditions)
    materialid = 2,
    pressure = 1e5,
    temperature = 343.85505,
    velocity = { 0.0, 0.0, 0.0 },

    -- left-side conditions
    box = {
      {
        materialid = 1,
        xmin = -1e-10, xmax = 0.5,
        ymin = -1.0, ymax = 1.0,
        zmin = -1.0, zmax = 1.0,
        pressure = 5e9,
        temperature = 343.85505,
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
    interval = 2,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 10,
    elemvar = {
      "F1",
      "density",
      "pressure",
      "specific_total_energy",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "g_tensor"
    }
  }

}
