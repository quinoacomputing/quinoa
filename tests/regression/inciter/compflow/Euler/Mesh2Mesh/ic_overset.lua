inciter = {

  title = "Coupled Euler equations",

  nstep = 1,
  ttyi = 1,
  cfl = 0.5,
  scheme = "oversetfe",

  partitioning = "rcb",

  compflow = {
    physics = "euler",
    problem = "rayleigh_taylor",
    alpha = 1.0,
    betax = 1.0,
    betay = 1.0,
    betaz = 1.0,
    p0 = 2.0,
    r0 = 2.0,
    kappa = 1.0
  },

  mesh = {
    -- depvars are automatically assigned and can be referenced
    -- downstream to request output variables
    { filename = "freestream_BGmesh_15k.exo" }, -- depvar: 'a'
    { filename = "sphere_OSmesh_12k.exo" } -- depvar: 'b' ...
  },

  material = {
    {
      gamma = { 1.66666666666667 }
    }
  },

  ic = {
    density = 1.0,
    velocity = { 0.0, 0.0, 0.0 },
    pressure = 10.0
  },

  bc = {
    {
      mesh = { 1 },
      dirichlet = { 1, 2, 3, 4, 5, 6 }
    },
    {
      mesh = { 2 },
      farfield = { 101 },
      pressure = 10.0,
      density = 1.0,
      velocity = { 0.0, 0.0, 0.0 }
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 5,
    nodevar = {
      "A1", "A2", "A3", "A4", "A5",
      "B1", "B2", "B3", "B4", "B5"
    }
  }

}
