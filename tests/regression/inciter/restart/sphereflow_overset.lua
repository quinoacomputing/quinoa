-- This is a comment
-- Keywords are case-sensitive

inciter = {

  title = "Overset: flow over sphere",

  nstep = 1,
  ttyi = 1,
  cfl = 0.5,
  scheme = "oversetfe",

  partitioning = "rcb",

  compflow = {
    physics = "euler",
    problem = "user_defined"
  },

  material = {
    {
      gamma = { 1.66666666666667 }
    }
  },

  mesh = {
    -- depvars are automatically assigned and can be referenced
    -- downstream to request output variables
    {
      filename = "freestream_BGmesh_15k.exo",
    }, -- depvar: 'a'
    {
      filename = "sphere_OSmesh_12k.exo",
    } -- depvar: 'b' ...
  },

  ic = {
    density = 1.0,
    pressure = 1.0,
    velocity = { 1.0, 0.0, 0.0 }
  },

  bc = {
    {
      mesh = { 1 },
      farfield = { 1, 2, 3, 4, 5, 6 },
      density = 1.0,
      pressure = 1.0,
      velocity = { 1.0, 0.0, 0.0 }
    },

    {
      mesh = { 2 },
      symmetry = { 102 },
      farfield = { 101 },
      density = 1.0,
      pressure = 1.0,
      velocity = { 1.0, 0.0, 0.0 }
    }
  },

  diagnostics = {
    interval = 25,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 1,
    nodevar = {
      "density",
      "pressure"
    }
  }

}
