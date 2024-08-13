inciter = {

  title = "Moving overset mesh",

  nstep = 50,
  ttyi = 5,
  cfl = 0.75,
  scheme = "oversetfe",

  partitioning = "rcb",

  compflow = {
    physics = "euler",
    problem = "user_defined"
  },

  mesh = {
    -- depvars are automatically assigned and can be referenced
    -- downstream to request output variables
    { filename = "freestream_BGmesh_15k.exo", velocity = { 0, 0, 0 } }, -- depvar: 'a'
    { filename = "sphere_OSmesh_12k.exo", velocity = { 5, 0, 0 } } -- depvar: 'b' ...
  },

  ic = {
    density = 1.0,
    velocity = { 0.0, 0.0, 0.0 },
    pressure = 1.0,

    box = {
      {
      xmin = -0.5, xmax = 0.525,
      ymin = -0.5, ymax = 1,
      zmin = -0.5, zmax = 1,

      density = 10.0,
      velocity = { 0, 0, 0 },
      pressure = 50.0
      }
    }
  },

  material = {
    {
      gamma = { 1.66666666666667 }
    }
  },

  bc = {
    {
      mesh = { 1 },
      dirichlet = { 4, 6 },
      symmetry = { 1, 2, 3, 5 }
    },

    {
      mesh = { 2 },
      symmetry = { 102 },
      farfield = { 101 },
      pressure = 1.0,
      density = 1.0,
      velocity = { 0.0, 0.0, 0.0 }
    }
  },

  diagnostics = {
    interval = 10,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    nodevar = {
      "A1", "A2", "A3", "A4", "A5",
      "B1", "B2", "B3", "B4", "B5"
    },
    interval = 25
  },

  history_output = {
    interval = 10,
    point = {
      {id = "p1", coord = { 0.5, 0.25, 0.25 }}
    }
  }

}
