-- This is a comment
-- Keywords are case-sensitive

inciter = {

  title = "Sod shock-tube",

  nstep = 5,  -- Max number of time steps
  cfl = 0.8,
  ttyi = 1,  -- TTY output interval
  scheme = "fv",
  limiter = "vertexbasedp1",
  flux = "ausm",

  partitioning = "mj",

  multimat = {
    physics = "euler",
    prelax = 0,
    intsharp = 1,
    intsharp_param = 2.5,
    nmat = 2
  },

  material = {
    {
      eos = "stiffenedgas",
      id = { 1, 2 },
      gamma = { 1.4, 1.4 }  -- ratio of specific heats
    }
  },

  ic = {
    -- background state is left state
    materialid = 1,
    temperature = 3.484321e-3,
    pressure = 1.0,
    velocity = { 0, 0, 0 },

    -- right state
    box = {
      {
        materialid = 2,
        xmin = 0.5, xmax = 1.1,
        ymin = -0.5, ymax = 0.5,
        zmin = -0.5, zmax = 0.5,
        temperature = 2.7874568e-3,
        pressure = 0.1,
        velocity = { 0, 0, 0 }
      }
    }
  },

  bc = {
    {
      extrapolate = { 1, 3 },
      symmetry = { 2, 4, 5, 6 },
    }
  },

  field_output = {
    interval = 1,
    elemvar = {
      "F1",
      "density",
      "pressure",
      "x-velocity"
    }
  }

}
