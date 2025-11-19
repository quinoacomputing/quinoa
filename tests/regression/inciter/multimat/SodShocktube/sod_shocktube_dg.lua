inciter = {

  title = "Multi-material Sod shock-tube",

  nstep = 25,   -- Max number of time steps
  dt = 1.0e-3,  -- Time step size
  ttyi = 10,    -- TTY output interval
  scheme = "dgp0",
  limsol_projection = false,
  lowspeed_kp = 1.0,

  multimat = {
    physics = "euler",
    prelax = 0,
    nmat = 2
  },

  material = {
    {
      id = { 1, 2 },
      gamma = { 1.4, 1.4 },  -- ratio of specific heats
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
      symmetry = { 2, 4, 5, 6 }
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 5,
    elemvar = {
      "F1",
      "F2",
      "density",
      "pressure",
      "specific_total_energy",
      "x-velocity",
      "y-velocity",
      "z-velocity"
    }
  }

}
