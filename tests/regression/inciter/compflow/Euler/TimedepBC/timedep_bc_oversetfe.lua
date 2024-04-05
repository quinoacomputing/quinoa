inciter = {

  title = "Time dependent BC test",

  term = 0.5,     -- Max physical time, s
  ttyi = 10,      -- TTY output interval
  cfl = 0.9,
  scheme = "oversetfe",

  partitioning = "mj",

  compflow = {
    physics = "euler",
  },

  material = {
    {
      gamma = { 1.4 },
      cv = { 717.5 }
    }
  },

  ic = {
    density = 1.0,
    velocity = { 0.0, 0.0, 0.0 },
    pressure = 1.0
  },

  bc = {
    {
      timedep = {
        {
          sideset = { 1 },
          -- The pressure, density, and velocity components are specified here as
          -- a discrete function in time (tabular form). This table is sampled
          -- to apply time dep},ent boundary conditions on the above side sets.
          fn = {
             -- t    p     rho           u             v     w
            0.100,   1.0,  1.0        ,  0.0        ,  0.0,  0.0,
            0.101,   5.0,  2.818181818,  1.606438658,  0.0,  0.0,
            0.300,   1.0,  1.0        ,  0.0        ,  0.0,  0.0
          }
        }
      },
      farfield = { 2, 3, 4, 5, 6 },
      pressure = 1.0,
      density = 1.0,
      velocity = { 0.0, 0.0, 0.0 }
    }
  },

  field_output = {
    time_interval = 0.1,
    nodevar = {
      "density",
      "pressure",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy"
    }
  },

  diagnostics = {
    interval = 10,
    format = "scientific",
    error = "l2"
  }

}
