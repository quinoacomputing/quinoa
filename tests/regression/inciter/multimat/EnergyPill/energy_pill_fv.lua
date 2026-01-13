inciter = {

  title = "Energy pill",

  nstep = 50,
  ttyi = 5,      -- TTY output interval
  cfl = 0.9,

  scheme = "fv",
  lowspeed_kp = 1.0,
  limiter = "vertexbasedp1",

  partitioning = "mj",

  multimat = {
    physics = "energy_pill",
    nmat = 2
  },

  material = {
    {
      id = { 1, 2 },
      gamma = { 1.4, 1.4 },
      cv = { 717.5, 717.5 }
    }
  },

  -- units: m,kg,s,K
  ic = {
    materialid = 1,
    temperature = 300.0,
    velocity = { 0.0, 0.0, 0.0 },
    pressure = 77.85e+3,

    meshblock = {
      {
        blockid = 2,
        materialid = 2,
        mass = 1024.0,
        volume = 0.64,
        energy_content = 9.0e+9,
        initiate = "impulse"
      }
    }
  },

  bc = {
    {
      symmetry = { 1, 3 }
    }
  },

  diagnostics = {
    interval = 5,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 10,
    sideset = { 1 },
    elemvar = {
      "density",
      "pressure",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy"
    }
  }

}
