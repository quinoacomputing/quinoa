inciter = {

  title = "Water-air shock-tube",

  nstep = 25,
  cfl = 0.7,
  ttyi = 10,    -- TTY output interval
  scheme = "dgp1",
  limiter = "vertexbasedp1",
  shock_detector_coeff = 1.0,
  limsol_projection = false,
  lowspeed_kp = 1.0,

  partitioning = "mj",

  multimat = {
    physics = "euler",
    prelax = 1,
    prelax_timescale = 0.25,
    intsharp = 1,
    intsharp_param = 2.5,
    nmat = 2
  },

  material = {
    {
      id = { 1, 2 },
      gamma = { 4.4, 1.4 },  -- ratio of specific heats
      cv = { 951.36, 717.5 },  -- specific heat at const volume
      pstiff = { 6.0e8, 0.0 }  -- sg-eos stiffness parameter
    }
  },

  ic = {
    -- background (right-side conditions)
    materialid = 2,
    pressure = 1.0e5,
    temperature = 34.844,
    velocity = { 0.0, 0.0, 0.0 },

    -- left-side conditions
    box = {
      {
        materialid = 1,
        xmin = -1e-10, xmax = 0.75,
        ymin = -1.0, ymax = 1.0,
        zmin = -1.0, zmax = 1.0,
        pressure = 1.0e9,
        temperature = 494.646,
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
    interval = 1,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 25,
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
