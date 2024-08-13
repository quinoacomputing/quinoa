inciter = {

  title = "Sod shock-tube",

  nstep = 10,    -- Max number of time steps
  term = 0.2,    -- Max physical time
  ttyi = 1,      -- TTY output interval
  cfl = 0.5,
  scheme = "alecg",

  partitioning = "mj",

  compflow = {
    physics = "euler",
  },

  material = {
    {
      gamma = { 1.4 }
    }
  },

  bc = {
    {
      symmetry = { 2, 4, 5, 6 }
    }
  },

  -- userdef problem setup for sod_shocktube
  ic = {
    density = -1.0,                      -- overwritten by boxes
    velocity = { 100.0, 100.0, 100.0 },  -- overwritten by boxes
    pressure = 10000.0,                  -- overwritten by boxes
    box = {
      {
        xmin = -0.5, xmax = 0.5,
        ymin = -0.5, ymax = 0.5,
        zmin = -0.5, zmax = 0.5,
        density = 1.0,
        pressure = 1.0
      },
      {
        xmin =  0.5, xmax = 1.5,
        ymin = -0.5, ymax = 0.5,
        zmin = -0.5, zmax = 0.5,
        density = 0.125,
        pressure = 0.1
      }
    }
  },

  field_output = {
    interval = 10000,
    nodevar = {
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
   },
    nodealias = {
      "density_numerical",
      "x-velocity_numerical",
      "y-velocity_numerical",
      "z-velocity_numerical",
      "specific_total_energy_numerical",
      "pressure_numerical"
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  },

}
