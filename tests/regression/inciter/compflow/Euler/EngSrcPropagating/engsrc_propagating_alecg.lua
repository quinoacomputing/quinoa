inciter = {

  title = "Blast tube",

  nstep = 100,
  ttyi = 5,      -- TTY output interval
  cfl = 0.9,
  scheme = "alecg",

  partitioning = "mj",

  compflow = {
    physics = "euler"
  },

  -- units: m,kg,s,
  ic = {
    density = 0.912, -- kg/m^3
    velocity = { 0.0, 0.0, 0.0 }, -- m/s
    pressure = 77.85e+3, -- Pa = N/m^2 = kg/m/s^2

    box = {
      {
        xmin = -0.200005, xmax = 1.0e-14,
        ymin = -0.200005, ymax = 0.200005,
        zmin = -0.000005, zmax = 6.000005,
        mass = 768.48,
        energy_content = 9.0e+9,

        -- linear propagation of energy source
        initiate = "linear",
        point = { 1.0e-14, 0.0, 6.0 },
        front_speed = 8.2e+3, -- detonation velocity: 0.82 cm/us = 82e+3 m/s
        front_width = 0.5
      }
    }
  },

  material = {
    {
      gamma = { 1.4 },
      cv = { 717.5 }
    }
  },

  bc = {
    {
      symmetry = { 1, 3 }
    }
  },

  field_output = {
    interval = 25,
    nodevar = {
      "density",
      "specific_total_energy",
      "pressure"
    },
    nodealias = {
      "density_numerical",
      "specific_total_energy_numerical",
      "pressure_numerical"
    }
  },

  diagnostics = {
    interval = 5,
    format = "scientific",
    error = "l2"
  }

}
