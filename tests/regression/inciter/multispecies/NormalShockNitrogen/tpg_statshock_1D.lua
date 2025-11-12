inciter = {

  title = "1D stationary normal shock (TPG)",

  nstep = 25,
  dt = 1.0e-7,
  --cfl = 0.5,
  ttyi = 25,  -- TTY output interval
  scheme = "p0p1",
  limiter = "vertexbasedp1",

  partitioning = "mj",

  multispecies = {
    nspec = 1,
  },

  material = {
    {
      eos    = "thermallyperfectgas",
    }
  },

  species = {
    {
      -- Values for Nitrogen (N2)
      id       = { 1 },
      R        = { 2.968030520448514e+02 },  -- specific gas constant
      cp_coeff = { { {2.210371497E+04, -3.818461820E+02, 6.082738360E+00, -8.530914410E-03, 1.384646189E-05, -9.625793620E-09, 2.519705809E-12, 7.108460860E+02},
                     {5.877124060E+05, -2.239249073E+03, 6.066949220E+00, -6.139685500E-04, 1.491806679E-07, -1.923105485E-11, 1.061954386E-15, 1.283210415E+04},
                     {8.310139160E+08, -6.420733540E+05, 2.020264635E+02, -3.065092046E-02, 2.486903333E-06, -9.705954110E-11, 1.437538881E-15, 4.938707040E+06} } }, -- Coefficients for cP / R (NASA GLENN polynomial)
      t_range  = { {200, 1000, 6000, 20000} }, -- Temperature range (K) over which the polynomials are valid
      dH_ref   = { 3.094984543111511e+05 } -- Reference enthalpy [J / kg]
    }
  },

  ic = {
    -- background (left-side conditions)
    mass_fractions = { 1.0 },
    pressure = 1.113011445168193e+05,
    temperature = 300,
    velocity = { 2.82474e+3, 0.0, 0.0 },

    -- right-side conditions
    box = {
      {
        mass_fractions = { 1.0 },
        xmin = 0.5, xmax = 1.0,
        ymin = -1.0, ymax = 1.0,
        zmin = -1.0, zmax = 1.0,
        pressure = 8.602273326645358e+06,
        temperature = 3.447516576090411e+03,
        velocity = { 4.200006676318108e+02, 0.0, 0.0 }
      }
    }
  },

  bc = {
    {
      extrapolate = { 3 },
      dirichlet = { 1 },
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
      "D1", "D2",
      "density",
      "specific_total_energy",
      "x-velocity",
      "y-velocity",
      "z-velocity"
    },
  }

}
