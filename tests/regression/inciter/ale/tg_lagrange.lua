inciter = {

  title = "Taylor-Green with pure Lagrangian mesh motion",

  nstep = 10,
  ttyi = 1,
  cfl = 0.5,
  scheme = "alecg",

  partitioning = "mj",

  ale = {
    dvcfl = 1.0,
    mesh_velocity = "fluid"
  },

  compflow = {
    physics = "euler",
    problem = "taylor_green"
  },

  mesh = {
    { filename = "unitcube_1k.exo" }
  },

  material = {
    {
      gamma = { 1.4 }
    }
  },

  bc = {
    {
      dirichlet = { 1, 2, 3, 4, 5, 6 }
    }
  },

  field_output = {
    interval = 10,
    nodevar = {
      "density",
      "x-velocity",
      "y-velocity",
      "z-velocity",
      "specific_total_energy",
      "pressure"
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2"
  }

}
