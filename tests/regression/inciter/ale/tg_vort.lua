inciter = {

  title = "Taylor-Green with ALE, Laplacian smoothing, and vorticity scaling",

  nstep = 10,
  ttyi = 1,
  cfl = 0.5,
  scheme = "alecg",

  partitioning = "mj",

  ale = {
    dvcfl = 1.0,
    mesh_velocity = "fluid",
    smoother = "laplace",
    vortmult = 1.0,
    maxit = 10,
    tolerance = 1.0,
    dirichlet = { 1, 2, 3, 4, 5, 6 }
  },

  compflow = {
    physics = "euler",
    problem = "taylor_green"
  },

  -- if only a single filename is listen in a mesh ... },
  -- block, it is not really a coupling, only a way to define the
  -- input mesh filename in the control file, however if more
  -- than one filename is given, the first one is considered the
  -- background mesh (src) and the rest are destinations
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

  amr = {
    t0ref = true,
    initial = {"uniform"}
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
