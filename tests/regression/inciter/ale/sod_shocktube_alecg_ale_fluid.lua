inciter = {

  title = "Sod shock-tube, ALECG, ALE",

  -- This ALE config for this mesh yields close to Lagrangian mesh velocity and
  -- should run until t=0.2 without a problem.

  nstep = 20,    -- Max number of time steps
  --term = 0.2,   -- Max physical time
  ttyi = 1,      -- TTY output interval
  cfl = 0.9,
  scheme = "alecg",

  partitioning = "rcb",

  compflow = {
    physics = "euler",
    problem = "sod_shocktube"
  },

  mesh = {
    { filename = "tube.exo" }
  },

  ale = {
    dvcfl = 1.0,
    mesh_velocity = "fluid",
    smoother = "laplace",
    vortmult = 0.0,      -- no vorticity scaling in Laplace smoother
    maxit = 20,
    tolerance = 2.0e-2,
    mesh_motion = { 0 }, -- mesh moves only in x
    dirichlet = { 1, 3 }
  },

  material = {
    {
      gamma = { 1.4 }
    }
  },

  bc = {
    {
      symmetry = { 2 }
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
