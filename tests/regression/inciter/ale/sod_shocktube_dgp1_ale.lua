inciter = {

  title = "Sod shock-tube, DG, ALE",

  -- This exercises the DG-ALE path with mesh_velocity = "fluid", i.e. the mesh
  -- velocity is computed from the fluid velocity in the elements surrounding
  -- each node, which include chare-boundary ghost elements. Run on more than a
  -- single mesh partition this therefore tests that the communicated solution
  -- ghost data are applied before the mesh velocity is evaluated.

  --term = 0.2,
  nstep = 100,   -- Max number of time steps
  cfl = 0.8,
  ttyi = 10,     -- TTY output interval
  scheme = "dgp1",
  limiter = "vertexbasedp1",
  flux = "ausm",
  lowspeed_ku = 0.0,
  lowspeed_kp = 0.0,

  partitioning = "mj",

  ale = {
    mesh_velocity = "fluid",
    symmetry = { 2, 4, 5, 6 },
    dvcfl = 1.0
  },

  multimat = {
    physics = "euler",
    nmat = 2,
    prelax = 1,
    intsharp = 1,
    intsharp_param = 2.5
  },

  material = {
    {
      eos = "stiffenedgas",
      id = { 1, 2 },
      gamma = { 1.4, 1.4 }  -- ratio of specific heats
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
    interval = 10,
    format = "scientific",
    error = "l2"
  },

  field_output = {
    interval = 50,
    elemvar = {
      "material_indicator",
      "density",       -- bulk density
      "pressure",      -- bulk pressure
      "x-velocity",
      "y-velocity",
      "z-velocity"
    }
  }

}
