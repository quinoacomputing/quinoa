-- This is a comment
-- Keywords are case-sensitive

inciter = {

  title = "Sphere in a cube",

  nstep = 10,
  cfl = 0.8,
  ttyi = 1,  -- TTY output interval
  scheme = "oversetfe",

  partitioning = "mj",

  compflow = {
    problem = "user_defined"
  },

  rigid_body_motion = {
    rigid_body_dof = 3,
    symmetry_plane = {0, 0, 1},
  },

  mesh = {
    {
      filename = "halfcube_BGmesh_266k.exo",
    },
    {
      filename = "sphere_inhalfcube_OSmesh_11k.exo",
      mass = 0.01,
      center_of_mass = {-0.2, 0.2, 0.0},
      moment_of_inertia =
        { {0.0000000666667, 0, 0},
          {0, 0.0000000666667, 0},
          {0, 0, 0.0000000666667} }
    }
  },

  material = {
    -- ideal gas
    {
      id = { 1 },
      gamma = { 1.4 },
      cv = { 717.5 }
    }
  },

  -- units: m,kg,s,
  ic = {
    materialid = 1,
    density = 1.0,
    pressure = 1.0,
    velocity = { 0.0, 0.0, 0.0 }, -- m/s

    box = {
      {
        xmin = -1.0, xmax = -0.4,
        ymin =  0.1, ymax =  1.0,
        zmin = -0.5, zmax = 1e8,
        materialid = 1,
        density = 10.0,
        pressure = 1000.0,
      }
    }
  },

  bc = {
    {
      mesh = { 1 },
      symmetry = { 1 }
    },
    {
      mesh = { 2 },
      slipwall = { 101 },
      symmetry = { 103 },
      farfield = { 102 },
      density = 1.0,
      pressure = 1.0,
      velocity = { 0.0, 0.0, 0.0 } -- m/s
    }
  },

  diagnostics = {
    interval = 1,
    format = "scientific",
    error = "l2",
  },

  field_output = {
    time_interval = 5,
    nodevar = {
      "density",
      "specific_total_energy",
      "pressure"
    }
  }

}
