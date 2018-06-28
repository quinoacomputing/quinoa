# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Coupled position, velocity, dissipation joint PDF computing a
       homogeneous shear flow using the simplifed Langevin model"

walker

  #nstep 2     # Max number of time steps
  term  15.0    # Max time
  dt    0.2   # Time step size
  npar  50000 # Number of particles
  ttyi  1    # TTY output interval

  rngs
    r123_philox end
  end

  position      # configure a position equation
    depvar x
    solve fluctuation
    velocity u  # couple a velocity model with dependent variable u
    #init zero
    init jointgaussian
    icgaussian
      gaussian 0.0 1.0 end
      gaussian 0.0 1.0 end
      gaussian 0.0 1.0 end
    end
    coeff const_shear
  end

  velocity      # configure a velocity equation
    depvar u
    solve fluctuation
    variant glm # select the simplified Langevin model (SLM/GLM)
    position x  # couple a position model with dependent variable x
    dissipation o  # couple a dissipation model with dependent variable o
    init jointgaussian
    icgaussian # unit kinetic energy, isotropic Reynolds stress at t=0
      gaussian 0.0 0.666667 end
      gaussian 0.0 0.666667 end
      gaussian 0.0 0.666667 end
    end
    #C0 2.1     # optional
    coeff const_shear
    rng r123_philox
  end

  dissipation
    depvar o
    velocity u  # couple a velocity model with dependent variable u
    init jointgamma
    icgamma
      gammapdf 4.0 0.25 end  # mean = 1.0, variance = 0.25
    end
    coeff const
    #C3 1.0     # optional
    #C4 0.25    # optional
    #COM1 0.44  # optional
    #COM2 0.9   # optional
    rng r123_philox
  end

  statistics
    interval 1
    # estimate and save mean velocity
    #<U1> <U2> <U3>
    # estimate and save Reynolds stress
    <U1U1> <U2U2> <U3U3> <U1U2> <U1U3> <U2U3>
    # estimate and save mean turbulence frequency
    <O>
  end

  pdfs
    interval  10
    filetype  txt
    policy    multiple
    centering elem
    format    scientific
    precision 6
    # save marginal PDFs of all velocity components
    U1( U1 : 1.0e-2 )
    U2( U2 : 1.0e-2 )
    U3( U3 : 1.0e-2 )
    O( O : 1.0e-2 )     # dissipation (turbulence frequency) PDF
  end
end
