# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Langevin velocity PDF model keeping the joint
       PDF constat in time using constant coefficients"

walker

  #nstep 2     # Max number of time steps
  term  1.0    # Max time
  dt    0.01   # Time step size
  npar  50000 # Number of particles
  ttyi  10    # TTY output interval

  rngs
    r123_philox end
  end

  velocity      # configure a velocity equation
    depvar u
    solve fullvar # fluctuation
    init jointgaussian
    icgaussian # unit kinetic energy, isotropic Reynolds stress at t=0
      gaussian 0.0 0.666667 end
      gaussian 0.0 0.666667 end
      gaussian 0.0 0.666667 end
    end
    coeff stationary
    rng r123_philox
  end

  statistics
    interval 1
    # estimate and save mean velocity
    <U1> <U2> <U3>
    # estimate and save Reynolds stress
    <u1u1> <u2u2> <u3u3> <u1u2> <u1u3> <u2u3>
  end

  pdfs
    interval  100
    filetype  txt
    policy    overwrite #multiple
    centering elem
    format    scientific
    precision 12
    # save marginal PDFs of all velocity components
    u1( u1 : 5.0e-2 ; -4 4 )
    u2( u2 : 5.0e-2 ; -4 4 )
    u3( u3 : 5.0e-2 ; -4 4 )
  end
end
