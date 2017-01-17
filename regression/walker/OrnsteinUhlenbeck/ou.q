title "Example problem"

walker

  #nstep 1     # Max number of time steps
  term  5.0    # Max time
  dt    0.01   # Time step size
  npar  10000 # Number of particles
  ttyi  100    # TTY output interval

  rngs
    r123_threefry end
  end

  ornstein-uhlenbeck
    depvar r
    init raw
    coeff const
    ncomp 3
    theta 1.0 2.0 3.0 end
    mu 0.0 0.5 1.0 end
    sigmasq
      4.0  2.5   1.1
          32.0   5.6
                23.0
    end
    rng r123_threefry
  end

  statistics
    interval 2
    <R> <rr> <R2> <r2r2> <R3> <r3r3> <r1r2> <r1r3> <r2r3>
  end

end
