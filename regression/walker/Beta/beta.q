title "Example problem"

walker

  #nstep 1      # Max number of time steps
  term  50.0    # Max time
  dt    0.005   # Time step size
  npar  1000   # Number of particles
  ttyi  1000    # TTY output interval

  rngs
    rngsse_gm31 end
  end

  beta
    depvar x
    ncomp 5

    init jointdelta
    icdelta
      spike 0.01 0.5 0.99 0.5 end
      spike 0.01 0.9 0.99 0.1 end
      spike 0.01 0.1 0.99 0.9 end
      spike 0.01 0.5 0.99 0.5 end
      spike 0.01 0.5 0.99 0.5 end
    end

    coeff const
    # alpha = Sb/kappa, beta = (1-S)b/kappa
    # S = 1/(1+\beta/alpha), delta = S/alpha = kappa/b
    kappa 2.0  0.76923  0.5  0.15873  0.5 end
    b     0.4  1.0      1.0  1.0      8.0 end
    S     0.5  0.53846  0.5  0.39683  0.5 end
    rng rngsse_gm31
  end

  statistics
    <X1> <X2> <X3> <X4> <X5>
    <x1x1> <x1x2> <x1x3> <x1x4> <x1x5>
           <x2x2> <x2x3> <x2x4> <x2x5>
                  <x3x3> <x3x4> <x3x5>
                         <x4x4> <x4x5>
                                <x5x5>
  end

end
