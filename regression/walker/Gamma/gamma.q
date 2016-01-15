title "Example problem"

walker

  term  35.0    # Max time
  dt    0.01    # Time step size
  npar  5000  # Number of particles
  ttyi  100     # TTY output interval

  rngs
    mkl_mcg59 seed 1 end
  end

  gamma
    depvar l
    init zero
    coeff const
    ncomp 2
    # k = bS/kappa, 1/theta = b(1-S)/kappa
    # <Y> = S/(1-S), <y^2> = kappa/b * <Y>/(1-S)
    b     1.5            2.5 end
    kappa 1.0            1.0 end
    S     0.666666666666 0.8 end
    rng mkl_mcg59
  end

  statistics
    <l1l1> <l2l2> <l1l2>
  end

end
