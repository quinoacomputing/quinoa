title "Test Ray's closure ideas for <y^2> and <rho v>"

walker

  #nstep 1      # Max number of time steps
  term  25.0    # Max time
  dt    0.002    # Time step size
  npar  500   # Number of particles
  ttyi  500      # TTY output interval

  rngs
    mkl_r250 end
  end

  massfracbeta
    depvar Y
    ncomp 15
    init zero
    coeff const
    # alpha = Sb/kappa, beta = (1-S)b/kappa
    kappa 2.0  0.76923  0.5  0.15873  0.1 end
    b     0.4  1.0      1.0  1.0    100.0 end
    S     0.5  0.53846  0.5  0.39683  0.5 end
    rng mkl_r250
    rho2 1.0 1.0 1.0 1.0 1.0 end
    #r 0.0101 0.0101 0.0101 0.0101 0.0101 end # low-A
    r 9.0 9.0 9.0 9.0 9.0 end   # high-A
  end

  statistics
    #precision 5
    #format    default
    # <Y>, mass fraction means
    <Y1>        # 3
    <Y2>        # 4
    <Y3>        # 5
    <Y4>        # 6
    <Y5>        # 7
    # <rho>, mean densities
    <Y6>        # 8
    <Y7>        # 9
    <Y8>        # 10
    <Y9>        # 11
    <Y10>       # 12
    # <V>, mean specific volumes
    <Y11>       # 13
    <Y12>       # 14
    <Y13>       # 15
    <Y14>       # 16
    <Y15>       # 17
    # <y^2>, mass fraction variances
    <y1y1>      # 23
    <y2y2>      # 24
    <y3y3>      # 25
    <y4y4>      # 26
    <y5y5>      # 27
     # <rho^2>, density variances
    <y6y6>      # 28
    <y7y7>      # 30
    <y8y8>      # 32
    <y9y9>      # 34
    <y10y10>    # 36
    # <v^2>, specific volume variances
    <y11y11>    # 38
    <y12y12>    # 39
    <y13y13>    # 40
    <y14y14>    # 41
    <y15y15>    # 42
     # <rho v>, density-specific-volume covariances
    <y6y11>     # 29
    <y7y12>     # 31
    <y8y13>     # 33
    <y9y14>     # 35
    <y10y15>    # 37
     # <rho v^2>
    <Y6y11y11>  # 18
    <Y7y12y12>  # 19
    <Y8y13y13>  # 20
    <Y9y14y14>  # 21
    <Y10y15y15> # 22
  end

end
