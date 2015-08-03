title "Test Ray's closure ideas for <v^2> and <v^3>
       for the number-fraction beta SDE"

walker

  term  25.0    # Max time
  dt    0.002   # Time step size
  npar  100   # Number of particles
  ttyi  100     # TTY output interval

  rngs
    mkl_r250 end
  end

  numfracbeta   # Select the number-fraction beta SDE
    depvar x    # Dependent variable: X, denoting mole or number fractions
    ncomp 15    # = 3x5 = 5 systems each computing 3 variables:
                #   X - number fraction,
                #   R - density,
                #   V - specific volume
    init zero
    coeff const
    # alpha = Sb/kappa, beta = (1-S)b/kappa
    kappa 2.0  0.76923  0.5  0.15873  0.1 end
    b     0.4  1.0      1.0  1.0    100.0 end
    S     0.5  0.53846  0.5  0.39683  0.5 end
    rng mkl_r250
    rho2 1.0 1.0 1.0 1.0 1.0 end
    # low-A
    rcomma 1.0e-2 1.0e-2 1.0e-2 1.0e-2 1.0e-2 end
    # high-A
    #rcomma 0.9 0.9 0.9 0.9 0.9 end
  end

  statistics    # column numbers in output file
    # <X>, mole fraction means
    <X1>        # 3
    <X2>        # 4
    <X3>        # 5
    <X4>        # 6
    <X5>        # 7
    # <rho>, mean densities
    <X6>        # 8
    <X7>        # 9
    <X8>        # 10
    <X9>        # 11
    <X10>       # 12
    # <V>, mean specific volumes
    <X11>       # 13
    <X12>       # 14
    <X13>       # 15
    <X14>       # 16
    <X15>       # 17
    # <x^2>, mole fraction variances
    <x1x1>      # 18
    <x2x2>      # 19
    <x3x3>      # 20
    <x4x4>      # 21
    <x5x5>      # 22
     # <rho v>, density-specific-volume covariances
    <x6x11>     # 25
    <x7x12>     # 30
    <x8x13>     # 25
    <x9x14>     # 40
    <x10x15>    # 45
    # <rho v^2>
    <x6x11x11>  # 26
    <x7x12x12>  # 31
    <x8x13x13>  # 36
    <x9x14x14>  # 41
    <x10x15x15> # 46
    # <rho^2 v>
    <x6x6x11>   # 23
    <x7x7x12>   # 28
    <x8x8x13>   # 33
    <x9x9x14>   # 38
    <x10x10x15> # 43
    # <rho^2v^2>
    <x6x6x11x11>   # 24
    <x7x7x12x12>   # 29
    <x8x8x13x13>   # 34
    <x9x9x14x14>   # 39
    <x10x10x15x15> # 43
    # <rho v^3>
    <x6x11x11x11>  # 27
    <x7x12x12x12>  # 32
    <x8x13x13x13>  # 37
    <x9x14x14x14>  # 42
    <x10x15x15x15> # 47
    # <v^2>, specific volume variances
    <x11x11>    # 48
    <x12x12>    # 50
    <x13x13>    # 52
    <x14x14>    # 54
    <x15x15>    # 56
    # <v^3>, specific volume third central moments
    <x11x11x11> # 49
    <x12x12x12> # 51
    <x13x13x13> # 53
    <x14x14x14> # 55
    <x15x15x15> # 57
  end

  pdfs
    interval  100
    filetype  txt
    policy    overwrite
    centering elem
    format    scientific
    precision 4
    p1( X1 : 1.0e-2 )
    p2( X2 : 1.0e-2 )
    p3( X3 : 1.0e-2 )
    p4( X4 : 1.0e-2 )
    p5( X5 : 1.0e-2 )
  end
end
