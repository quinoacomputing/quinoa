title "Example problem"

walker

  term  10.0    # Max time
  dt    0.001   # Time step size
  npar  1000   # Number of particles
  ttyi  1000    # TTY output interval

  rngs
   rngsse_mrg32k3a end
  end

  diag_ou
    depvar o
    init raw
    coeff const
    ncomp 2
    sigmasq 0.25 1.0 end
    theta 1.0 1.0 end
    mu 0.0 1.5 end
    rng rngsse_mrg32k3a
  end

  statistics
   interval 2
   <o1o1> <o2o2> <o1o2>
  end

end
