title "Example problem"

walker

  term 10.0 # Max time
  dt 0.001 # Time step size
  npar 1000 # Number of particles
  ttyi 1000 # TTY output interval

  rngs
   rngsse_gq58.3 end
  end

  skew-normal
    depvar m
    init zero
    coeff const
    ncomp 2
    T 1.0 3.5 end
    sigmasq 0.04 0.25 end
    lambda 100.0 -50.0 end
    rng rngsse_gq58.3
  end

  statistics
    interval 2
    <m1m1> <m2m2>
  end

end
