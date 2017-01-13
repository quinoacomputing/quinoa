# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Dirichlet for the IJSA paper using Random123's ThreeFry RNG"

walker
  term  140.0   # Max time
  dt    0.05    # Time step size
  npar  10000   # Number of particles
  ttyi  1000    # TTY output interval

  rngs
    r123_threefry end
  end

  dirichlet     # Select Dirichlet SDE
    depvar y
    init zero
    coeff const
    ncomp 2  # = K = N-1
    b     0.1    1.5 end
    S     0.625  0.4 end
    kappa 0.0125 0.3 end
    rng r123_threefry
  end

  statistics
    <Y1>
    <Y2>
    <y1y1>
    <y2y2>
    <y1y2>
  end
end
