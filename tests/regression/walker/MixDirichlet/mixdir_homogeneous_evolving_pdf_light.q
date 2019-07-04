# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Mixture Dirichlet forcing d<R>/dt=0 constaining S"

walker
  nstep 300
  #term  10.0   # Max time
  dt    0.001    # Time step size
  npar  50000   # Number of particles
  ttyi  100     # TTY output interval

  rngs
    r123_threefry end
  end

  mixdirichlet
    depvar y
    ncomp 5  # = N+2 = K+3

    init jointdirichlet
    icdirichlet
      dirichletpdf
        0.011834   0.035503   0.106509
        #5.0 2.0 3.0     # IJSA paper
      end
    end

    coeff homogeneous

    # IJSA paper
    b     0.1    1.5 end
    S     0.625  0.4 end
    kappaprime 0.0125 0.3 end

    #kappaprime
    #  1.0   1.0
    #end
    #b
    #  0.11834   0.14201
    #end
    #S
    #  0.1   0.25
    #end

    normalization light
    # rho in arbitrary order, will be sorted
    rho 1.0 3.0 9.0 end # need N = K+1

    rng r123_threefry
  end

  statistics
    interval  1
    precision 8 #max
    # mean mass fractions
    <Y1> <Y2> <Y3>
    # mean density
    <Y4>
    # mean specific volume
    <Y5>
    # mass fractions covariance matrix
    <y1y1> <y2y2> <y1y2>
    # density variance
    <y4y4>
    # density third moment
    <y4y4y4>
    # density-specific-volume covariance
    <y4y5>                     # <rv> for b = -<rv>
    <y4y1> <y4y2> <y4y3>       # <ry.> for b. = -<ry.>/<R>
    # other statistics, needed for the model
    <Y4Y1> <Y4Y2> <Y4Y3>       # <RY.>
    <Y4Y4>                     # <R^2>
    <Y4Y4Y1> <Y4Y4Y2> <Y4Y4Y3> # <R^2Y.>
    <Y4Y4Y4Y1> <Y4Y4Y4Y2>      # <R^3Y.>
    <Y4Y4Y4Y1Y1> <Y4Y4Y4Y1Y2> <Y4Y4Y4Y1Y3> # <R^3YiYj>
    <Y4Y4Y4Y2Y1> <Y4Y4Y4Y2Y2> <Y4Y4Y4Y2Y3>
    <Y4Y4Y4Y3Y1> <Y4Y4Y4Y3Y2> <Y4Y4Y4Y3Y3>
    <y4y4y1> <y4y4y2>          # <r^2y>
    <y4y4y5> <y4y4y5>          # <r^2v>
  end

 #pdfs
 #  interval  1000
 #  filetype  txt
 #  policy    overwrite
 #  centering elem
 #  format    scientific
 #  precision 8
 #  p1( Y1 : 1.0e-2 ; 0.01 0.99 )
 #  #p2( Y2 : 1.0e-2 ; 0.01 0.99 )
 #end
end
