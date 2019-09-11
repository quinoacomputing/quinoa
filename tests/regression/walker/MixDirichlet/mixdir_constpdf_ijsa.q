# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Mixture Dirichlet holding Dirichlet IC constant in time"

walker
  nstep 300
  #term  1.0   # Max time
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
        5.0 2.0 3.0     # IJSA paper
        #0.011834   0.035503   0.106509
      end
    end

    coeff const_coeff

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
    <Y1> <Y2>
    # mass fractions covariance matrix
    <y1y1> <y2y2> <y1y2>
  end

 pdfs
   interval  1000
   filetype  txt
   policy    overwrite
   centering elem
   format    scientific
   precision 8
   p1( Y1 : 1.0e-2 ; 0.01 0.99 )
   #p2( Y2 : 1.0e-2 ; 0.01 0.99 )
  end
end
