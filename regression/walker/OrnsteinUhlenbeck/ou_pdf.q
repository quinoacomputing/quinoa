title "Example problem"

walker

  term  5.0    # Max time
  dt    0.01   # Time step size
  npar  100000 # Number of particles
  ttyi  100    # TTY output interval

  rngs
    mkl_mrg32k3a seed 0 end
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
    rng mkl_mrg32k3a
  end

  pdfs
    interval          500
    filetype          txt
    policy            overwrite
    centering         elem
    f1( r1 : 2.0e-1 ; -6.0 6.0 )
    f2( R1 R2 : 2.0e-1 2.0e-1 )                 # output but unverified
    f3o( R1 R2 R3 : 5.0e-1 5.0e-1 5.0e-1 )      # output but unverified
    f3c( r1 r2 r3 : 5.0e-1 5.0e-1 5.0e-1 )      # output but unverified
  end
end
