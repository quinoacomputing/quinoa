title "Example problem"

walker

  #nstep 1     # Max number of time steps
  term  5.0    # Max time
  dt    0.01   # Time step size
  npar  10000 # Number of particles
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

  statistics
    interval 2
    <R> <rr> <R2> <r2r2> <R3> <r3r3> <r1r2> <r1r3> <r2r3>
  end

  pdfs
    interval          100
    filetype          gmshbin
    policy            overwrite
    centering         node
    #format            scientific
    #precision         4
    f2( r1 r2 : 2.0e-1 2.0e-1 ) #; -2 2 -2 2 )
    f3( r1 r2 r3 : 5.0e-1 5.0e-1 5.0e-1 ) #; 0 1 0 1 -0.5 0.5 )
  end
end
