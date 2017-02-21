title "Example problem"

walker

  term  5.0    # Max time
  dt    0.01   # Time step size
  npar  100000 # Number of particles
  ttyi  100    # TTY output interval

  rngs
    r123_threefry end
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
    rng r123_threefry
  end

  pdfs
    interval          500
    filetype          exodusii
    policy            overwrite
    centering         elem
    f2( R1 R2 : 2.0e-1 2.0e-1 ; -6.6 7 -12.4 13 )
    f3o( R1 R2 R3 : 5.0e-1 5.0e-1 5.0e-1 ; -6.5 7.5 -12.5 13 -8 10 )
    f3c( r1 r2 r3 : 5.0e-1 5.0e-1 5.0e-1 ; -6.5 7.5 -13 12.5 -9 9 )
  end
end
