# vim: filetype=sh:

title "Test evolution of some stats in mass fraction mixing"

walker

  #nstep 1      # Max number of time steps
  term  5.0    # Max time
  dt    0.01    # Time step size
  npar  1000   # Number of particles
  ttyi  100      # TTY output interval

  rngs
    mkl_mcg59
      beta_method cja
    end
  end

  mixmassfracbeta
    depvar y
    ncomp 20
    init jointbeta
    icbeta
      betapdf 0.01 0.01 0.0 1.0 end       # light = heavy
      betapdf 0.2 0.8 0.0 1.0 end       # light > heavy
      betapdf 0.8 0.2 0.0 1.0 end       # light < heavy
      betapdf 0.116517612 0.2 0.0 1.0 end # <r^3> approximately 0
      betapdf 2.0 5.0 0.0 1.0 end
    end
#    init jointdelta
#    icdelta
#      spike 0.01 0.5 0.99 0.5 end
#      spike 0.01 0.9 0.99 0.1 end
#      spike 0.01 0.1 0.99 0.9 end
#      spike 0.01 0.5 0.99 0.5 end
#      spike 0.01 0.5 0.99 0.5 end
#    end
    coeff homdecay
    # alpha = Sb/kappa, beta = (1-S)b/kappa
    # decay in <y^2> if bprime/kappaprime > 1/4
    kappaprime 1.0  1.0  1.0  1.0  1.0 end
    bprime     1.9  1.9  1.9  1.9  1.9 end
    S          0.5  0.5  0.5  0.5  0.5 end
    rng mkl_mcg59
    rho2 1.0 1.0 1.0 1.0 1.0 end
    #r 2.0 2.0 2.0 2.0 2.0 end
    #r 0.0101 0.0101 0.0101 0.0101 0.0101 end    # low-A
    r 9.0 9.0 9.0 9.0 9.0 end                  # high-A
    #r 4.0 4.0 4.0 4.0 4.0 end
  end

  statistics
    format    scientific
    precision 12
    # <Y>, mass fraction means
    <Y1> <Y2> <Y3> <Y4> <Y5>
    # <R>, mean densities
    <Y6> <Y7> <Y8> <Y9> <Y10>
    # <V>, mean specific volumes
    <Y11> <Y12> <Y13> <Y14> <Y15>
    # <y^2>, mass fraction variances
    <y1y1> <y2y2> <y3y3> <y4y4> <y5y5>

    # stats required for homdecay
    # <y^3>, mass fraction third moments
    <y1y1y1> <y2y2y2> <y3y3y3> <y4y4y4> <y5y5y5>
    # <r^2>, density variances
    <y6y6> <y7y7> <y8y8> <y9y9> <y10y10>
    # <r^3>, density third moments
    <y6y6y6> <y7y7y7> <y8y8y8> <y9y9y9> <y10y10y10>
    # <rv>, density-specific-volume covariances
    <y6y11> <y7y12> <y8y13> <y9y14> <y10y15>
    # <RY>
    <Y6Y1> <Y7Y2> <Y8Y3> <Y9Y4> <Y10Y5>
    # <Rv^2>
    <Y6y11y11> <Y7y12y12> <Y8y13y13> <Y9y14y14> <Y10y15y15>
    # <rv^2>
    <y6y11y11> <y7y12y12> <y8y13y13> <y9y14y14> <y10y15y15>
    # <v^2>, specific volume variances
    <y11y11> <y12y12> <y13y13> <y14y14> <y15y15>

#    # additional stats required for mchomdecay
#    # <R^2>
#    <Y6Y6> <Y7Y7> <Y8Y8> <Y9Y9> <Y10Y10>
#    # <YR^2>
#    <Y1Y6Y6> <Y2Y7Y7> <Y3Y8Y8> <Y4Y9Y9> <Y5Y10Y10>
#    # <Y(1-Y)R^3>
#    <Y1Y16Y6Y6Y6>
#    <Y2Y17Y7Y7Y7>
#    <Y3Y18Y8Y8Y8>
#    <Y4Y19Y9Y9Y9>
#    <Y5Y20Y10Y10Y10>
  end

end
