# vim: filetype=sh:

title "Test evolution of some stats in mass fraction mixing"

walker

  nstep 10      # Max number of time steps
  #term  24.0    # Max time

  #dt    0.46         # Time step size, # A = 0.05
  #dt    0.1999999    # Time step size, # A = 0.5
  #dt    0.01        # Time step size, # A = 0.75
  dt 0.025

  npar  300000 #200000   # Number of particles
  ttyi  1      # TTY output interval

  rngs
    rngsse_mrg32k3a end
#    mkl_mcg59
#      beta_method cja
#    end
  end

  mixmassfracbeta
    depvar y
    ncomp 20
    init jointbeta
    solve fullvar

    icbeta
      # IC for A = 0.05, light = heavy: rho1=1.0, rho2=1.105263158, r=rho2/rho1-1=0.105263158
      #     <R>(t=0) = rhobar = 0.10521051556947E+01 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.479991014835299
      #     b(t=0) = -<rv>(t=0) = 0.21670652947218E-02 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 0.00250601493996136
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.215840181119751
      #     y2tnm = Yt(1-Yt) = 0.249599640512679
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.135254439163397
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.0750751648290741
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.0813343564092243
      # --------------------
      #     Vb = <V> = 0.95253507681071E+00 (from DNS)
      #     v2b = <v^2> = 0.19610183239403E-02 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.501618306963373
      #     y2b = (rho2/r)^2 * <v^2> = 0.216202269823195
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.135181861158038
      #     alpha = Yb*theta/(1-theta) = 0.0784091975881237
      #     beta = (1-Yb)*theta/(1-theta) = 0.0779032744641552
      #betapdf 0.0784091975881237 0.0779032744641552 0.0 1.0 end # A = 0.05, light = heavy

      # IC for A = 0.5, light = heavy: rho1=1.0, rho2=3.0, r=rho2/rho1-1=2.0
      #     <R>(t=0) = rhobar = 1.9899979572286 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.25377
      #     b(t=0) = -<rv>(t=0) = 0.28451153288492 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 0.33330
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.16165
      #     y2tnm = Yt(1-Yt) = 0.18937
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.14638
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.031709
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.093244
      # --------------------
      #     Vb = <V> = 0.64548384495498 (from DNS)
      #     v2b = <v^2> = 0.96076085517480E-01 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.468225767432470
      #     y2b = (rho2/r)^2 * <v^2> = 0.216171192414330
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.131809121857526
      #     alpha = Yb*theta/(1-theta) = 0.0535816391581342
      #     beta = (1-Yb)*theta/(1-theta) = 0.0608538380945397
      #betapdf 0.0535816391581342 0.0608538380945397 0.0 1.0 end # A = 0.5, light = heavy

      # IC for A = 0.75, light = heavy: rho1=1.0, rho2=7.0, r=rho2/rho1-1=6.0
      #     <R>(t=0) = rhobar = 0.39699938716878E+01 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.12720
      #     b(t=0) = -<rv>(t=0) = 0.10683383472286E+01 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 1.2856
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.092262
      #     y2tnm = Yt(1-Yt) = 0.11102
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.16896
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.025861
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.17745
      # --------------------
      #     Vb = <V> = 0.5209928312424 (from DNS)
      #     v2b = <v^2> = 0.15838400809351 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.44116
      #     y2b = (rho2/r)^2 * <v^2> = 0.21558
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.12557
      #     alpha = Yb*theta/(1-theta) = 0.063352
      #     beta = (1-Yb)*theta/(1-theta) = 0.080251
      betapdf 0.063352 0.080251 0.0 1.0 end # A = 0.75, light = heavy

      # IC for A = 0.05, light << heavy: rho1=1.0, rho2=1.105263158, r=rho2/rho1-1=0.105263158
      #     <R>(t=0) = rhobar = 0.10772925661975E+01 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.24666
      #     b(t=0) = -<rv>(t=0) = 0.16676130752122E-02 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 0.0019560
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.15842
      #     y2tnm = Yt(1-Yt) = 0.18582
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.14745
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.042660
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.13029
      # --------------------
      #     Vb = <V> = 0.92980091434993 (from DNS)
      #     v2b = <v^2> = 0.15034676534196E-02 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.26291
      #     y2b = (rho2/r)^2 * <v^2> = 0.16576
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.14463
      #     alpha = Yb*theta/(1-theta) = 0.044454
      #     beta = (1-Yb)*theta/(1-theta) = 0.12463
      #betapdf 0.044454 0.12463 0.0 1.0 end # A = 0.05, light << heavy

      # IC for A = 0.5, light << heavy: rho1=1.0, rho2=3.0, r=rho2/rho1-1=2.0
      #     <R>(t=0) = rhobar = 0.24685587563022E+01 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.107642008184135
      #     b(t=0) = -<rv>(t=0) = 0.21078901519139 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 0.260150897297512
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.0778293772633859
      #     y2tnm = Yt(1-Yt) = 0.0960552062582216
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.189743270612942
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.0252072533894947
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.208969475702874
      # --------------------
      #     Vb = <V> = 0.49048417911882E+00 (from DNS)
      #     v2b = <v^2> = 0.68422463030763E-01 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.235726268678230
      #     y2b = (rho2/r)^2 * <v^2> = 0.153950541819217
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.145475916611286
      #     alpha = Yb*theta/(1-theta) = 0.0401305190479043
      #     beta = (1-Yb)*theta/(1-theta) = 0.130111513258996
      #betapdf 0.0401305190479043 0.130111513258996 0.0 1.0 end # A = 0.5, light << heavy

      # IC for A = 0.75, light << heavy: rho1=1.0, rho2=7.0, r=rho2/rho1-1=6.0
      #     <R>(t=0) = rhobar = 0.54056762689193E+01 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.049156
      #     b(t=0) = -<rv>(t=0) = 0.76593537640565 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 1.0034
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.035677
      #     y2tnm = Yt(1-Yt) = 0.046740
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.23669
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.015242
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.29484
      # --------------------
      #     Vb = <V> = 0.32668167469715 (from DNS)
      #     v2b = <v^2> = 0.10620660817596 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.21446
      #     y2b = (rho2/r)^2 * <v^2> = 0.14456
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.14191
      #     alpha = Yb*theta/(1-theta) = 0.035467
      #     beta = (1-Yb)*theta/(1-theta) = 0.12991
      betapdf 0.035467 0.12991 0.0 1.0 end # A = 0.75, light << heavy

      # IC for A = 0.05, light >> heavy: rho1=1.0, rho2=1.105263158, r=rho2/rho1-1=0.105263158
      #     <R>(t=0) = rhobar = 0.10270889331314E+01 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.72307
      #     b(t=0) = -<rv>(t=0) = 0.16483632625279E-02 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 0.0019160
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.17227
      #     y2tnm = Yt(1-Yt) = 0.20024
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.13968
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.11740
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.044962
      # --------------------
      #     Vb = <V> = 0.97523041184457E+00 (from DNS)
      #     v2b = <v^2> = 0.14970707243051E-02 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.73992
      #     y2b = (rho2/r)^2 * <v^2> = 0.16505
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.14232
      #     alpha = Yb*theta/(1-theta) = 0.12278
      #     beta = (1-Yb)*theta/(1-theta) = 0.043157
      #betapdf 0.12278 0.043157 0.0 1.0 end # A = 0.05, light >> heavy

      # IC for A = 0.5, light >> heavy: rho1=1.0, rho2=3.0, r=rho2/rho1-1=2.0
      #     <R>(t=0) = rhobar = 0.15146897290206E+01 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.490301823047220
      #     b(t=0) = -<rv>(t=0) = 0.22418713440749E+00 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 0.254824646960634
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.219859807246818
      #     y2tnm = Yt(1-Yt) = 0.249905945363793
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.120229785142710
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.0670048631387328
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.0696555774084762
      # --------------------
      #     Vb = <V> = 0.80820983397097E+00 (from DNS)
      #     v2b = <v^2> = 0.78773502980312E-01 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.712314750956455
      #     y2b = (rho2/r)^2 * <v^2> = 0.177240381705702
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.135085566709960
      #     alpha = Yb*theta/(1-theta) = 0.111251978352117
      #     beta = (1-Yb)*theta/(1-theta) = 0.0449317567210851
      #betapdf 0.111251978352117 0.0449317567210851 0.0 1.0 end # A = 0.5, light >> heavy

      # IC for A = 0.75, light >> heavy: rho1=1.0, rho2=7.0, r=rho2/rho1-1=6.0
      #     <R>(t=0) = rhobar = 0.25440691870620E+01 (from DNS)
      #     Yt = <RY>/<R> = (rho2/<R>-1)/r = 0.29192
      #     b(t=0) = -<rv>(t=0) = 0.86607028308464 (from DNS)
      #     bnm = r^2/(1+r*Yt)^2 * Yt(1-Yt) = 0.98289
      #     y2t = <Ry"^2>/<R> = b*(r<R>/rho2)^{-2} = 0.18213
      #     y2tnm = Yt(1-Yt) = 0.20670
      #     thetaF = 1 - y2t/Yt/(1-Yt) = 0.11888
      #     alphaF = Yt*thetaF/(1-thetaF) = 0.039386
      #     betaF = (1-Yt)*thetaF/(1-thetaF) = 0.095534
      # --------------------
      #     Vb = <V> = 0.73349824469289 (from DNS)
      #     v2b = <v^2> = 0.13752428088617 (from DNS)
      #     Yb = <Y> = (rho2<V>-1)/r = 0.68908
      #     y2b = (rho2/r)^2 * <v^2> = 0.18719
      #     theta = 1 - y2b/Yb/(1-Yb) = 0.12630
      #     alpha = Yb*theta/(1-theta) = 0.099612
      #     beta = (1-Yb)*theta/(1-theta) = 0.044946
      betapdf 0.099612 0.044946 0.0 1.0 end # A = 0.75, light >> heavy

      betapdf 0.1167 0.2 0.0 1.0 end     # <r^3> = 0
      betapdf 0.1 0.1 0.0 1.0 end
    end
#    init jointdelta
#    icdelta
#      spike 0.01 0.5 0.99 0.5 end
#      spike 0.01 0.9 0.99 0.1 end
#      spike 0.01 0.1 0.99 0.9 end
#      spike 0.01 0.5 0.99 0.5 end
#      spike 0.01 0.5 0.99 0.5 end
#    end
    coeff hydrotimescale

    #hydrotimescales eq_A005S eq_A005H eq_A005L eq_A005S eq_A005S end # A = 0.05
    #hydrotimescales eq_A05S eq_A05H eq_A05L eq_A05S eq_A05S end # A = 0.5
    hydrotimescales eq_A075S eq_A075H eq_A075L eq_A075S eq_A075S end # A = 0.75

    #hydroproductions prod_A005S prod_A005H prod_A005L prod_A005S prod_A005S end # A = 0.05
    #hydroproductions prod_A05S prod_A05H prod_A05L prod_A05S prod_A05S end # A = 0.5
    hydroproductions prod_A075S prod_A075H prod_A075L prod_A075S prod_A075S end # A = 0.75

    #coefficients constraints for homdecay
    # alpha = Sb/kappa, beta = (1-S)b/kappa
    # decay in <y^2> if bprime/kappaprime > 1/4
    kappaprime 0.1    0.1      0.1   1.0  1.0 end
    bprime     1.0    1.0      1.0   1.0  1.0 end
    #S          2.2857 2.2857   2.2857  0.5  0.5 end
    #S          4.0 4.0 4.0  0.5  0.5 end
    #S          2.0 1.0 0.0  0.0  0.0 end
    S          1.0 1.0 4.0  0.0  0.0 end
    #rng mkl_mcg59
    rng rngsse_mrg32k3a

    #rho2 1.105263158 1.105263158 1.105263158 1.105263158 1.105263158 end # A = 0.05
    #rho2 3.0 3.0 3.0 3.0 3.0 end # A = 0.5
    rho2 7.0 7.0 7.0 7.0 7.0 end # A = 0.75

    #r 0.105263158 0.105263158 0.105263158 0.105263158 0.105263158 end # A = 0.05
    #r 2.0 2.0 2.0 2.0 2.0 end           # A = 0.5
    r 6.0 6.0 6.0 6.0 6.0 end           # A = 0.75
  end

  statistics
    interval  1
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
    # <rv^3>
    <y6y11y11y11> <y7y12y12y12> <y8y13y13y13> <y9y14y14y14> <y10y15y15y15>
    # <v^2>, specific volume variances
    <y11y11> <y12y12> <y13y13> <y14y14> <y15y15>
    # <v^3>, specific volume third moments
    <y11y11y11> <y12y12y12> <y13y13y13> <y14y14y14> <y15y15y15>
    # <r^2v>
    <y6y6y11> <y7y7y12> <y8y8y13> <y9y9y14> <y10y10y15>

#    # <R^2>
#    <Y6Y6> <Y7Y7> <Y8Y8> <Y9Y9> <Y10Y10>
#    # <R^3>
#    <Y6Y6Y6> <Y7Y7Y7> <Y8Y8Y8> <Y9Y9Y9> <Y10Y10Y10>

#    # <r^2v^2>
#    <y6y6y11y11> <y7y7y12y12> <y8y8y13y13> <y9y9y14y14> <y10y10y15y15>
#    # <Ry^3>
#    <Y6y1y1y1> <Y7y2y2y2> <Y8y3y3y3> <Y9y4y4y4> <Y10y5y5y5>
#    # <Ry^2>
#    <Y6y1y1> <Y7y2y2> <Y8y3y3> <Y9y4y4> <Y10y5y5>

#    <Y1y1> <y2y2> <y3y3> <y4y4> <y5y5>
#    <y11y11> <y12y12> <y13y13> <y14y14> <y15y15>
#    <Y6y1y1y1> <Y7y2y2y2> <Y8y3y3y3> <Y9y4y4y4> <Y10y5y5y5>
#    <Y6y1y1> <Y7y2y2> <Y8y3y3> <Y9y4y4> <Y10y5y5>
#    <Y6Y1> <Y7Y2> <Y8Y3> <Y9Y4> <Y10Y5>
#    <Y1> <Y2> <Y3> <Y4> <Y5>
#    <Y6> <Y7> <Y8> <Y9> <Y10>
#    <y6y6> <y7y7> <y8y8> <y9y9> <y10y10>
#    <y6y6y6> <y7y7y7> <y8y8y8> <y9y9y9> <y10y10y10>
#    <y6y11> <y7y12> <y8y13> <y9y14> <y10y15>
#    <Y6y11y11y11> <Y7y12y12y12> <Y8y13y13y13> <Y9y14y14y14> <Y10y15y15y15>

#    <y1y1y1y1>
#    <y2y2y2y2>
#    <y3y3y3y3>
#    <y4y4y4y4>
#    <y5y5y5y5>

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

#  pdfs
#    interval  1
#    filetype  txt
#    policy    multiple
#    centering elem
#    format    scientific
#    precision 6
##    p1( Y1 : 1.0e-2 )   # mass fraction PDFs
##    p2( Y2 : 1.0e-2 )
##    p3( Y3 : 1.0e-2 )
##    p4( Y4 : 1.0e-2 )
##    p5( Y5 : 1.0e-2 )
#    p6( Y6 : 1.0e-2 )   # density PDFs
#    p7( Y7 : 1.0e-2 )
#    p8( Y8 : 1.0e-2 )
#    p9( Y9 : 1.0e-2 )
#   p10( Y10 : 1.0e-2 )
#  end
end
