# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Test all available MKL RNGs with SmallCrush"

smallcrush

  mkl_mcg31 end
  #mkl_r250 end         # no leapfrog support
  #mkl_mrg32k3a end     # no leapfrog support
  mkl_mcg59 end
  mkl_wh end
  #mkl_mt19937 end      # no leapfrog support
  #mkl_mt2203 end       # no leapfrog support
  #mkl_sfmt19937 end    # no leapfrog support
  #mkl_sobol end        # no leapfrog support
  #mkl_niederr end      # no leapfrog support
  #mkl_nondeterm end    # no leapfrog support

end
