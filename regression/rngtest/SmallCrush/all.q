# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Test all available RNGs with SmallCrush"

smallcrush

  rngsse_gm19 end
  rngsse_gm29 end
  rngsse_gm31 end
  rngsse_gm55 end
  rngsse_gm61 end
  rngsse_gq58.1 end
  rngsse_gq58.3 end
  rngsse_gq58.4 end
  rngsse_mt19937 end
  rngsse_lfsr113 end
  rngsse_mrg32k3a end

  mkl_mcg31 end
  mkl_r250 end
  mkl_mrg32k3a end
  mkl_mcg59 end
  mkl_wh end
  mkl_mt19937 end
  mkl_mt2203 end
  mkl_sfmt19937 end
  mkl_sobol end
  mkl_niederr end
  mkl_nondeterm end

end
