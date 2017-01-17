# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Example RNG test suite exercising RNGs from three RNG libraries"

smallcrush

  mkl_mcg31 seed 0 end
  mkl_mcg59 seed 2 end

  rngsse_gm55
    seed 0
    seqlen long
  end
  rngsse_mrg32k3a seed 0 end

  r123_threefry seed 2 end
  r123_philox seed 21234 end
end
