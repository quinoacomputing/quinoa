# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Example RNG test suite exercising RNGs from three RNG libraries"

smallcrush

  mkl_mcg31 seed 0 end

  rngsse_gm55
    seed 0
    seqlen long
  end

  r123_philox seed 21234 end
end
