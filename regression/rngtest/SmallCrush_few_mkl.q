# vim: filetype=sh:
# This is a comment
# Keywords are case-sensitive

title "Example RNG test suite"

smallcrush

  mkl_mcg31 seed 0 end
  mkl_mcg59 seed 2 end

  rngsse_gm55
    seed 0
    seqlen long
  end
  rngsse_mrg32k3a seed 0 end

end
