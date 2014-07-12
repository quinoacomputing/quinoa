//******************************************************************************
/*!
  \file      src/Control/RNGSSEKeywords.h
  \author    J. Bakosi
  \date      Fri 27 Dec 2013 08:16:32 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     RNG keywords for RNGSSE lib
  \details   Random number generator selector keywords for those generators
  available in the RNGSSE library's. The keywords are defined by specializing
  struct 'keyword', defined in Control/Keyword.h. Introducing a new keyword
  requires a more human readable (but still short) name as well as a short,
  few-line, help-like description.
*/
//******************************************************************************
#ifndef Keywords
#error "RNGSSEKeywords.h should only be included within a *Keywords.h"
#endif

namespace tk {
namespace kw {

using namespace pegtl::ascii;

// Keyword 'rngsse_gm19'
struct rngsse_gm19_info {
  static const char* name() { return "RNGSSE GM19"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM19 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm19 =
  keyword< rngsse_gm19_info, r,n,g,s,s,e,'_',g,m,'1','9' >;

// Keyword 'rngsse_gm29'
struct rngsse_gm29_info {
  static const char* name() { return "RNGSSE GM29"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM29 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm29 =
  keyword< rngsse_gm29_info, r,n,g,s,s,e,'_',g,m,'2','9' >;

// Keyword 'rngsse_gm31'
struct rngsse_gm31_info {
  static const char* name() { return "RNGSSE GM31"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM31 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm31 =
  keyword< rngsse_gm31_info, r,n,g,s,s,e,'_',g,m,'3','1' >;

// Keyword 'rngsse_gm55'
struct rngsse_gm55_info {
  static const char* name() { return "RNGSSE GM55"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM55 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm55 =
  keyword< rngsse_gm55_info, r,n,g,s,s,e,'_',g,m,'5','5' >;

// Keyword 'rngsse_gm61'
struct rngsse_gm61_info {
  static const char* name() { return "RNGSSE GM61"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM61 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm61 =
  keyword< rngsse_gm61_info, r,n,g,s,s,e,'_',g,m,'6','1' >;

// Keyword 'rngsse_gq58.1'
struct rngsse_gq581_info {
  static const char* name() { return "RNGSSE GQ58.1"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GQ58.1 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gq581 =
  keyword< rngsse_gq581_info, r,n,g,s,s,e,'_',g,q,'5','8','.','1' >;

// Keyword 'rngsse_gq58.3'
struct rngsse_gq583_info {
  static const char* name() { return "RNGSSE GQ58.3"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GQ58.3 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gq583 =
  keyword< rngsse_gq583_info, r,n,g,s,s,e,'_',g,q,'5','8','.','3' >;

// Keyword 'rngsse_gq58.4'
struct rngsse_gq584_info {
  static const char* name() { return "RNGSSE GQ58.4"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GQ58.4 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gq584 =
  keyword< rngsse_gq584_info, r,n,g,s,s,e,'_',g,q,'5','8','.','4' >;

// Keyword 'rngsse_mt19937'
struct rngsse_mt19937_info {
  static const char* name() { return "RNGSSE MT19937"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's MT19937 generator, a "
    "Mersenne Twister generator.";
  }
};

using rngsse_mt19937 =
  keyword< rngsse_mt19937_info, r,n,g,s,s,e,'_',m,t,'1','9','9','3','7' >;

// Keyword 'rngsse_lfsr113'
struct rngsse_lfsr113_info {
  static const char* name() { return "RNGSSE LFSR113"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's LFSR113 generator.";
  }
};

using rngsse_lfsr113 =
  keyword< rngsse_lfsr113_info, r,n,g,s,s,e,'_',l,f,s,r,'1','1','3' >;

// Keyword 'rngsse_mrg32k3a'
struct rngsse_mrg32k3a_info {
  static const char* name() { return "RNGSSE MRG32K3A"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's MRG32K3A generator, a "
    "combined multiple recursive random number generator with two components "
    "of order 3.";
  }
};

using rngsse_mrg32k3a =
  keyword< rngsse_mrg32k3a_info, r,n,g,s,s,e,'_',m,r,g,'3','2',k,'3',a >;

// Keyword 'seqlen'
struct seqlen_info {
  static const char* name() { return "RNGSSE sequence length"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's generator sequence "
    "length. Valid options are 'short', 'medium', and 'long'.";
  }
};

using seqlen =
  keyword< seqlen_info, s,e,q,l,e,n >;

// Keyword 'short'
struct seq_short_info {
  static const char* name() { return "short"; }
  static const char* help() { return
    "This keyword is used to select the short sequence length for RNGSSE "
    "library's generator sequence.";
  }
};

using seq_short =
  keyword< seq_short_info, s,h,o,r,t >;

// Keyword 'medium'
struct seq_med_info {
  static const char* name() { return "medium"; }
  static const char* help() { return
    "This keyword is used to select the medium sequence length for RNGSSE "
    "library's generator sequence.";
  }
};

using seq_med =
  keyword< seq_med_info, m,e,d,i,u,m >;

// Keyword 'long'
struct seq_long_info {
  static const char* name() { return "long"; }
  static const char* help() { return
    "This keyword is used to select the long sequence length for RNGSSE "
    "library's generator sequence.";
  }
};

using seq_long =
  keyword< seq_long_info, l,o,n,g >;

} // kw::
} // tk::
