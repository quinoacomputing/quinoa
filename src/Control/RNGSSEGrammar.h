//******************************************************************************
/*!
  \file      src/Control/RNGSSEGrammar.h
  \author    J. Bakosi
  \date      Fri 27 Dec 2013 07:53:11 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGSSE grammar
  \details   RNGSSE grammar
*/
//******************************************************************************
#ifndef RNGSSEGrammar_h
#define RNGSSEGrammar_h

namespace tk {
//! RNGSSE grammar definition: state, actions, grammar
namespace rngsse {

  //! rng: match any one of the RNGSSE random number generators
  struct rng :
         pegtl::sor< kw::rngsse_gm19::pegtl_string,
                     kw::rngsse_gm29::pegtl_string,
                     kw::rngsse_gm31::pegtl_string,
                     kw::rngsse_gm55::pegtl_string,
                     kw::rngsse_gm61::pegtl_string,
                     kw::rngsse_gq581::pegtl_string,
                     kw::rngsse_gq583::pegtl_string,
                     kw::rngsse_gq584::pegtl_string,
                     kw::rngsse_mt19937::pegtl_string,
                     kw::rngsse_lfsr113::pegtl_string,
                     kw::rngsse_mrg32k3a::pegtl_string > {};

  //! RNGSSE seed
  template< typename Stack, typename sel, typename vec, typename... tags >
  struct seed :
         tk::grm::process< Stack,
                           tk::kw::seed::pegtl_string,
                           tk::grm::Insert_field< Stack,
                                                  quinoa::ctr::seed,
                                                  sel, vec, tags... > > {};
  //! RNGSSE sequence length
  template< typename Stack, typename sel, typename vec, typename... tags >
  struct seqlen :
         grm::rng_option< Stack,
                          tk::kw::seqlen,
                          quinoa::ctr::RNGSSESeqLen,
                          quinoa::ctr::seqlen,
                          sel, vec, tags... > {};

  //! rngs blocks
  template< typename Stack, typename sel, typename vec, typename... tags >
  struct rngs :
         pegtl::ifmust<
           tk::grm::scan< rng,
                          tk::grm::store_back_option< Stack,
                                                      quinoa::ctr::RNG,
                                                      sel, vec > >,
           tk::grm::block< Stack,
                           seed< Stack, sel, vec, tags... >,
                           seqlen< Stack, sel, vec, tags... > > > {};

} // rngsse::
} // tk::

#endif // RNGSSEGrammar_h
