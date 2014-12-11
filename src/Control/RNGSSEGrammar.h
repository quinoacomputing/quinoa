//******************************************************************************
/*!
  \file      src/Control/RNGSSEGrammar.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 02:21:00 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     RNGSSE grammar
  \details   RNGSSE grammar
*/
//******************************************************************************
#ifndef RNGSSEGrammar_h
#define RNGSSEGrammar_h

namespace tk {
//! Toolkit, grammar definition for the RNGSSE library
namespace rngsse {

  // RNGSSE PEGTL actions

  //! convert and insert sequence option value to map at position given by tags
  template< class Stack, class Option, typename field, typename sel,
            typename vec, typename tag, typename... tags >
  struct insert_seq : pegtl::action_base< insert_seq< Stack, Option, field, sel,
                                                      vec, tag, tags... > > {
    static void apply( const std::string& value, Stack& stack ) {
      ctr::RNGSSESeqLen opt;
      using EnumType = ctr::RNGSSESeqLen::EnumType;
      // get recently inserted key from <sel,vec>
      using key_type =
        typename Stack::template nT< sel >::template nT< vec >::value_type;
      const key_type& key = stack.template get< sel, vec >().back();
      // Error out if RNG does not support option specified
      if ( !ctr::RNG().supportsOpt( key, opt.value(value) ) ) {
        grm::Message< Stack, grm::ERROR, grm::MsgKey::UNSUPPORTED >
                    ( stack, value );
      }
      stack.template insert_opt< key_type, field, EnumType, tag, tags... >
                               ( key, opt.value(value) );
    }
  };

  // RNGSSE PEGTL grammar

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
                           kw::seed::pegtl_string,
                           tk::grm::Insert_field< Stack,
                                                  tag::seed,
                                                  sel, vec, tags... > > {};

  //! RNG sequence length parameter
  template< typename Stack, typename keyword, typename option, typename field,
            typename sel, typename vec, typename... tags >
  struct rngsse_seq :
         grm::process< Stack,
                       typename keyword::pegtl_string,
                       insert_seq< Stack, option, field, sel, vec, tags... >,
                       pegtl::alpha > {};

  //! RNGSSE sequence length
  template< typename Stack, typename sel, typename vec, typename... tags >
  struct seqlen :
         rngsse_seq< Stack,
                     kw::seqlen,
                     ctr::RNG,
                     tag::seqlen,
                     sel, vec, tags... > {};

  //! rngs blocks
  template< typename Stack, typename sel, typename vec, typename... tags >
  struct rngs :
         pegtl::ifmust<
           tk::grm::scan< rng,
                          tk::grm::store_back_option< Stack,
                                                      ctr::RNG,
                                                      sel, vec > >,
           tk::grm::block< Stack,
                           kw::end,
                           seed< Stack, sel, vec, tags... >,
                           seqlen< Stack, sel, vec, tags... > > > {};

} // rngsse::
} // tk::

#endif // RNGSSEGrammar_h
