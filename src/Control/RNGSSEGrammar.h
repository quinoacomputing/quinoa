//******************************************************************************
/*!
  \file      src/Control/RNGSSEGrammar.h
  \author    J. Bakosi
  \date      Wed 14 Jan 2015 02:21:40 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     RNGSSE-related grammar
  \details   This file defines RNGSSE2 library related grammar, (re-)used by
     several executables.
*/
//******************************************************************************
#ifndef RNGSSEGrammar_h
#define RNGSSEGrammar_h

namespace tk {
//! Toolkit, grammar definition for the RNGSSE library
namespace rngsse {

  // RNGSSE PEGTL actions

  //! \brief convert and insert sequence option value to map at position given
  //!    by tags
  //! \author J. Bakosi
  template< class Stack, template < class > class use, class Option,
            typename field, typename sel, typename vec, typename tag,
            typename... tags >
  struct insert_seq : pegtl::action_base<
                        insert_seq< Stack, use, Option, field, sel, vec, tag,
                                    tags... > > {
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
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      boost::mpl::for_each< typename Option::keywords >
                          ( tk::grm::is_keyword< use >() );
    }
  };

  // RNGSSE PEGTL grammar

  //! \brief rng: match any one of the RNGSSE random number generators
  //! \author J. Bakosi
  template< template< class > class use >
  struct rng :
         pegtl::sor< typename use< kw::rngsse_gm19 >::pegtl_string,
                     typename use< kw::rngsse_gm29 >::pegtl_string,
                     typename use< kw::rngsse_gm31 >::pegtl_string,
                     typename use< kw::rngsse_gm55 >::pegtl_string,
                     typename use< kw::rngsse_gm61 >::pegtl_string,
                     typename use< kw::rngsse_gq581 >::pegtl_string,
                     typename use< kw::rngsse_gq583 >::pegtl_string,
                     typename use< kw::rngsse_gq584 >::pegtl_string,
                     typename use< kw::rngsse_mt19937 >::pegtl_string,
                     typename use< kw::rngsse_lfsr113 >::pegtl_string,
                     typename use< kw::rngsse_mrg32k3a >::pegtl_string > {};

  //! \brief Match and set RNGSSE RNG seed
  //! \author J. Bakosi
  template< typename Stack, template< class > class use, typename sel,
            typename vec, typename... tags >
  struct seed :
         tk::grm::process< Stack,
                           use< kw::seed >,
                           tk::grm::Insert_field< Stack,
                                                  tag::seed,
                                                  sel, vec, tags... > > {};

  //! \brief Match and set RNG sequence length parameter
  //! \author J. Bakosi
  template< typename Stack, template < class > class use, typename keyword,
            typename option, typename field, typename sel, typename vec,
            typename... tags >
  struct rngsse_seq :
         grm::process< Stack,
                       keyword,
                       insert_seq< Stack, use, option, field, sel, vec,
                                   tags... >,
                       pegtl::alpha > {};

  //! \brief Match and set RNGSSE sequence length
  //! \author J. Bakosi
  template< typename Stack, template< class > class use, typename sel,
            typename vec, typename... tags >
  struct seqlen :
         rngsse_seq< Stack, use,
                     use< kw::seqlen >,
                     ctr::RNG,
                     tag::seqlen,
                     sel, vec, tags... > {};

  //! \brief Match RNGSSE RNGs in an rngs ... end block
  //! \see walker::deck::rngs
  //! \author J. Bakosi
  template< typename Stack, template< class > class use, typename sel,
            typename vec, typename... tags >
  struct rngs :
         pegtl::ifmust<
           tk::grm::scan< Stack, rng< use >,
                          tk::grm::store_back_option< Stack,
                                                      use,
                                                      ctr::RNG,
                                                      sel, vec > >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           seed< Stack, use, sel, vec, tags... >,
                           seqlen< Stack, use, sel, vec, tags... > > > {};

} // rngsse::
} // tk::

#endif // RNGSSEGrammar_h
