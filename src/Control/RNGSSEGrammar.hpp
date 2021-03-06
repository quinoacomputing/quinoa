// *****************************************************************************
/*!
  \file      src/Control/RNGSSEGrammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     RNGSSE-related grammar
  \details   This file defines RNGSSE2 library related grammar, (re-)used by
     several executables.
*/
// *****************************************************************************
#ifndef RNGSSEGrammar_h
#define RNGSSEGrammar_h

#include <brigand/algorithms/for_each.hpp>

#include "CommonGrammar.hpp"

namespace tk {
namespace grm {

  // Note that PEGTL action specializations must be in the same namespace as the
  // template being specialized. See http://stackoverflow.com/a/3052604.

  // RNGSSE PEGTL actions

  //! Rule used to trigger action
  template< template < class > class use, class Option,
            typename field, typename sel, typename vec, typename tag,
            typename... tags > struct insert_seq : pegtl::success {};
  //! \brief Convert and insert RNGSSE sequence option value to map at position
  //!   given by tags
  template< template < class > class use, class Option,
            typename field, typename sel, typename vec, typename tag,
            typename... tags >
  struct action< insert_seq< use, Option, field, sel, vec, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      ctr::RNGSSESeqLen opt;
      using EnumType = ctr::RNGSSESeqLen::EnumType;
      // get recently inserted key from <sel,vec>
      const auto& key = stack.template get< sel, vec >().back();
      // Error out if RNG does not support option specified
      if ( !ctr::RNG().supportsOpt( key, opt.value(in.string()) ) )
        Message< Stack, ERROR, MsgKey::UNSUPPORTED >( stack, in );
      stack.template
        insert_field< field, EnumType, tag, tags... >
                    ( key, opt.value(in.string()) );
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      brigand::for_each< typename Option::keywords >( is_keyword< use >() );
    }
  };

} // ::grm

//! Toolkit, grammar definition for the RNGSSE library
namespace rngsse {

  using namespace tao;

  // RNGSSE PEGTL grammar

  //! \brief rng: match any one of the RNGSSE random number generators
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
  template< template< class > class use, typename sel,
            typename vec, typename... tags >
  struct seed :
         tk::grm::process< use< kw::seed >,
                           tk::grm::insert_seed< sel, vec, tags... > > {};

  //! \brief Match and set RNG sequence length parameter
  template< template < class > class use, typename keyword,
            typename option, typename field, typename sel, typename vec,
            typename... tags >
  struct rngsse_seq :
         grm::process<
           keyword,
           tk::grm::insert_seq< use, option, field, sel, vec, tags... >,
           pegtl::alpha > {};

  //! \brief Match and set RNGSSE sequence length
  template< template< class > class use, typename sel,
            typename vec, typename... tags >
  struct seqlen :
         rngsse_seq< use,
                     use< kw::seqlen >,
                     ctr::RNG,
                     tag::seqlen,
                     sel, vec, tags... > {};

  //! \brief Match RNGSSE RNGs in an rngs ... end block
  //! \see walker::deck::rngs
  template< template< class > class use, typename sel,
            typename vec, typename... tags >
  struct rngs :
         pegtl::if_must<
           tk::grm::scan< rng< use >,
                          tk::grm::store_back_option< use,
                                                      ctr::RNG,
                                                      sel, vec > >,
           tk::grm::block< use< kw::end >,
                           seed< use, sel, vec, tags... >,
                           seqlen< use, sel, vec, tags... > > > {};

} // rngsse::
} // tk::

#endif // RNGSSEGrammar_h
