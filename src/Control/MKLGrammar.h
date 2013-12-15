//******************************************************************************
/*!
  \file      src/Control/MKLGrammar.h
  \author    J. Bakosi
  \date      Sun 15 Dec 2013 04:11:07 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL grammar
  \details   MKL grammar
*/
//******************************************************************************
#ifndef MKLGrammar_h
#define MKLGrammar_h

namespace tk {
//! MKL grammar definition: state, actions, grammar
namespace mkl {

  //! rng: match any one of the MKL random number generators
  struct rng :
         pegtl::sor< kw::mkl_mcg31::pegtl_string,
                     kw::mkl_r250::pegtl_string,
                     kw::mkl_mrg32k3a::pegtl_string,
                     kw::mkl_mcg59::pegtl_string,
                     kw::mkl_wh::pegtl_string,
                     kw::mkl_mt19937::pegtl_string,
                     kw::mkl_mt2203::pegtl_string,
                     kw::mkl_sfmt19937::pegtl_string,
                     kw::mkl_sobol::pegtl_string,
                     kw::mkl_niederr::pegtl_string,
                     kw::mkl_iabstract::pegtl_string,
                     kw::mkl_dabstract::pegtl_string,
                     kw::mkl_sabstract::pegtl_string,
                     kw::mkl_nondeterm::pegtl_string > {};

  //! insert MKL parameter
  template< typename Stack, typename keyword, typename option, typename field,
            typename sel, typename vec, typename... tags >
  struct rng_option :
         tk::grm::process< Stack,
                           typename keyword::pegtl_string,
                           tk::grm::insert_option< Stack,
                                                   option,
                                                   field,
                                                   sel, vec, tags... >,
                           pegtl::alpha > {};

  //! MKL RNG seed
  template< typename Stack, typename sel, typename vec, typename... tags >
  struct seed :
         tk::grm::process< Stack,
                           tk::kw::seed::pegtl_string,
                           tk::grm::Insert_field< Stack,
                                                  quinoa::ctr::seed,
                                                  sel, vec, tags... > > {};

  //! MKL uniform method
  template< typename Stack, typename sel, typename vec, typename... tags >
  struct uniform_method :
         rng_option< Stack,
                     tk::kw::uniform_method,
                     quinoa::ctr::MKLUniformMethod,
                     quinoa::ctr::uniform_method,
                     sel, vec, tags... > {};

  //! MKL Gaussian method
  template< typename Stack, typename sel, typename vec, typename... tags >
  struct gaussian_method :
         rng_option< Stack,
                     tk::kw::gaussian_method,
                     quinoa::ctr::MKLGaussianMethod,
                     quinoa::ctr::gaussian_method,
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
                           uniform_method< Stack, sel, vec, tags... >,
                           gaussian_method< Stack, sel, vec, tags... > > > {};

} // mkl::
} // tk::

#endif // MKLGrammar_h
