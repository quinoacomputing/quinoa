//******************************************************************************
/*!
  \file      src/Control/MKLGrammar.h
  \author    J. Bakosi
  \date      Tue 13 Jan 2015 12:35:04 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Intel MKL-related grammar
  \details   This file defines Intel Math Kernel Library related grammar,
    (re-)used by several executables.
*/
//******************************************************************************
#ifndef MKLGrammar_h
#define MKLGrammar_h

namespace tk {
//! Toolkit, grammar definition for Intel's Math Kernel Library
namespace mkl {

  //! \brief rng: match any one of the MKL random number generators
  //! \author J. Bakosi
  template< template< class > class use >
  struct rng :
         pegtl::sor< typename use< kw::mkl_mcg31 >::pegtl_string,
                     typename use< kw::mkl_r250 >::pegtl_string,
                     typename use< kw::mkl_mrg32k3a >::pegtl_string,
                     typename use< kw::mkl_mcg59 >::pegtl_string,
                     typename use< kw::mkl_wh >::pegtl_string,
                     typename use< kw::mkl_mt19937 >::pegtl_string,
                     typename use< kw::mkl_mt2203 >::pegtl_string,
                     typename use< kw::mkl_sfmt19937 >::pegtl_string,
                     typename use< kw::mkl_sobol >::pegtl_string,
                     typename use< kw::mkl_niederr >::pegtl_string,
                     //typename use< kw::mkl_iabstract >::pegtl_string,
                     //typename use< kw::mkl_dabstract >::pegtl_string,
                     //typename use< kw::mkl_sabstract >::pegtl_string,
                     typename use< kw::mkl_nondeterm >::pegtl_string > {};

  //! \brief Match and set MKL RNG seed
  //! \author J. Bakosi
  template< typename Stack, template< class > class use, typename sel,
            typename vec, typename... tags >
  struct seed :
         tk::grm::process< Stack,
                           use< kw::seed >,
                           tk::grm::Insert_field< Stack,
                                                  tag::seed,
                                                  sel, vec, tags... > > {};

  //! \brief Match and set MKL uniform method algorithm
  //! \author J. Bakosi
  template< typename Stack, template< class > class use, typename sel,
            typename vec, typename... tags >
  struct uniform_method :
         grm::rng_option< Stack, use,
                          use< kw::uniform_method >,
                          ctr::MKLUniformMethod,
                          tag::uniform_method,
                          sel, vec, tags... > {};
     
  //! \brief Match and set MKL Gaussian method algorithm
  //! \author J. Bakosi
  template< typename Stack, template< class > class use, typename sel,
            typename vec, typename... tags >
  struct gaussian_method :
         grm::rng_option< Stack, use,
                          use< kw::gaussian_method >,
                          ctr::MKLGaussianMethod,
                          tag::gaussian_method,
                          sel, vec, tags... > {};

  //! \brief Match MKL RNGs in an rngs ... end block
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
                           kw::end,
                           seed< Stack, use, sel, vec, tags... >,
                           uniform_method< Stack, use, sel, vec, tags... >,
                           gaussian_method< Stack, use, sel, vec, tags... > > >
  {};

} // mkl::
} // tk::

#endif // MKLGrammar_h
