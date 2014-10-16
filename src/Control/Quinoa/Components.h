//******************************************************************************
/*!
  \file      src/Control/Quinoa/Components.h
  \author    J. Bakosi
  \date      Wed 20 Aug 2014 09:08:38 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Storage for number of components
  \details   Storage for number of components
*/
//******************************************************************************
#ifndef QuinoaComponents_h
#define QuinoaComponents_h

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/remove.hpp>
#include <boost/mpl/at.hpp>

#include <make_list.h>

namespace quinoa {
namespace ctr {

//! Number of components storage: All this trickery with boost::mpl allows the
//! code below to be generic. As a result, adding a new component requires
//! adding a single line (a tag and its type) to the already existing list, see
//! the templated typedef 'ncomps' below. The member functions, doing
//! initialization, computing the number of total components, and the offset for
//! a given tag needs no changes -- even if the order of the number of
//! components change.
template< typename... Tags >
class ncomponents : public tk::tuple::tagged_tuple< Tags... > {

  private:
    // Create mpl::list from variadic pack, storing no. of comps with types
    using tagsncomps = typename tk::make_list< Tags... >::type;

  public:
    // Find out 'number of components' type
    using ncomp = typename boost::mpl::at< tagsncomps, boost::mpl::int_<1> >
                  ::type::value_type;

    // Remove std::vector< ncomp > types (i.e., keep only the no. of comps,
    // i.e., the tags)
    using tags = typename boost::mpl::remove< tagsncomps, std::vector< ncomp > >
                 ::type;

  private:
    //! Function object for zeroing all number of components
    struct zero {
      ncomponents* const m_host;
      zero( ncomponents* const host ) : m_host( host ) {}
      template< typename U > void operator()( U ) {
        for (auto& c : m_host->template get< U >()) c = 0;
      }
    };

    //! Function object for computing the total number of components (i.e., the
    //! sum of all of the number of components)
    struct addncomp {
      const ncomponents* const m_host;
      ncomp& m_nprop;
      addncomp( const ncomponents* const host, ncomp& nprop ) :
        m_host( host ), m_nprop( nprop = 0 ) {}
      template< typename U > void operator()( U ) {
        for (const auto& c : m_host->template get< U >()) m_nprop += c;
      }
    };

    //! Function object for computing the offset for a given tag (i.e., the sum
    //! of the number of components up to a given tag)
    template< typename tag >
    struct addncomp4tag {
      const ncomponents* const m_host;
      ncomp& m_offset;
      const ncomp m_c;
      bool m_found;
      addncomp4tag( const ncomponents* const host, ncomp& offset, ncomp c ) :
        m_host( host ), m_offset( offset = 0 ), m_c( c ), m_found( false ) {}
      template< typename U > void operator()( U ) {
        if (std::is_same< tag, U >::value) {
          for (ncomp c=0; c<m_c; ++c) m_offset += m_host->template get<U>()[c];
          m_found = true;
        } else if (!m_found) {
          for (const auto& c : m_host->template get< U >()) m_offset += c;
        }
      }
    };

  public:
    //! Default constructor: set defaults as zero for all number of components
    ncomponents() { boost::mpl::for_each< tags >( zero( this ) ); }

    //! Return total number of components
    ncomp nprop() const noexcept {
      ncomp n;
      boost::mpl::for_each< tags >( addncomp( this, n ) );
      return n;
    }

    //! Return offset for tag
    template< typename tag >
    ncomp offset( ncomp c ) const noexcept {
      ncomp n;
      boost::mpl::for_each< tags >( addncomp4tag< tag >( this, n, c ) );
      return n;
    }
};

//! Number of components of models and equations
template< typename ncomp >
using ncomps = ncomponents<
  tag::position,   std::vector< ncomp >,      //!< Position models
  tag::mass,       std::vector< ncomp >,      //!< Mass models
  tag::hydro,      std::vector< ncomp >,      //!< Hydro models
  tag::mix,        std::vector< ncomp >,      //!< Material mix models
  tag::frequency,  std::vector< ncomp >,      //!< Turbulent frequency models

  tag::dirichlet,  std::vector< ncomp >,      //!< Dirichlet SDEs
  tag::gendir,     std::vector< ncomp >,      //!< Generalized Dirichlet SDEs
  tag::ou,         std::vector< ncomp >,      //!< Ornstein-Uhlenbeck SDEs
  tag::lognormal,  std::vector< ncomp >,      //!< Log-normal SDEs
  tag::skewnormal, std::vector< ncomp >,      //!< Skew-normal SDEs
  tag::gamma,      std::vector< ncomp >,      //!< Gamma SDEs
  tag::beta,       std::vector< ncomp >       //!< Beta SDEs
>;

} // ctr::
} // quinoa::

#endif // QuinoaComponents_h
