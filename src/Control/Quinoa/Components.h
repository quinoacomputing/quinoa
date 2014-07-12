//******************************************************************************
/*!
  \file      src/Control/Quinoa/Components.h
  \author    J. Bakosi
  \date      Sat 22 Feb 2014 02:26:00 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
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
//! components are changed.
template< typename... Tags >
class ncomponents : public tk::tuple::tagged_tuple< Tags... > {

  private:
    // Create mpl::list from variadic pack, storing no. of comps with types
    using tagsncomps = typename tk::make_list< Tags... >::type;

  public:
    // Find out 'number of components' type
    using ncomp =
      typename boost::mpl::at< tagsncomps, boost::mpl::int_<1> >::type;

    // Remove ncomp types (i.e., keep only the no. of comps, i.e., the tags)
    using tags = typename boost::mpl::remove< tagsncomps, ncomp >::type;

  private:
    //! Function object for zeroing all number of components
    struct zero {
      ncomponents* const m_host;
      zero( ncomponents* const host ) : m_host( host ) {}
      template< typename U > void operator()( U ) {
        m_host->template set< U >( 0 );
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
        m_nprop += m_host->template get< U >();
      }
    };

    //! Function object for computing the offset for a given tag (i.e., the sum
    //! of the number of components up to a given tag)
    template< typename tag >
    struct addncomp4tag {
      const ncomponents* const m_host;
      ncomp& m_offset;
      bool m_found;
      addncomp4tag( const ncomponents* const host, ncomp& offset ) :
        m_host( host ), m_offset( offset = 0 ), m_found( false ) {}
      template< typename U > void operator()( U ) {
        if (std::is_same< tag, U >::value) {
          m_found = true;
        } else if (!m_found) {
          m_offset += m_host->template get< U >();
        }
      }
    };

  public:
    //! Default constructor: set defaults as zero for all number of components
    ncomponents() {
      boost::mpl::for_each< tags >( zero( this ) );
    }

    //! Return total number of components
    ncomp nprop() const noexcept {
      ncomp n;
      boost::mpl::for_each< tags >( addncomp( this, n ) );
      return n;
    }

    //! Return offset for tag
    template< typename tag >
    ncomp offset() const noexcept {
      ncomp n;
      boost::mpl::for_each< tags >( addncomp4tag< tag >( this, n ) );
      return n;
    }
};

//! Number of components, storing using type given
template< typename ncomp >
using ncomps = ncomponents<
  tag::position,   ncomp,      //!< Position model
  tag::mass,       ncomp,      //!< Mass model
  tag::hydro,      ncomp,      //!< Hydro model
  tag::mix,        ncomp,      //!< Material mix model
  tag::frequency,  ncomp,      //!< Turbulent frequency model

  tag::dirichlet,  ncomp,      //!< Dirichlet SDE
  tag::gendir,     ncomp,      //!< Generalized Dirichlet SDE
  tag::ou,         ncomp,      //!< Ornstein-Uhlenbeck SDE
  tag::lognormal,  ncomp,      //!< Log-normal SDE
  tag::skewnormal, ncomp,      //!< Skew-normal SDE
  tag::gamma,      ncomp,      //!< Gamma SDE
  tag::beta,       ncomp       //!< Beta SDE
>;

} // ctr::
} // quinoa::

#endif // QuinoaComponents_h
