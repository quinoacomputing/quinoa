//******************************************************************************
/*!
  \file      src/Control/Quinoa/Components.h
  \author    J. Bakosi
  \date      Thu 13 Feb 2014 07:36:04 PM CET
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
//! the templated typedef 'comps' below. The member functions, doing
//! initialization, computing the number of total components, and the offset for
//! a given tag needs no changes -- even if the order of the number of
//! components are changed.
template< typename... tags >
class components : public tk::tuple::tagged_tuple< tags... > {

  private:
    // Create mpl::list from variadic pack, storing no. of comps with types
    using ncompsWtypes = typename tk::make_list< tags... >::type;
    // Find out 'number of components' type
    using ncomp =
      typename boost::mpl::at< ncompsWtypes, boost::mpl::int_<1> >::type;
    // Remove ncomp types (i.e., keep only the no. of comps, i.e., the tags)
    using vectags = typename boost::mpl::remove< ncompsWtypes, ncomp >::type;

    //! Function object for zeroing all number of components
    template< class Host >
    struct zero {
      Host* const m_host;
      zero( Host* const host ) : m_host( host ) {}
      template< typename U > void operator()( U ) {
        m_host->template set< U >( 0 );
      }
    };

    //! Function object for computing the total number of components (i.e., the
    //! sum of all of the number of components)
    template< class Host >
    struct addncomp {
      const Host* const m_host;
      ncomp& m_nprop;
      addncomp( const Host* const host, ncomp& nprop ) :
        m_host( host ), m_nprop( nprop = 0 ) {}
      template< typename U > void operator()( U ) {
        m_nprop += m_host->template get< U >();
      }
    };

    //! Function object for computing the offset for a given tag (i.e., the sum
    //! of the number of components up to a given tag)
    template< class Host, typename tag >
    struct addncomp4tag {
      const Host* const m_host;
      ncomp& m_offset;
      bool m_found;
      addncomp4tag( const Host* const host, ncomp& offset ) :
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
    components() {
      boost::mpl::for_each< vectags >( zero< components >( this ) );
    }

    //! Return total number of components
    ncomp nprop() const noexcept {
      ncomp n;
      boost::mpl::for_each< vectags >( addncomp< components >( this, n ) );
      return n;
    }

    //! Return offset for tag
    template< typename tag >
    ncomp offset() const noexcept {
      ncomp n;
      boost::mpl::for_each< vectags >( addncomp4tag<components,tag>(this,n) );
      return n;
    }
};

template< typename ncomp >
using comps = components<
  tag::nposition,   ncomp, //!< N. of position components in position model
  tag::ndensity,    ncomp, //!< N. of density components in mass model
  tag::nvelocity,   ncomp, //!< N. of velocity components in hydro model
  tag::nscalar,     ncomp, //!< N. of mixing scalars in material mix model
  tag::nfrequency,  ncomp, //!< N. of frequencies in turb. frequency model
  tag::ndirichlet,  ncomp, //!< N. of components in Dirichlet SDE
  tag::ngendir,     ncomp, //!< N. of components in gen. Dirichlet SDE
  tag::nou,         ncomp, //!< N. of components in Ornstein-Uhlenbeck SDE
  tag::nlognormal,  ncomp, //!< N. of components in Log-normal SDE
  tag::nskewnormal, ncomp, //!< N. of components in Skew-normal SDE
  tag::ngamma,      ncomp, //!< N. of components in gamma SDE
  tag::nbeta,       ncomp  //!< N. of components in beta SDE
>;

} // ctr::
} // quinoa::

#endif // QuinoaComponents_h
