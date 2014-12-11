//******************************************************************************
/*!
  \file      src/Control/Components.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 09:34:17 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Storage for number of components
  \details   Storage for number of components
*/
//******************************************************************************
#ifndef Components_h
#define Components_h

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/remove.hpp>
#include <boost/mpl/at.hpp>

#include <make_list.h>

#include <ControlTypes.h>

namespace tk {
//! Toolkit control, general purpose user input to internal data transfer
namespace ctr {

// Map associating a dependend variable to its offset in the particle data.
// We use a case-insensitive character comparison functor for the offset
// map, since the keys (dependent variables) in the offset map are only used
// to indicate the equation's dependent variable, however, queries (using
// find) are fired up for both ordinary and central moments (which are upper
// and lower case) for which the offset (for the same depvar) should be the
// same.
using OffsetMap = std::map< char, unsigned int, CaseInsensitiveCharLess >;

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
    // Remove std::vector< unsigned int > types, i.e., keep only the tags
    using tags = typename
      boost::mpl::remove< tagsncomps, std::vector< unsigned int > >::type;

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
      unsigned int& m_nprop;
      addncomp( const ncomponents* const host, unsigned int& nprop ) :
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
      unsigned int& m_offset;
      const unsigned int m_c;
      bool m_found;
      addncomp4tag( const ncomponents* const host, unsigned int& offset,
                    unsigned int c ) :
        m_host( host ), m_offset( offset = 0 ), m_c( c ), m_found( false ) {}
      template< typename U > void operator()( U ) {
        if (std::is_same< tag, U >::value) {
          for (unsigned int c=0; c<m_c; ++c)
            m_offset += m_host->template get<U>()[c];
          m_found = true;
        } else if (!m_found) {
          for (const auto& c : m_host->template get< U >()) m_offset += c;
        }
      }
    };

    //! Function object for creating a run-time std::map< depvar, offset >
    struct map {
      const ncomponents* const m_host;
      OffsetMap& m_map;
      const std::vector< std::vector< char > >& m_depvars;
      map( const ncomponents* const host,
           const std::vector< std::vector< char > >& depvars,
           OffsetMap& map ) :
        m_host( host ), m_map( map ), m_depvars( depvars ) {}
      template< typename U > void operator()( U ) {
        // Loop through the vector of dependent variables for a component type
        // given by type U and insert an entry to the offset map for each, i.e.,
        // handle all potentially multiple components of the same type.
        unsigned int c = 0;    // offset counter
        for (const auto& vec : m_depvars)
          for (const auto& v : vec)
            m_map[ v ] = m_host->template offset< U >( c++ );
      }
    };

  public:
    //! Default constructor: set defaults as zero for all number of components
    ncomponents() { boost::mpl::for_each< tags >( zero( this ) ); }

    //! Return total number of components
    unsigned int nprop() const noexcept {
      unsigned int n;
      boost::mpl::for_each< tags >( addncomp( this, n ) );
      return n;
    }

    //! Return offset for tag
    template< typename tag >
    unsigned int offset( unsigned int c ) const noexcept {
      unsigned int n;
      boost::mpl::for_each< tags >( addncomp4tag< tag >( this, n, c ) );
      return n;
    }

    //! Return map of offsets associated to dependent variables
    OffsetMap offsetmap( const std::vector< std::vector< char > >& depvars )
    const {
      OffsetMap offset;
      boost::mpl::for_each< tags >( map( this, depvars, offset ) );
      return offset;
    }
};

} // ctr::
} // tk::

#endif // Components_h
