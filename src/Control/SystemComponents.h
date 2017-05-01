// *****************************************************************************
/*!
  \file      src/Control/SystemComponents.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Operations on numbers of scalar components of systems of equations
  \details   Operations on numbers of scalar components of systems of equations,
    e.g. multiple equation sets of a physical model or a set of (stochastic
    ordinary or partial differential) equations of different kinds.

    _Problem:_ We are given a type that is a tk::tuple::tagged_tuple that
    contains an arbitrary number of std::vectors of integers. The number of
    vectors are fixed at compile-time (accessed via tags) but their (differing)
    length is only known at run-time (after parsing user input). What we need
    are functions that operate on this data structure and return, e.g., the
    total number of components in the whole system or the offset in the whole
    data structure for a given tag. The functions should thus be able to operate
    on a list of types, i.e., a double for loop over all tags and associated
    vectors - one at compile-time and the other one at run-time.

    _Example:_ An example is to define storage for systems of stochastic
    differential equations as in walker::ctr::ncomps. In
    [walker](walker.html) various types of stochastic differential
    equations are implemented and numerically integrated in time in order to
    compute and learn about their statistical behavior. Since different types of
    systems can be parameterized differently, there can be several systems of
    the same type (the number is user-defined so only known at run-time). For
    example, it is possible to numerically integrate two systems of Dirichlet
    equations, see DiffEq/Dirichlet.h, one with 4, the other one with 5 scalar
    variables, parameterized differently, and estimate their arbitrary coupled
    statistics and joint PDFs. This requires that each type of system
    has a vector of integers storing the number of scalar variables.

    _Solution:_ Looping through elements of a tuple is done via [Boost]
    (http://www.boost.org)'s [MetaProgramming Library (MPL)]
    (http://www.boost.org/doc/libs/release/libs/mpl)'s [for_each]
    (http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/for-each.html).
    Such operations on types happen at compile-time, i.e., the code runs inside
    the compiler and only its result gets compiled into code to be run at
    run-time. Advantages are abstraction and generic code that is independent
    of the size and order of the tags in the tuple (and the associated vectors).
    Though it is not of primary concern for this data structure, this solution
    also generates efficient code, since part of the algorithm (for computing
    the offset, the total number of components, etc.) runs inside the compiler.
    All of this is kept generic, so it does not need to be changed when a new
    pair of tag and system of equations are added to the underlying tuple.
    Thus, _the main advantage is that changing the order of the vectors in the
    tuple or adding new ones at arbitrary locations will not require change to
    client-code._

    _An alternative approach:_ Of course this could have been simply done with a
    purely run-time vector of vector of integers. However, that would not have
    allowed a tag-based access to the different systems of equations. (Sure, one
    can define an enum for equation indexing, but that just adds code to what is
    already need to be changed if a change happens to the underlying data
    structure.) As a consequence, a component in a vector of vector would have
    had to be accessed as something like, ncomps[2][3], which is more
    error-prone to read than ncomps< tag::dirichlet >( 3 ). Furthermore, that
    design would have resulted in double-loops at run-time, instead of letting
    the compiler loop over the types at compile-time and do the run-time loops
    only for what is only known at run-time. Additionally, and most importantly,
    _any time a new system is introduced (or the order is changed), all
    client-code would have to be changed._

    For example client-code, see walker::Dirichlet::Dirichlet() in
    DiffEq/Dirichlet.h, or walker::Integrator::Integrator in Walker/Integrator.h.
*/
// *****************************************************************************
#ifndef SystemComponents_h
#define SystemComponents_h

#include "NoWarning/for_each.h"
#include <boost/mpl/remove.hpp>
#include <boost/mpl/at.hpp>

#include "Make_list.h"
#include "TaggedTuple.h"
#include "StatCtr.h"
#include "Keywords.h"

namespace tk {
//! Toolkit control, general purpose user input to internal data transfer
namespace ctr {

//! Inherit type of number of components from keyword 'ncomp'
using ncomp_type = kw::ncomp::info::expect::type;

//! \brief Map associating a dependent variable of a differential equation to
//!   its offset in its data array storing its dependent variables
//! \details We use a case-insensitive character comparison functor for the
//!   offset map, since the keys (dependent variables) in the offset map are
//!   only used to indicate the equation's dependent variable, however, queries
//!   (using std::map's member function, find) can be fired up for both ordinary
//!   and central moments (which are denoted by upper and lower case,
//!   characters, respectively) for which the offset (for the same dependent
//!   variable) should be the same.
using OffsetMap = std::map< char, ncomp_type, CaseInsensitiveCharLess >;

//! \brief Number of components storage
//! \details All this trickery with boost::mpl allows the code below to be
//!   generic. As a result, adding a new component requires adding a single line
//!   (a tag and its type) to the already existing list, see typedefs 'ncomps'.
//!   The member functions, doing initialization, computing the number of total
//!   components, the offset for a given tag, and computing the offset map, need
//!   no change -- even if the order of the number of components change.
template< typename... Tags >
class ncomponents : public tk::tuple::tagged_tuple< Tags... > {

  private:
    //! Create mpl::list from variadic pack, storing no. of comps with types
    using tagsncomps = typename tk::make_list< Tags... >::type;

  public:
    //! Remove std::vector< ncomp_type > types, i.e., keep only the tags
    using tags = typename
      boost::mpl::remove< tagsncomps, std::vector< ncomp_type > >::type;

  private:
    //! Function object for zeroing all number of components
    //! \author J. Bakosi
    struct zero {
      //! Need to store reference to host class whose data we operate on
      ncomponents* const m_host;
      //! Constructor: store host object pointer
      zero( ncomponents* const host ) : m_host( host ) {}
      //! Function call operator templated on the type that does the zeroing
      template< typename U > void operator()( U ) {
        //! Loop through and zero all elements of the vector for this system
        for (auto& c : m_host->template get< U >()) c = 0;
      }
    };

    //! \brief Function object for computing the total number of components
    //!   (i.e., the sum of all of the number of components)
    //! \author J. Bakosi
    struct addncomp {
      //! Need to store reference to host class whose data we operate on
      const ncomponents* const m_host;
      //! Internal reference used to return the total number of components
      ncomp_type& m_nprop;
      //! \brief Constructor: store host object pointer and initially zeroed
      //!   counter reference
      addncomp( const ncomponents* const host, ncomp_type& nprop ) :
        m_host( host ), m_nprop( nprop = 0 ) {}
      //! Function call operator templated on the type that does the counting
      template< typename U > void operator()( U ) {
        //! Loop through and add up all elements of the vector for this system
        for (const auto& c : m_host->template get< U >()) m_nprop += c;
      }
    };

    //! \brief Function object for computing the offset for a given tag (i.e.,
    //!   the sum of the number of components up to a given tag)
    //! \details This is used to index into the data array (the equation systems
    //!   operate on during the numerical solution) and get to the beginning of
    //!   data for a given differential equation system.
    //! \author J. Bakosi
    template< typename tag >
    struct addncomp4tag {
      //! Need to store reference to host class whose data we operate on
      const ncomponents* const m_host;
      //! Internal reference used to return the offset for the tag given
      ncomp_type& m_offset;
      //! \brief Internal storage for the index of a system within systems
      //! \details Example: I want the second Dirichlet system: m_c = 1.
      //! \see offset().
      const ncomp_type m_c;
      //! Indicates whether the tag (eq system) was found, so it is time to quit
      bool m_found;
      //! \brief Constructor: store host object pointer, initially zeroed offset
      //!   reference, and system index we are looking for
      addncomp4tag( const ncomponents* const host, ncomp_type& offset,
                    ncomp_type c ) :
        m_host( host ), m_offset( offset = 0 ), m_c( c ), m_found( false ) {}
      //! \brief Function call operator templated on the type that does the
      //!   offset calculation
      template< typename U > void operator()( U ) {
        if (std::is_same< tag, U >::value) {
          // Make sure we are not trying to index beyond the length for this U
          Assert( m_host->template get<U>().size() >= m_c,
                  "Indexing out of bounds in addncomp4tag!" );
          // If we have found the tag we are looking for, we count up to the
          // given system index (passed in via the constructor) and add those to
          // the offset, then quit
          for (ncomp_type c=0; c<m_c; ++c)
            m_offset += m_host->template get<U>()[c];
          m_found = true;
        } else if (!m_found) {
          // If we have not found the tag we are looking for, we add all the
          // number of scalars for that tag to the offset
          for (const auto& c : m_host->template get< U >()) m_offset += c;
        }
      }
    };

  public:
    //! Default constructor: set defaults to zero for all number of components
    //! \author J. Bakosi
    ncomponents() { boost::mpl::for_each< tags >( zero( this ) ); }

    //! \return Total number of components
    //! \author J. Bakosi
    ncomp_type nprop() const noexcept {
      ncomp_type n;
      boost::mpl::for_each< tags >( addncomp( this, n ) );
      return n;
    }

    //! \return offset for tag
    //! \param[in] c Index for system given by template argument tag
    //! \author J. Bakosi
    template< typename tag >
    ncomp_type offset( ncomp_type c ) const noexcept {
      ncomp_type n;
      boost::mpl::for_each< tags >( addncomp4tag< tag >( this, n, c ) );
      return n;
    }

    //! \return Map of offsets associated to dependent variables
    //! \param[in] depvars Vector of vectors of dependent variables
    //! \author J. Bakosi
    OffsetMap offsetmap( const std::vector< std::vector< char > >& depvars )
    const {
      OffsetMap o;
      // This is now simple, as the heavy lifting has already been done by
      // extracting the dependent variables to a vector of vectors that is
      // passed in. Note that for now, we simply assign linearly increasing
      // offsets. This may be changed in the future if/when we decide to
      // implement a user-configurable order of the data layout for the various
      // equations.
      ncomp_type c = 0;
      for (const auto& vec : depvars) for (const auto& v : vec) o[ v ] = c++;
      return o;
    }
};

} // ctr::
} // tk::

#endif // SystemComponents_h
