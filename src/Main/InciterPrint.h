// *****************************************************************************
/*!
  \file      src/Main/InciterPrint.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
// *****************************************************************************
#ifndef InciterPrint_h
#define InciterPrint_h

#include <iostream>
#include <string>

#include "NoWarning/format.h"
#include "NoWarning/for_each.h"

#include "Print.h"
#include "ContainerUtil.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/Options/Physics.h"
#include "Inciter/Options/Problem.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! InciterPrint : tk::Print
class InciterPrint : public tk::Print {

  public:
    //! Constructor
    //! \param[in,out] str Verbose stream
    //! \param[in,out] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    //! \author J. Bakosi
    explicit InciterPrint( std::ostream& str = std::clog,
                           std::ostream& qstr = std::cout ) :
      Print( str, qstr ) {}

    //! Print control option: 'group : option'
    //! \author J. Bakosi
    template< typename Option, typename... tags >
    void Item() const {
      Option opt;
      m_stream << m_item_name_value_fmt
                  % m_item_indent % opt.group()
                  % opt.name( g_inputdeck.get< tags... >() );
    }

    // Helper class for compact output of PDE policies
    //! \author J. Bakosi
    class Policies {
      public:
        // Default constructor
        explicit Policies() : phys(), prob() {}
        // Initializer constructor
        explicit Policies( const std::string& h, const std::string& r ) :
          phys(h), prob(r) {}
        // Operator += for adding up two Policies structs
        Policies& operator+= ( const Policies& p ) {
          phys += p.phys;
          prob += p.prob;
          return *this;
        }
        // Output unique policies to output stream
        friend std::ostream& operator<< ( std::ostream& os, const Policies& p )
        {
          Policies copy( p );     // copy policies
          copy.unique();          // get rid of duplicate policies
          os << "h:" << copy.phys << ", r:" << copy.prob;
          return os;
        }

      private:
        // Make all policies unique
        void unique() { tk::unique( phys ); tk::unique( prob ); }

        std::string phys;
        std::string prob;
    };

    //! Print equation list with policies
    //! \param[in] t Section title
    //! \param[in] factory Factory to get equation data from
    //! \param[in] ntypes Unique equation types
    //! \author J. Bakosi
    template< class Factory >
    void eqlist( const std::string& t,
                 const Factory& factory,
                 std::size_t ntypes ) const
    {
      if (!factory.empty()) {
        section( t );
        item( "Unique equation types", ntypes );
        item( "With all policy combinations", factory.size() );
        raw( '\n' );
        raw( m_item_indent + "Legend: equation name : supported policies\n" );
        raw( '\n' );
        raw( m_item_indent + "Policy codes:\n" );
        static_assert( tk::HasTypedefCode< kw::physics::info >::value,
                       "Policy code undefined for keyword" );
        static_assert( tk::HasTypedefCode< kw::problem::info >::value,
                       "Policy code undefined for keyword" );
        raw( m_item_indent + " * " + *kw::physics::code() + ": "
                           + kw::physics::name() + ":\n" );
        boost::mpl::for_each< ctr::Physics::keywords >( echoPolicies( this ) );
        raw( m_item_indent + " * " + *kw::problem::code() + ": "
                           + kw::problem::name() + ":\n" );
        boost::mpl::for_each< ctr::Problem::keywords >( echoPolicies( this ) );
        raw( '\n' );
        // extract eqname and supported policies for output
        const auto h = ctr::Physics();
        const auto r = ctr::Problem();
        std::map< std::string, Policies > eqs;      // eqname : policies
        for (const auto& f : factory)
          eqs[ PDEName( f.first ) ] +=
            Policies( h.code( f.first.template get< tag::physics >() ),
                      r.code( f.first.template get< tag::problem >() ) );
        // output eqname and supported policies
        for (const auto& e : eqs)
          m_stream << m_item_name_value_fmt % m_item_indent
                                            % e.first % e.second;
      }
    }

    //! Print configuration of a stack of partial differential equations
    void pdes( const std::string& t,
      const std::vector< std::vector< std::pair< std::string, std::string > > >&
        info ) const;

    //! Print time integration header
    void inthead( const std::string& t, const std::string& name,
                  const std::string& legend, const std::string& head ) const;

  private:
    //! Return partial differential equation name
    //! \param[in] key Equation key
    //! \return Partial differential equation name based on key
    template< class Key >
    std::string PDEName ( const Key& key ) const
    { return ctr::PDE().name( key.template get< tag::pde >() ); }
};

} // inciter::

#endif // InciterPrint_h
