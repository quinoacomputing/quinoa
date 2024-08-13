// *****************************************************************************
/*!
  \file      src/Main/InciterPrint.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
// *****************************************************************************
#ifndef InciterPrint_h
#define InciterPrint_h

#include <iostream>
#include <string>

#include "NoWarning/format.hpp"

#include "Print.hpp"
#include "ContainerUtil.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "Inciter/Options/Physics.hpp"
#include "Inciter/Options/Problem.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;
extern ctr::InputDeck g_inputdeck;

//! InciterPrint : tk::Print
class InciterPrint : public tk::Print {

  public:
    //! Constructor
    //! \param[in] screen Screen output filename
    //! \param[in,out] str Verbose stream
    //! \param[in] mode Open mode for screen output file, see
    //!   http://en.cppreference.com/w/cpp/io/ios_base/openmode
    //! \param[in,out] qstr Quiet stream
    //! \see tk::RNGPrint::RNGPrint and tk::Print::Print
    explicit InciterPrint( const std::string& screen,
                           std::ostream& str = std::clog,
                           std::ios_base::openmode mode = std::ios_base::out,
                           std::ostream& qstr = std::cout ) :
      Print( screen, str, mode, qstr ) {}

    //! Print control option: 'group : option'
    template< typename Option, typename... tags >
    void Item() const {
      Option opt;
      m_stream << m_item_name_value_fmt
                  % m_item_indent % opt.group()
                  % opt.name( g_inputdeck.get< tags... >() );
    }

    // Helper class for compact output of PDE policies
    class Policies {
      public:
        // Default constructor
        explicit Policies() : phys(), prob() {}
        // Initializer constructor
        explicit Policies( const std::string& p, const std::string& t ) :
          phys(p), prob(t) {}
        // Operator += for adding up two Policies structs
        Policies& operator+= ( const Policies& p ) {
          phys += p.phys;
          prob += p.prob;
          return *this;
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
    template< class Factory >
    void eqlist( const std::string& t,
                 const Factory& factory,
                 std::size_t ntypes ) const
    {
      if (!factory.empty()) {
        section( t );
        item( "Unique equation types", ntypes );
        item( "With all policy combinations", factory.size() );
      }
    }

    //! Print configuration of a stack of partial differential equations
    void pdes( const std::string& t,
      const std::vector< std::vector< std::pair< std::string, std::string > > >&
        info ) const;

    //! Print out info on solver coupling
    void couple( const std::vector< Transfer >& transfer,
                 const std::vector< char >& depvar ) const;

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
