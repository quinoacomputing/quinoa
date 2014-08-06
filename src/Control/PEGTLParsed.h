//******************************************************************************
/*!
  \file      src/Control/PEGTLParsed.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 10:08:08 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Class to equip parsed classes with PEGTL instruments
  \details   Class to equip parsed classes with PEGTL instruments
*/
//******************************************************************************
#ifndef PEGTLParsed_h
#define PEGTLParsed_h

namespace tk {
namespace ctr {

struct unused {};

//! PEGTLParsed
template< class Parsed,
          typename Input,
          typename cmdtag = unused,
          class Cmd = unused >
class PEGTLParsed : public Parsed {

  public:
    //! Constructor
    explicit PEGTLParsed( const Input& input ) : m_input(input) {}

    //! Constructor setting command line
    explicit PEGTLParsed( const Input& input, const Cmd& cl ) : m_input(input)
    { Parsed::template set< cmdtag >( cl ); }

    //! PEGTL location accessor
    const typename Input::location_type location() const
    { return m_input.location(); }

  private:
    const Input& m_input;      //!< Reference to PEGTL input parsed
};

} // ctr::
} // tk::

#endif // PEGTLParsed_h
