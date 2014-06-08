//******************************************************************************
/*!
  \file      src/Control/PEGTLParsed.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 08:08:32 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Class to equip parsed classes with PEGTL instruments
  \details   Class to equip parsed classes with PEGTL instruments
*/
//******************************************************************************
#ifndef PEGTLParsed_h
#define PEGTLParsed_h

namespace quinoa {
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
    { Parsed::template set< cmdtag >( cl );}

    //! PEGTL location accessor
    const typename Input::location_type location() const
    { return m_input.location(); }

  private:
    //! Don't permit copy constructor
    PEGTLParsed(const PEGTLParsed&) = delete;
    //! Don't permit copy assigment
    PEGTLParsed& operator=(const PEGTLParsed&) = delete;
    //! Don't permit move constructor
    PEGTLParsed(PEGTLParsed&&) = delete;
    //! Don't permit move assigment
    PEGTLParsed& operator=(PEGTLParsed&&) = delete;

    const Input& m_input;      //!< Reference to PEGTL input parsed
};

} // ctr::
} // quinoa::

#endif // PEGTLParsed_h
