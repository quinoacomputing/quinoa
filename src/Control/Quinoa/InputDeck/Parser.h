//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Sat 19 Oct 2013 07:20:53 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck file parser
  \details   Quinoa's input deck file parser
*/
//******************************************************************************
#ifndef QuinoaInputDeckParser_h
#define QuinoaInputDeckParser_h

#include <FileParser.h>
#include <Quinoa/CmdLine/CmdLine.h>
#include <Quinoa/InputDeck/InputDeck.h>

namespace quinoa {

//! InputDeckParser : FileParser
class InputDeckParser : public tk::FileParser {

  public:
    //! Constructor
    explicit InputDeckParser(const tk::Print& print,
                             std::unique_ptr< ctr::CmdLine > cmdline,
                             std::unique_ptr< ctr::InputDeck >& inputdeck);

    //! Destructor
    ~InputDeckParser() noexcept override = default;

  private:
    //! Don't permit copy constructor
    InputDeckParser(const InputDeckParser&) = delete;
    //! Don't permit copy assigment
    InputDeckParser& operator=(const InputDeckParser&) = delete;
    //! Don't permit move constructor
    InputDeckParser(InputDeckParser&&) = delete;
    //! Don't permit move assigment
    InputDeckParser& operator=(InputDeckParser&&) = delete;

    //! Make requested statistics unique
    void unique(std::vector<ctr::Product>& statistics);
};

} // namespace quinoa

#endif // QuinoaInputDeckParser_h
