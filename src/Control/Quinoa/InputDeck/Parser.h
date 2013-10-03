//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Parser.h
  \author    J. Bakosi
  \date      Thu Oct  3 11:23:43 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck file parser
  \details   Quinoa's input deck file parser
*/
//******************************************************************************
#ifndef QuinoaInputDeckParser_h
#define QuinoaInputDeckParser_h

#include <vector>

#include <FileParser.h>
#include <Base.h>

namespace quinoa {

//! InputDeckParser : FileParser
class InputDeckParser : public FileParser {

  public:
    //! Constructor
    explicit InputDeckParser(Base& base);

    //! Destructor
    ~InputDeckParser() noexcept override = default;

    //! Parse quinoa control file
    void parse() override;

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
