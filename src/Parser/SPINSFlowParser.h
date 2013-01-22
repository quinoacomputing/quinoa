//******************************************************************************
/*!
  \file      src/Parser/SPINSFlowParser.h
  \author    J. Bakosi
  \date      Mon 21 Jan 2013 09:10:16 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SPINSFlow parser
  \details   SPINSFlow parser
*/
//******************************************************************************
#ifndef SPINSFlowParser_h
#define SPINSFlowParser_h

using namespace std;

namespace Quinoa {

//! SPINSFlowParser : Parser
class SPINSFlow : private Parser {

  public:
    //! Constructor
    SPINSFlowParser();

    //! Destructor
    virtual ~SPINSFlowParser();

  private:
    //! Don't permit copy constructor
    SPINSFlowParser(const SPINSFlowParser&) = delete;
    //! Don't permit copy assigment
    SPINSFlowParser& operator=(const SPINSFlowParser&) = delete;
    //! Don't permit move constructor
    SPINSFlowParser(SPINSFlowParser&&) = delete;
    //! Don't permit move assigment
    SPINSFlowParser& operator=(SPINSFlowParser&&) = delete;

    // Include token IDs
    #include <SPINSFlow.def.h>

//     template <typename Lexer>
//     struct word_count_tokens : lex::lexer<Lexer> {
//       word_count_tokens() {
//           // define tokens (the regular expression to match and the corresponding
//           // token id) and add them to the lexer 
//           this->self.add
//             ("[^ \t\n]+", ID_WORD) // words (anything except ' ', '\t' or '\n')
//             ("\n", ID_EOL)         // newline characters
//             (".", ID_CHAR)         // anything else is a plain character
//           ;
//       }
//     };
};

} // namespace Quinoa

#endif // SPINSFlowParser_h
