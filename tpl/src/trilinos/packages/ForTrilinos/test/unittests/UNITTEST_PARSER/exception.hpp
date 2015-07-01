#ifndef PARSER_EXCEPTION_HPP
#define PARSER_EXCEPTION_HPP

class ParserException
{
  public:
    ParserException(std::string err, size_t linenum = 0,
        std::string file = "(unknown)", std::string code = "");
    ~ParserException() { }

    std::string getErr();
    std::string getFile();
    size_t getLine();
    std::string getCode();
    std::string getMsg();

  private:
    std::string _err;
    size_t _linenum;
    std::string _file;
    std::string _code;
    std::string _msg;
};

#endif /* PARSER_EXCEPTION_HPP */
