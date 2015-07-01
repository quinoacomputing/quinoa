#include <sstream>
#include "exception.hpp"

ParserException::ParserException(std::string err, size_t linenum,
    std::string file, std::string code)
  : _err(err), _linenum(linenum), _file(file), _code(code)
{
  std::stringstream ss;
  ss << _err << " at " << _file << ":" << _linenum;
  if (_code.size() > 0)
    ss << ":" << std::endl << _code;
  _msg = ss.str();
}

std::string ParserException::getErr()
{
  return _err;
}

std::string ParserException::getFile()
{
  return _file;
}

size_t ParserException::getLine()
{
  return _linenum;
}

std::string ParserException::getCode()
{
  return _code;
}

std::string ParserException::getMsg()
{
  return _msg;
}

