#include <iostream>
#include <sstream>
#include "parser.hpp"
#include "exception.hpp"

namespace {

void proc_module_def(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_MODULE_DEF(CLASS_BEING_TESTED)
  if (args.size() != 1)
    throw ParserException("Incorrect number of arguments");
  std::stringstream ss;
  ss << "module " << args[0] << "_test_impls";
  result.push_back(ss.str());
}

void proc_module_begin(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_MODULE_BEGIN(CLASS_BEING_TESTED)
  if (args.size() != 1)
    throw ParserException("Incorrect number of arguments");
  result.push_back("implicit none");
  result.push_back("public");
  result.push_back("contains");
  result.push_back("logical function Does_Nothing_Returns_True() result(success)");
  result.push_back("success = .TRUE.");
  result.push_back("end function");
}

void proc_module_end(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_MODULE_END(CLASS_BEING_TESTED)
  if (args.size() != 1)
    throw ParserException("Incorrect number of arguments");
  std::stringstream ss;
  ss << "end module " << args[0] << "_test_impls";
  result.push_back(ss.str());
}

void proc_switch_begin(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_SWITCH_BEGIN(CLASS_BEING_TESTED)
  if (args.size() != 1)
    throw ParserException("Incorrect number of arguments");
}

void proc_switch_end(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_SWITCH_END(CLASS_BEING_TESTED)
  if (args.size() != 1)
    throw ParserException("Incorrect number of arguments");
}

void proc_unittest_def(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_DEF(CLASS_BEING_TESTED, ModuleName)
  if (args.size() != 2)
    throw ParserException("Incorrect number of arguments");
  std::stringstream ss;
  ss << "logical function " << args[0] << "_test_" << args[1] << "_UnitTest() result(success)";
  result.push_back(ss.str());
}

void proc_unittest_begin(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_BEGIN
  if (args.size() != 0)
    throw ParserException("Incorrect number of arguments");
}

void proc_unittest_end(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_END
  if (args.size() != 0)
    throw ParserException("Incorrect number of arguments");
  result.push_back("end function");
}

}

//////////////////////////////////////////////////////////////////////

void ImplParser::init()
{
  add_match("MODULE_DEF", &proc_module_def);
  add_match("MODULE_BEGIN", &proc_module_begin);
  add_match("MODULE_END", &proc_module_end);
  add_match("SWITCH_BEGIN", &proc_switch_begin);
  add_match("SWITCH_END", &proc_switch_end);
  add_match("DEF", &proc_unittest_def);
  add_match("BEGIN", &proc_unittest_begin);
  add_match("END", &proc_unittest_end);
}

