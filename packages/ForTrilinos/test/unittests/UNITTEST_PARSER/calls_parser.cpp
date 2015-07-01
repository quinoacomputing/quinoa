#include <iostream>
#include <sstream>
#include "parser.hpp"
#include "exception.hpp"

namespace {

static bool ignore_now = false;

void proc_module_def(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_MODULE_DEF(CLASS_BEING_TESTED)
  if (args.size() != 1)
    throw ParserException("Incorrect number of arguments");
  {
    std::stringstream ss;
    ss << "module " << args[0] << "_test_calls";
    result.push_back(ss.str());
  }
  {
    std::stringstream ss;
    ss << "use " << args[0] << "_test_impls";
    result.push_back(ss.str());
  }
  ignore_now = true;
}

void proc_module_begin(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_MODULE_BEGIN(CLASS_BEING_TESTED)
  if (args.size() != 1)
    throw ParserException("Incorrect number of arguments");
  result.push_back("use iso_fortran_env ,only : error_unit ,output_unit");
  result.push_back("implicit none");
  result.push_back("public");
  result.push_back("contains");
  result.push_back("logical function select_test(which_test) result(success)");
  result.push_back("character(len=50),intent(in) :: which_test");
  result.push_back("success=.FALSE.");
  result.push_back("if (.FALSE.) then");
  result.push_back("! ignore me");
}

void proc_module_end(const std::string &line, const std::string &term,
  const std::vector<std::string> &args, std::vector<std::string> &result)
{
//FORTRILINOS_UNITTEST_MODULE_END(CLASS_BEING_TESTED)
  if (args.size() != 1)
    throw ParserException("Incorrect number of arguments");
  result.push_back("else");
  result.push_back("write(output_unit,fmt='(a)') \"Missing test. TEST FAILED\"");
  result.push_back("endif");
  result.push_back("end function");
  std::stringstream ss;
  ss << "end module " << args[0] << "_test_calls";
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
  {
    std::stringstream ss;
    ss << "else if (which_test==\"" << args[1] << "\") then";
    result.push_back(ss.str());
  }
  {
    std::stringstream ss;
    ss << "success=" << args[0] << "_test_" << args[1] << "_UnitTest()";
    result.push_back(ss.str());
  }
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
}

}

//////////////////////////////////////////////////////////////////////

void CallsParser::init()
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

void CallsParser::proc_plain(const std::string &line, std::vector<std::string> &result)
{
  if (!ignore_now) {
    if (!is_blank(line))
      result.push_back(line);
  }
}

