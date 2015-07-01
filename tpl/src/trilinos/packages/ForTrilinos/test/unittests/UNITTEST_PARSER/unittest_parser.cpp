#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <string>

#include "parser.hpp"
#include "exception.hpp"

bool read_template(const std::string &fname, std::vector<std::string> &input)
{
  bool success = true;
  const size_t buffer_size = 2048;
  char buffer[buffer_size];
  std::ifstream file;
  file.open(fname.c_str());
  if (!file.is_open()) {
    std::cerr << "Unable to open input file: " << fname << std::endl;
    success = false;
    return success;
  }
  try {
    while (!file.eof()) {
      file.getline(buffer, buffer_size);
      std::string line(buffer);
      input.push_back(line);
    }
  } catch (...) {
    std::cerr << "Error reading input: " << fname << std::endl;
    success = false;
  }
  file.close();
  return success;
}

bool replace_modname(const std::string &modname, std::vector<std::string> &input)
{ // don't change the line numbering here
  std::string modvar("CLASS_BEING_TESTED");
  for (size_t i=0; i<input.size(); i++) {
    std::string tmp = input[i];
    size_t loc = tmp.find(modvar);
    if (loc != tmp.npos) {
      tmp.replace(loc, modvar.size(), modname);
      input[i] = tmp;
    }
  }
  return true;
}

bool write_output(const std::string &fname, const std::vector<std::string> &output)
{
  bool success = true;
  std::ofstream file;
  file.open(fname.c_str());
  if (!file.is_open()) {
    std::cerr << "Unable to open input file: " << fname << std::endl;
    success = false;
    return success;
  }
  try {
    for (size_t i=0; i<output.size(); i++) {
//      std::cout << output[i] << std::endl;
      file << output[i] << std::endl;
    }
  } catch (...) {
    std::cerr << "Error writing output: " << fname << std::endl;
    success = false;
  }
  file.close();
  return success;
}

bool parse_impl(const std::string &ifname, const std::vector<std::string> &input, std::vector<std::string> &output)
{
  ImplParser p("FORTRILINOS_UNITTEST");
  p.init();
  std::vector<std::string> result;
  for (size_t i=0; i<input.size(); i++) {
    try {
      result = p.process_line(input[i]);
    } catch (ParserException &ex) {
      throw ParserException(ex.getErr(), i, ifname, input[i]);
    }
    for (size_t j=0; j<result.size(); j++) {
      output.push_back(result[j]);
    }
    result.clear();
  }
  return true;
}

bool parse_calls(const std::string &ifname, const std::vector<std::string> &input, std::vector<std::string> &output)
{
  CallsParser p("FORTRILINOS_UNITTEST");
  p.init();
  std::vector<std::string> result;
  for (size_t i=0; i<input.size(); i++) {
    try {
      result = p.process_line(input[i]);
    } catch (ParserException &ex) {
      throw ParserException(ex.getErr(), i, ifname, input[i]);
    }
    for (size_t j=0; j<result.size(); j++) {
      output.push_back(result[j]);
    }
    result.clear();
  }
  return true;
}

bool parse_list(const std::string &ifname, const std::vector<std::string> &input, std::vector<std::string> &output)
{
  ListParser p("FORTRILINOS_UNITTEST");
  p.init();
  std::vector<std::string> result;
  for (size_t i=0; i<input.size(); i++) {
    try {
      result = p.process_line(input[i]);
    } catch (ParserException &ex) {
      throw ParserException(ex.getErr(), i, ifname, input[i]);
    }
    for (size_t j=0; j<result.size(); j++) {
      output.push_back(result[j]);
    }
    result.clear();
  }
  return true;
}

bool parse_mpilist(const std::string &ifname, const std::vector<std::string> &input, std::vector<std::string> &output)
{
  MpiListParser p("FORTRILINOS_UNITTEST");
  p.init();
  std::vector<std::string> result;
  for (size_t i=0; i<input.size(); i++) {
    try {
      result = p.process_line(input[i]);
    } catch (ParserException &ex) {
      throw ParserException(ex.getErr(), i, ifname, input[i]);
    }
    for (size_t j=0; j<result.size(); j++) {
      output.push_back(result[j]);
    }
    result.clear();
  }
  return true;
}

int main(int argc, char *argv[])
{
  bool success = true;

  if (argc < 3) {
    std::cerr << "Syntax error: use" << std::endl;
    std::cerr << argv[0] << " classname inputpath" << std::endl;
    return 1;
  }

  std::string modname(argv[1]);
  std::string inputpath(argv[2]);
  std::string ifname = inputpath + std::string("/") + modname + std::string("_tests.un");
  std::string impl_fname = modname + std::string("_test_impls-tmp.F90");
  std::string calls_fname = modname + std::string("_test_calls.F90");
  std::string list_fname = modname + std::string("_tests.tests");
  std::string mpilist_fname = modname + std::string("_tests.mpitests");

  std::vector<std::string> input, output;

  try {
    success = success && read_template(ifname, input);
    success = success && replace_modname(modname, input);

    success = success && parse_impl(ifname, input, output);
    success = success && write_output(impl_fname, output);
    output.clear();

    success = success && parse_calls(ifname, input, output);
    success = success && write_output(calls_fname, output);
    output.clear();

    success = success && parse_list(ifname, input, output);
    success = success && write_output(list_fname, output);
    output.clear();

    success = success && parse_mpilist(ifname, input, output);
    success = success && write_output(mpilist_fname, output);
    output.clear();
  } catch (ParserException &ex) {
    std::cerr << "Caught exception: " << ex.getMsg() << std::endl;
    success = false;
  } catch (std::string &ex) {
    std::cerr << "Caught exception: " << ex << std::endl;
    success = false;
  } catch (...) {
    std::cerr << "Caught unknown exception." << std::endl;
    success = false;
  }

  return (success ? 0 : 1);
}

