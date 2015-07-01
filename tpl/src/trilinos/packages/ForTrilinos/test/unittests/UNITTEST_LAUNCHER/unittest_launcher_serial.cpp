#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>

#include "ForTrilinos_config.h"

//#define DEBUG

std::vector<std::string> get_tests(const char *filename)
{
  std::vector<std::string> names;

  const size_t maxlen = 1000;
  char *tmp = new char[maxlen];
  char *tmp1 = new char[maxlen];
  char *tmp2 = new char[maxlen];

#ifdef DEBUG
  std::cout << "Opening test list..." << std::endl;
#endif

  std::ifstream f(filename);
  if (!f.is_open())
    throw std::string("Error opening test list: ") + std::string(filename);

#ifdef DEBUG
  std::cout << "Loading test list..." << std::endl;
#endif

  int i=0;
  while (!f.eof()) {
    i++;
    f.getline(tmp, maxlen);
    int cnt = sscanf(tmp, "%s %s\n", tmp1, tmp2);
    if (cnt == 2) {
      if (tmp1[0] != '#') {
        std::string stmp = std::string(tmp1) + std::string(" ") + std::string(tmp2);
        names.push_back(stmp);
#ifdef DEBUG
        std::cout << "Adding test: " << stmp << std::endl;
#endif
      }
#ifdef DEBUG
      else {
        std::cout << "Skipping line " << i << ": " << tmp << std::endl;
      }
#endif
    }
#ifdef DEBUG
    else {
      std::cout << "Skipping line " << i << ": " << tmp << std::endl;
    }
#endif
  }

#ifdef DEBUG
  std::cout << "Done loading." << std::endl;
#endif

  f.close();

  delete [] tmp;
  delete [] tmp1;
  delete [] tmp2;

  return names;
}

bool run_unittest(const std::string & name)
{
  std::cout << "######################################################################" << std::endl;

  std::string cmd(name);

  std::cout << "Running test " << cmd << " ..." << std::endl;
  int ret = system(cmd.c_str());

  std::cout << "######################################################################" << std::endl;
  std::cout << std::endl;

  return (ret == 0);
}

int run_all_unittests(std::vector<std::string> &names)
{
  int failed = 0;

  size_t cnt = names.size();
  for (size_t i=0; i<cnt; i++) {
    bool success = run_unittest(names[i]);
    if (!success) failed += 1;
  }

  return failed;
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  std::cerr << "The serial unittest launcher is not capable of running MPI tests."
            << std::endl
            << "If you are seeing this error, then CMake isn't behaving as expected!"
            << std::endl
            << "END RESULT: SOME TESTS FAILED"
            << std::endl;
  return 1;
#endif

  int failed = 0;

  if (argc < 2) {
    std::cerr << "Specify unit test input file on the command line." << std::endl;
    return 1;
  }

  std::vector<std::string> names;

  try {
    names = get_tests(argv[1]);
  } catch (std::string &ex) {
    std::cerr << ex << std::endl;
    std::cerr << "END RESULT: SOME TESTS FAILED" << std::endl;
    return 1;
  }

  int cnt = names.size();
  failed = run_all_unittests(names);
  std::cout << (cnt-failed) << " tests passed out of " << cnt << std::endl;

  if (failed > 0) {
    std::cerr << "END RESULT: SOME TESTS FAILED" << std::endl;
  } else {
    std::cout << "END RESULT: ALL TESTS PASSED" << std::endl;
  }

  return (failed == 0 ? 0 : 1);
}
