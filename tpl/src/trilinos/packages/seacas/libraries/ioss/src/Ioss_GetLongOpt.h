/* S Manoharan. Advanced Computer Research Institute. Lyon. France */

#ifndef _GetLongOption_h_
#define _GetLongOption_h_

#include <iostream>

namespace Ioss {
  /** \brief A database of program command line and environment variable options and methods for
   * manipulating them.
   *
   *  A collection of long command line option names for a program that uses the Ioss library.
   */
  class GetLongOption
  {
  public:
    enum OptType { NoValue, OptionalValue, MandatoryValue };

  private:
    struct Cell
    {
      const char *option;      // option name
      OptType     type;        // option type
      const char *description; // a description of option
      const char *value;       // value of option (string)
      const char *opt_value;   // If optional value and value not entered, assign opt_value to value
      Cell *      next;        // pointer to the next cell

      Cell()
      {
        option = description = value = opt_value = nullptr;
        next                                     = nullptr;
        type                                     = NoValue;
      }
    };

  private:
    Cell *      table;       // option table
    const char *ustring;     // usage message
    char *      pname;       // program basename
    Cell *      last;        // last entry in option table
    int         enroll_done; // finished enrolling
    char        optmarker;   // option marker

  private:
    int setcell(Cell *c, char *valtoken, char *nexttoken, const char *name);

  public:
    explicit GetLongOption(const char optmark = '-');
    ~GetLongOption();

    static char *basename(char *const pathname);

    int parse(int argc, char *const *argv);
    int parse(char *const str, char *const p);

    int enroll(const char *const opt, const OptType t, const char *const desc,
               const char *const val, const char *const optval = nullptr);
    const char *retrieve(const char *const opt) const;

    void usage(std::ostream &outfile = std::cout) const;

    /** \brief Set the program usage string.
     *
     *  The program usage string should define the command line
     *  syntax for program options and arguments and contain
     *  other helpful usage text.
     *  \param[in] str The usage string.
     */
    void usage(const char *str) { ustring = str; }
  };
}
#endif /* _GetLongOption_h_ */
