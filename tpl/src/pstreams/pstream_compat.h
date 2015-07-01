/* $Id: pstream_compat.h,v 1.7 2008/07/07 21:38:14 redi Exp $
PStreams - POSIX Process I/O for C++
Copyright (C) 2001,2002,2007,2008 Jonathan Wakely

This file is part of PStreams.

PStreams is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

PStreams is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @file pstream_compat.h
 * Declares non-standard implementations of the PStreams classes for older
 * compilers.
 */

#ifndef REDI_PSTREAM_COMPAT_H_SEEN
#define REDI_PSTREAM_COMPAT_H_SEEN

#if ! ( BACK_COMPAT == 1 || GCC_BACK_COMPAT == 1 )

#error \
    You must define either BACK_COMPAT or GCC_BACK_COMPAT to be 1 \
    to use this file. For gcc versions 2.7 to 2.9x (including egcs) \
    define GCC_BACK_COMPAT = 1. For other compilers define BACK_COMPAT = 1 \
    For ISO C++ conforming compilers use pstream.h instead.

#elif GCC_BACK_COMPAT == 1

#warning "PStreams needs an ISO C++ compliant compiler."
#warning "These classes are unsupported and largely untested."

// turn on non-standard extensions
#define _STREAM_COMPAT 1

#include <fstream>
#include <stdio.h>
#include <errno.h>
#include <string.h>

namespace redi
{
  class pstreambase : public fstreambase
  {
  public:
    void open(const char* command)
    {
      if ((fp_ = popen(command, mode_)))
      {
        attach(fileno(fp_));
        // errno is set to zero because fstreambase tries to lseek()
        // to the end of the stream which is not allowed for pipes
        // and sets errno to ESPIPE, so we clear it.
        errno=0;
      }
    }

    int close() { return fp_ ? rdbuf()->sync(), pclose(fp_) : 0; }

    bool is_open() { return fp_; }

  protected:
    pstreambase(const char* mode) : fstreambase(), mode_(mode), fp_(0) { }

    pstreambase(const char* command, const char* mode)
    : fstreambase(), mode_(mode), fp_(0)
    {
      open(command);
    }

    ~pstreambase() { close(); }

  private:
    const char* const mode_; // Points to static mode member of derived class.
    FILE*             fp_;   // Preserves return value of popen().
  };


  class ipstream : public pstreambase, public istream
  {
    static const char * const mode = "r";
  public:
    ipstream() : pstreambase(mode) { }
    ipstream(const char* command) : pstreambase(command, mode) { }
  };

  class opstream : public pstreambase, public ostream
  {
    static const char * const mode = "w";
  public:
    opstream() : pstreambase(mode) { }
    opstream(const char* command) : pstreambase(command, mode) { }
  };

#if REDI_PSTREAMS_POPEN_USES_BIDIRECTIONAL_PIPE
  class pstream :  protected pstreambase, public iostream {
    static const char * const mode = "r+";
  public:
    pstream() : pstreambase(mode) { }
    pstream(const char* command) : pstreambase(command, mode) { }
  };
#endif

}  // namespace redi

#elif BACK_COMPAT == 1

#warning "PStreams needs an ISO C++ compliant compiler!"
#warning "These classes are unsupported and untested."

#include <fstream>
#include <cstdio>
#include <cerrno>

  // very basic implementations
  // incomplete interfaces: no open(const char*)
  // use non-standard extensions
  // use at your own risk!

namespace redi
{
  struct FILE_wrapper {
    explicit FILE_wrapper(FILE* p) : fp(p) {}
    FILE* fp;
  };

  class opstream : private FILE_wrapper, public std::ofstream {
    static const char * const mode = "w";
  public:
    opstream(const char* command)
    : FILE_wrapper(popen(command, mode)), std::ofstream(fileno(fp))
    {
      errno=0;
    }

    ~opstream() { pclose(fp); }
  };

  class ipstream : private FILE_wrapper, public std::ifstream {
    static const char * const mode = "r";
  public:
    ipstream(const char* command)
    : FILE_wrapper(popen(command, mode)), std::ifstream(fileno(fp))
    {
      errno=0;
    }

    ~ipstream() { pclose(fp); }
  };

#if REDI_PSTREAMS_POPEN_USES_BIDIRECTIONAL_PIPE
  class pstream : private FILE_wrapper, public std::fstream {
    static const char * const mode = "r+";
  public:
    pstream(const char* command)
    : FILE_wrapper(popen(command, mode)), std::fstream(fileno(fp))
    {
      errno=0;
    }
    
    ~pstream() { pclose(fp); }
  };
#endif

}  // namespace redi


#endif  // __GNUC__

#endif // REDI_PSTREAM_COMPAT_H_SEEN

