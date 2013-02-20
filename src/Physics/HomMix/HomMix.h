//******************************************************************************
/*!
  \file      src/Physics/HomMix/HomMix.h
  \author    J. Bakosi
  \date      Tue 19 Feb 2013 09:41:24 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous material mix model
  \details   Homogeneous material mix model
*/
//******************************************************************************
#ifndef HomMix_h
#define HomMix_h

#include <map>
#include <sstream>

#include <Physics.h>
#include <ControlTypes.h>

namespace Quinoa {

class Memory;
class Paradigm;
class Mix;

//! HomMix : Physics
class HomMix : public Physics {

  public:
    //! Constructor
    HomMix(Memory* const memory,
           Paradigm* const paradigm,
           Control* const control);

    //! Destructor
    virtual ~HomMix();

    //! Echo informaion on model
    virtual void echo();

    //! Initialize model
    virtual void init();

    //! Solve model
    virtual void solve();

  private:
    //! Don't permit copy constructor
    HomMix(const HomMix&) = delete;
    //! Don't permit copy assigment
    HomMix& operator=(const HomMix&) = delete;
    //! Don't permit move constructor
    HomMix(HomMix&&) = delete;
    //! Don't permit move assigment
    HomMix& operator=(HomMix&&) = delete;

    //! Information needed when outputing a point in time history
    struct stamp {
      string filename;          //!< File name at time
      bool written;             //!< False if file is not yet written
      stamp(string f, bool w) : filename(f), written(w) {}
    };

    // History is a map of times and stamps
    using History = map<real,stamp>;

    //! Generate filenames based on times and return them as a map
    History history(const vector<real>& times, const string& base) {
      History n;
      for (auto& t : times) {
        stringstream ss;
        ss << base << "." << t;
        n.insert(pair<real,stamp>(t, stamp(ss.str(),false)));
      }
      return n;
    }

    //! Output joint scalar PDF
    void outJPDF();

    History m_JPDFHistory;            //!< Joint PDF save times and filenames
    Mix* m_mix;                       //!< Mix model object
};

} // namespace Quinoa

#endif // HomMix_h
