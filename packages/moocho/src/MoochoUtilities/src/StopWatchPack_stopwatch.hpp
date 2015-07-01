// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STPWATCH_H
#define STPWATCH_H

#include "Moocho_ConfigDefs.hpp"

namespace StopWatchPack {

/** @name namespace StopWatchPack
  *
  * @memo Package for CPU timing.
  */
//@{

double seconds(void);

/** \brief Simple stopwatch object.
  */
class stopwatch {
public:

  /// Initializes of not running.
  stopwatch() : running_(false), last_time_(0.0), total_(0.0)
  {}

  /// Returns true if <tt>this</tt> is currently timming.
  bool is_running() const {
    return running_;
  }

  /// Starts timing if it has already not been started.
  void start() {
    if (!running_) {
      last_time_ = seconds();
      running_ = true;
    }
  }

  /// Stops timing and returns the time (sec.) since start() was called
  double stop()  {
    if (running_) {
      total_ += seconds() - last_time_; 
      running_ = false;
    }
    //std::cout << "total = " << total_ << std::endl;
    return total_; 
  }

  /// Stops and resets the clock if it is running.
  void reset() {
    running_ = false;
    last_time_ = 0.0;
    total_ = 0.0;
  }

  /// Reads the elapsed time (sec.) and leaves the clock running.
  double read() {
    if (running_) {
      double curr_time = seconds();
      total_ += curr_time - last_time_;
      last_time_ = curr_time;
    }
    return total_;
  }
  
private:
  bool running_;
  double last_time_;
  double total_;
};

//	end namespace StopWatchPack 
//@}

}  // end namespace StopWatchPack 

#endif // STPWATCH_H
