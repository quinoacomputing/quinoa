// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <unistd.h>

#include <stk_util/diag/PrintTimer.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/util/Writer.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/use_cases/UseCaseEnvironment.hpp>

using namespace use_case;

enum {
  TIMER_DOMAIN		= 0x00001000,		///< Enable domain timers
  TIMER_REGION		= 0x00002000,		///< Enable region timers
  TIMER_PROCEDURE	= 0x00004000,		///< Enable procedure timers
  TIMER_MECHANICS	= 0x00008000,		///< Enable mechanics timers
  TIMER_ALGORITHM	= 0x00010000,		///< Enable algorithm timers
  TIMER_SOLVER		= 0x00020000,		///< Enable solver timers
  TIMER_CONTACT		= 0x00040000,		///< Enable contact timers
  TIMER_MATERIAL	= 0x00080000,		///< Enable material timers
  TIMER_SEARCH		= 0x00100000,		///< Enable search timers
  TIMER_TRANSFER	= 0x00200000,		///< Enable transfer timers
  TIMER_ADAPTIVITY	= 0x00400000,		///< Enable adaptivity
  TIMER_UNUSED_1	= 0x00800000,		///< Enable unused 1
  TIMER_PROFILE_1	= 0x01000000,		///< Enable profile 1 timers
  TIMER_PROFILE_2	= 0x02000000,		///< Enable profile 2 timers
  TIMER_PROFILE_3	= 0x04000000,		///< Enable profile 3 timers
  TIMER_PROFILE_4	= 0x08000000,		///< Enable profile 4 timers
  TIMER_APP_1		= 0x10000000,		///< Enable application defined 1
  TIMER_APP_2		= 0x20000000,		///< Enable application defined 2
  TIMER_APP_3		= 0x40000000,		///< Enable application defined 3
  TIMER_APP_4		= 0x80000000,		///< Enable application defined 4
  TIMER_ALL		= 0xFFFFF000		///< Enable all timers
};

namespace {

double
work()
{
  double x = 1.0;

  for (int i = 0; i < 100000; ++i)
//  for (int i = 0; i < 100; ++i)
    x += std::sin(i);

  return x;
}

stk::diag::TimerSet &
unitTestTimerSet()
{
  static stk::diag::TimerSet s_unitTestTimerSet(TIMER_REGION);

  return s_unitTestTimerSet;
}


struct UseCaseRootTimer : public stk::diag::Timer
{
  UseCaseRootTimer()
    : stk::diag::Timer(stk::diag::createRootTimer("Unit test timer", unitTestTimerSet()))
  {}

  ~UseCaseRootTimer() {
    stk::diag::deleteRootTimer(*this);
  }
};


stk::diag::Timer &unitTestTimer() {
  static UseCaseRootTimer s_unitTestTimer;

  return s_unitTestTimer;
}


struct RootObject
{
  RootObject()
    : m_timer("Root object", TIMER_REGION, unitTestTimer())
  {}

  stk::diag::Timer      m_timer;
};


struct Object
{
  Object(const std::string &name, RootObject &root_object)
    : m_id(0),
      m_name(name),
      m_timer(name, root_object.m_timer)
  {}

  Object(int id, const Object &parent)
    : m_id(id),
      m_name(id_name(id)),
      m_timer(m_name, parent.m_timer)
  {}

  static std::string id_name(int id) {
    std::ostringstream s;
    s << "Object id " << id << " run";
    return s.str();
  }

  void run() {
    stk::diag::TimeBlock _time(m_timer);
    m_x = work();
  }

  int                   m_id;
  double                m_x;
  std::string           m_name;
  stk::diag::Timer      m_timer;
};

} // namespace <empty>


int
use_case_timer(
  UseCaseEnvironment &  use_case_environment)
{
  int parallel_rank = stk::parallel_machine_rank(use_case_environment.m_comm);

  stk::diag::TimeBlock root_time_block(unitTestTimer());
  
  stk::diag::updateRootTimer(unitTestTimer());

  // Check lap count
  {
    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = unitTestTimer().getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    ThrowRequire(lap_count == 1);
  }
  
  // Create subtimer and run 100 laps
  {
    static stk::diag::Timer run_timer("Run 100 times twice", TIMER_REGION, unitTestTimer());

    for (int i = 0; i < 100; ++i) {
      stk::diag::TimeBlock _time(run_timer);
      work();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = run_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    ThrowRequire(lap_count == 100);
  }

  // Grab previous subtimer and run 100 laps
  {
    static stk::diag::Timer run_timer("Run 100 times twice", TIMER_REGION, unitTestTimer());

    for (int i = 0; i < 100; ++i) {
      stk::diag::TimeBlock _time(run_timer);
      work();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = run_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    ThrowRequire(lap_count == 200);
  }

  // Create root object
  RootObject root_object;

  {
    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = root_object.m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    ThrowRequire(lap_count == 0);
  }

  // Create object
  {
    Object time_object("One object", root_object);

    for (int i = 0; i < 100; ++i) {
      time_object.run();
    }

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = time_object.m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    ThrowRequire(lap_count == 100);
  }

  // Create object tree
  {
    std::vector<Object> object_vector;
    object_vector.push_back(Object("Object Tree", root_object));

    int id = 0;
    for (size_t i = 0; i < 2; ++i) {
      size_t ix = object_vector.size();
      object_vector.push_back(Object(id++, object_vector[0]));
      for (size_t j = 0; j < 2; ++j) {
        size_t jx = object_vector.size();
        object_vector.push_back(Object(id++, object_vector[ix]));
        for (int k = 0; k < 2; ++k) {
          object_vector.push_back(Object(id++, object_vector[jx]));
        }
      }
    }

    sierra::out() << "One object run 100 times" << std::endl;

    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, false);

    stk::diag::MetricTraits<stk::diag::LapCount>::Type lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j)
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    ThrowRequire(lap_count == stk::diag::MetricTraits<stk::diag::LapCount>::Type(0));

    for (size_t j = 0; j < object_vector.size(); ++j)
      object_vector[j].run();

    sierra::out() << "Object tree run once" << std::endl;

    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, false);

    lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j)
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    ThrowRequire(lap_count == stk::diag::MetricTraits<stk::diag::LapCount>::Type(object_vector.size()));

    for (size_t i = 1; i < 100; ++i)
      for (size_t j = 0; j < object_vector.size(); ++j)
        object_vector[j].run();

    sierra::out() << "Object tree 100 times (checkpointed)" << std::endl;

    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, true);

    lap_count = 0;
    for (size_t j = 0; j < object_vector.size(); ++j)
      lap_count += object_vector[j].m_timer.getMetric<stk::diag::LapCount>().getAccumulatedLap(false);

    ThrowRequire(lap_count == stk::diag::MetricTraits<stk::diag::LapCount>::Type(100*object_vector.size()));

    for (size_t i = 0; i < 100; ++i)
      for (size_t j = 0; j < object_vector.size(); ++j)
        object_vector[j].run();

    sierra::out() << "Object tree 100 more times (checkpointed)" << std::endl;

    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, true);

    sierra::out() << "Object tree (not checkpointed)" << std::endl;
    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, false);

    // Add a static function timer on odd processors
    if (parallel_rank%2 == 1) {
      static stk::diag::Timer s_functionTimer("Odd processor function timer", TIMER_REGION, unitTestTimer());
      {
        stk::diag::TimeBlock _time(s_functionTimer);
        work();
      }
      if (parallel_rank%3 == 1) {
        static stk::diag::Timer s_functionTimer2("Odd processor function timer 2", TIMER_REGION, s_functionTimer);
        {
          stk::diag::TimeBlock _time(s_functionTimer2);
          work();
        }
      }
    }

    sierra::out() << "Add asymetric timers on odd processors, parallel collected output (checkpointed)" << std::endl;
    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, true);

    for (size_t i = 0; i < 100; ++i)
      for (size_t j = 0; j < object_vector.size(); ++j)
        object_vector[j].run();

    sierra::out() << "Object tree 100 more times, parallel collected output (checkpointed)" << std::endl;
    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, true, MPI_COMM_WORLD);

    sierra::out() << "Object tree, parallel collected output (not checkpointed)" << std::endl;
    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);

    stk::diag::setTimerTimeFormat(stk::TIMEFORMAT_SECONDS);
    
    sierra::out() << "Object tree, parallel collected output (not checkpointed), seconds" << std::endl;
    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);

    stk::diag::setTimerTimeFormat(stk::TIMEFORMAT_SECONDS | stk::TIMEFORMAT_MILLIS);
    
    sierra::out() << "Object tree, parallel collected output (not checkpointed), seconds and milliseconds" << std::endl;
    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);

    stk::diag::setTimerTimeFormat(stk::TIMEFORMAT_HMS);
    
    sierra::out() << "Object tree, parallel collected output (not checkpointed), hh:mm:ss" << std::endl;
    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);

    stk::diag::setTimerTimeFormat(stk::TIMEFORMAT_HMS | stk::TIMEFORMAT_MILLIS);
    
    sierra::out() << "Object tree, parallel collected output (not checkpointed), hh:mm:ss.mmm" << std::endl;
    stk::diag::printTimersTable(sierra::out(), unitTestTimer(), stk::diag::METRICS_ALL, false, MPI_COMM_WORLD);
  }

  return 0;
}
