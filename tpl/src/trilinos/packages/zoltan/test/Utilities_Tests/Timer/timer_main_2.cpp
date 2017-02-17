/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include "zoltan_timer_cpp.h"
using namespace std;

/* Test program to exercise Timer capabilities. */

/****************************************************************************/

void first_test(Zoltan_Timer_Object *zt)
{
/* First test of Timer:  This test accrues times through
 * separate calls to first_test.  
 * The time for timer two should be roughly twice that of timer one.
 * The time for timer three should be roughly four times that of timer one.
 */

int i, j, me; 
static int firsttime=1;
const int LOOP1=1000,
          LOOP2=2000,
          LOOP3=4000;
const int MAINLOOP=100;
const int USE_BARRIER=1;
static int t1=-1, t2=-1, t3=-1;

  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  for (i = 0; i < MAINLOOP; i++) {

    if (firsttime)
      t1 = zt->Init(USE_BARRIER, "Loop 1");

    zt->Start(t1, MPI_COMM_WORLD);
    for (j = 0; j < LOOP1; j++) {
      double a;
      a = sqrt((double) (j * LOOP1));
    }
    zt->Stop(t1, MPI_COMM_WORLD);

    if (firsttime)
      t2 = zt->Init(USE_BARRIER, "Loop 2");

    zt->Start(t2, MPI_COMM_WORLD);
    for (j = 0; j < LOOP2; j++) {
      double a;
      a = sqrt((double) (j * LOOP2));
    }
    zt->Stop(t2, MPI_COMM_WORLD);

    if (firsttime)
      t3 = zt->Init(USE_BARRIER, "Loop 3");

    zt->Start(t3, MPI_COMM_WORLD);
    for (j = 0; j < LOOP3; j++) {
      double a;
      a = sqrt((double) (j * LOOP3));
    }
    zt->Stop(t3, MPI_COMM_WORLD);

    firsttime=0;
  }

  zt->PrintAll(me, MPI_COMM_WORLD, stdout);

}

/****************************************************************************/

void second_test(Zoltan_Timer_Object *zt)
{
/* Second test of Timer:  This test does not accrue times through
 * separate function calls.  It exercises the REALLOC when more timers
 * than INITLENGTH are requested. 
 * Computation is similar to first test.  Timer_flag used is different.
 */
int i, j, me; 
const int LOOP1=1000,
          LOOP2=2000,
          LOOP3=4000;
const int MAINLOOP=100;
const int USE_BARRIER=1;
int t1=-1, t2=-1, t3=-1;
char str[200];
static int cnt = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  sprintf(str, "STLoop 1 %d", cnt);
  t1 = zt->Init(USE_BARRIER, str);
  sprintf(str, "STLoop 2 %d", cnt);
  t2 = zt->Init(USE_BARRIER, str);
  sprintf(str, "STLoop 3 %d", cnt);
  t3 = zt->Init(USE_BARRIER, str);
  cnt++;

  for (i = 0; i < MAINLOOP; i++) {
    zt->Start(t1, MPI_COMM_WORLD);
    for (j = 0; j < LOOP1; j++) {
      double a;
      a = sqrt((double) (j * LOOP1));
    }
    zt->Stop(t1, MPI_COMM_WORLD);

    zt->Start(t2, MPI_COMM_WORLD);
    for (j = 0; j < LOOP2; j++) {
      double a;
      a = sqrt((double) (j * LOOP2));
    }
    zt->Stop(t2, MPI_COMM_WORLD);

    zt->Start(t3, MPI_COMM_WORLD);
    for (j = 0; j < LOOP3; j++) {
      double a;
      a = sqrt((double) (j * LOOP3));
    }
    zt->Stop(t3, MPI_COMM_WORLD);
  }

  zt->PrintAll(me, MPI_COMM_WORLD, stdout);
}

/****************************************************************************/
main(int argc, char *argv[])
{
Zoltan_Timer_Object zt1(ZOLTAN_TIME_WALL), zt2(ZOLTAN_TIME_CPU);
int i, me;
const int MAINLOOP=20;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  for (i = 0; i < MAINLOOP; i++) {
    cout << "\n\n\t****Beginning first test****\n";
    first_test(&zt1);

    cout << "\n\n\t****Beginning second test****\n";
    second_test(&zt2);
  }

  cout << "\n\nFINAL RESULTS -- FIRST TEST:\n";
  zt1.PrintAll(me, MPI_COMM_WORLD, stdout);
  cout << "\n\nFINAL RESULTS -- SECOND TEST:\n";
  zt2.PrintAll(me, MPI_COMM_WORLD, stdout);

  //  Copy tests
  Zoltan_Timer_Object zt3 = Zoltan_Timer_Object(zt1);
  Zoltan_Timer_Object zt4 = zt2;
  for (i = 0; i < MAINLOOP; i++) {
    cout << "\n\n\t****Beginning first copy test****\n";
    first_test(&zt3);

    cout << "\n\n\t****Beginning second copy test****\n";
    second_test(&zt4);
  }

  cout << "\n\nFINAL RESULTS -- FIRST TEST:\n";
  zt3.PrintAll(me, MPI_COMM_WORLD, stdout);
  cout << "\n\nFINAL RESULTS -- SECOND TEST:\n";
  zt4.PrintAll(me, MPI_COMM_WORLD, stdout);


  cout << "\n\nTHE END\n";

  MPI_Finalize();
}

