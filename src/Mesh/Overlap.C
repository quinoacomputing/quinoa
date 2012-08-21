// -----------------------------------------------------------------------------
// \file    src/Mesh/Overlap.C
// \author  jbakosi
// \date    Thu Aug 14 9:32:00 2012
// \brief   Triangle-triangle overlap test
// \note    Copyright 2012 Jozsef Bakosi
//          All rights reserved.
// -----------------------------------------------------------------------------
//
// Triangle-Triangle Overlap Test Routines
// July, 2002
// Updated December 2003
//
// This file contains C implementation of algorithms for
// performing two and three-dimensional triangle-triangle intersection test
// The algorithms and underlying theory are described in
//
//  "Fast and Robust Triangle-Triangle Overlap Test
//   Using Orientation Predicates"  P. Guigue - O. Devillers
//  Journal of Graphics Tools, 8(1), 2003
//
// Several geometric predicates are defined.  Their parameters are all
// points.  Each point is an array of two or three double precision
// floating point numbers. The geometric predicates implemented in
// this file are:
// 
//   int tri_tri_overlap_test_3d(p1,q1,r1,p2,q2,r2)
//   int tri_tri_overlap_test_2d(p1,q1,r1,p2,q2,r2)
//
//   int tri_tri_intersection_test_3d(p1,q1,r1,p2,q2,r2,
//                                    coplanar,source,target)
//
//      is a version that computes the segment of intersection when
//      the triangles overlap (and are not coplanar)
//
//   each function returns 1 if the triangles (including their
//   boundary) intersect, otherwise 0
//
//
// Other information are available from the Web page
// http:///www.acm.org/jgt/papers/GuigueDevillers03/
//
// -----------------------------------------------------------------------------
// Everything is taken out in order to compile faster, except
// tri_tri_overlap_test_2d() and related stuff.
// J, Bakosi, October, 7, 2006
// -----------------------------------------------------------------------------

#include "Overlap.h"

// Two dimensional Triangle-Triangle Overlap Test

#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))

#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
      if (ORIENT_2D(P1,P2,Q1) > 0.0f) {\
        if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1; \
        else return 0;} else {\
        if (ORIENT_2D(P1,P2,R1) >= 0.0f)\
          if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1; \
          else return 0;\
        else return 0;}\
    else \
      if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
        if (ORIENT_2D(R2,Q2,R1) <= 0.0f)\
          if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1; \
          else return 0;\
        else return 0;\
      else return 0;\
  else\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) \
      if (ORIENT_2D(Q1,R1,R2) >= 0.0f)\
        if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;\
        else return 0;\
      else \
        if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
          if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1; \
          else return 0; }\
        else return 0; \
    else  return 0; \
 };

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
        if (ORIENT_2D(P1,Q1,R2) >= 0.0f) return 1; \
        else return 0;} else { \
      if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
        if (ORIENT_2D(R1,P1,P2) >= 0.0f) return 1; else return 0;} \
      else return 0; } \
  } else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) {\
      if (ORIENT_2D(P1,P2,R1) >= 0.0f) {\
        if (ORIENT_2D(P1,R1,R2) >= 0.0f) return 1;  \
        else {\
          if (ORIENT_2D(Q1,R1,R2) >= 0.0f) return 1; else return 0;}}\
      else  return 0; }\
    else return 0; }}

static int ccw_tri_tri_intersection_2d(double p1[2], double q1[2], double r1[2],
                                       double p2[2], double q2[2], double r2[2])
{
  if ( ORIENT_2D(p2,q2,p1) >= 0.0f ) {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) return 1;
      else INTERSECTION_TEST_EDGE(p1,q1,r1,p2,q2,r2)
    } else {  
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
        INTERSECTION_TEST_EDGE(p1,q1,r1,r2,p2,q2)
      else INTERSECTION_TEST_VERTEX(p1,q1,r1,p2,q2,r2)}}
  else {
    if ( ORIENT_2D(q2,r2,p1) >= 0.0f ) {
      if ( ORIENT_2D(r2,p2,p1) >= 0.0f ) 
        INTERSECTION_TEST_EDGE(p1,q1,r1,q2,r2,p2)
      else  INTERSECTION_TEST_VERTEX(p1,q1,r1,q2,r2,p2)}
    else INTERSECTION_TEST_VERTEX(p1,q1,r1,r2,p2,q2)}
}

int tri_tri_overlap_test_2d(double p1[2], double q1[2], double r1[2], 
                            double p2[2], double q2[2], double r2[2])
{
  if ( ORIENT_2D(p1,q1,r1) < 0.0f )
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,r1,q1,p2,q2,r2);
  else
    if ( ORIENT_2D(p2,q2,r2) < 0.0f )
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,r2,q2);
    else
      return ccw_tri_tri_intersection_2d(p1,q1,r1,p2,q2,r2);
}
