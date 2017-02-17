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
// Program to generate uniform-grid input files for Zoltan's zdrive program.
// Compile:
//    c++ generate_miniFElike_grids.cpp
// Run:
//    a.out x y z filename
// where x,y,z are the number of grid points in the x, y, and z directions,
// respectively, and filename is the base filename to be used for the output.
// 
// Output files are in Chaco and Matrix-Market format.
// Output files include filename.coords, filename.graph and filename.mtx.
//
// Note that x, y, and z are the number of GRID POINTS, not elements.  So
// to get input equivalent to miniFE nx=100 ny=100 nz=100, you must specify
// x = y = z = 101 for this generator.
//
// The graph is node connectivity within elements:  assume a finite
// element interpretation, with linear elements; then grid point j is connected
// to each grid point in all cells containing j.

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;


void compute_neighbors(
  ofstream &gfp,
  ofstream &mfp,
  int64_t cnt,
  int i,  int j,  int k,
  int nx, int ny, int nz
)
{
#define TEST_CONVERT_PRINT(i,j,k) \
  if ((i) >= 0 && (i) < nx && (j) >= 0 && (j) < ny && (k) >= 0 && (k) < nz) { \
    int64_t tmp = 1 + (k)*face + (j)*nx + (i); \
    gfp << tmp << " "; \
    mfp << setw(10) << cnt << setw(10) << tmp << setw(10) << "1." << endl; \
    ecnt++; \
  }

  static int64_t ecnt = 0;
  int64_t face = nx * ny;
  // Self
  TEST_CONVERT_PRINT(i,   j,   k)
  TEST_CONVERT_PRINT(i-1, j,   k)
  TEST_CONVERT_PRINT(i+1, j,   k)
  TEST_CONVERT_PRINT(i,   j-1, k)
  TEST_CONVERT_PRINT(i-1, j-1, k)
  TEST_CONVERT_PRINT(i+1, j-1, k)
  TEST_CONVERT_PRINT(i,   j+1, k)
  TEST_CONVERT_PRINT(i-1, j+1, k)
  TEST_CONVERT_PRINT(i+1, j+1, k)
  TEST_CONVERT_PRINT(i,   j,   k-1)
  TEST_CONVERT_PRINT(i-1, j,   k-1)
  TEST_CONVERT_PRINT(i+1, j,   k-1)
  TEST_CONVERT_PRINT(i,   j-1, k-1)
  TEST_CONVERT_PRINT(i-1, j-1, k-1)
  TEST_CONVERT_PRINT(i+1, j-1, k-1)
  TEST_CONVERT_PRINT(i,   j+1, k-1)
  TEST_CONVERT_PRINT(i-1, j+1, k-1)
  TEST_CONVERT_PRINT(i+1, j+1, k-1)
  TEST_CONVERT_PRINT(i,   j,   k+1)
  TEST_CONVERT_PRINT(i-1, j,   k+1)
  TEST_CONVERT_PRINT(i+1, j,   k+1)
  TEST_CONVERT_PRINT(i,   j-1, k+1)
  TEST_CONVERT_PRINT(i-1, j-1, k+1)
  TEST_CONVERT_PRINT(i+1, j-1, k+1)
  TEST_CONVERT_PRINT(i,   j+1, k+1)
  TEST_CONVERT_PRINT(i-1, j+1, k+1)
  TEST_CONVERT_PRINT(i+1, j+1, k+1)
  gfp << endl;
#undef TEST_CONVERT_PRINT
}


int main(int narg, char *arg[])
{

  if (narg != 5) {
    cout << "Usage:  a.out x y z filename" << endl;
    exit(-1);
  }

  int nx = atoi(arg[1]);
  int ny = atoi(arg[2]);
  int nz = atoi(arg[3]);
  int64_t nvtx = nx * ny * nz;

  int nx_gt_one = (nx > 1);
  int ny_gt_one = (ny > 1);
  int nz_gt_one = (nz > 1);
  int nx_m_two = nx - 2;
  int ny_m_two = ny - 2;
  int nz_m_two = nz - 2;

  int64_t nedge = 
      // interior vertices
      27 * nx_gt_one * ny_gt_one * nz_gt_one * nx_m_two * ny_m_two * nz_m_two
      // face vertices
    + 9 * (1 + 3 * nz_gt_one) * nx_gt_one * ny_gt_one * nx_m_two * ny_m_two
    + 9 * (1 + 3 * ny_gt_one) * nx_gt_one * nz_gt_one * nx_m_two * nz_m_two
    + 9 * (1 + 3 * nx_gt_one) * ny_gt_one * nz_gt_one * ny_m_two * nz_m_two
      // edge vertices
    + 12 * (ny_gt_one+nz_gt_one) * (ny_gt_one+nz_gt_one) * nx_gt_one * nx_m_two
    + 12 * (nx_gt_one+nz_gt_one) * (nx_gt_one+nz_gt_one) * ny_gt_one * ny_m_two
    + 12 * (nx_gt_one+ny_gt_one) * (nx_gt_one+ny_gt_one) * nz_gt_one * nz_m_two
      // corner vertices
    + (1+nx_gt_one) * (1+nx_gt_one) * (1+ny_gt_one) * (1+ny_gt_one) 
       * (1+nz_gt_one) * (1+nz_gt_one);

  char *basename = arg[4];

  // Open output files
  char fname[127];

  ofstream cfp;  // coordinates file
  sprintf(fname, "%s.coords", basename);
  cfp.open(fname);

  ofstream gfp;  // graph file
  sprintf(fname, "%s.graph", basename);
  gfp.open(fname);
  gfp << nvtx << " " << nedge << endl;

  ofstream mfp;  // matrix-market file
  sprintf(fname, "%s.mtx", basename);
  mfp.open(fname);
  mfp << "%%MatrixMarket matrix coordinate real general\n%%\n";
  mfp << nvtx << " " << nvtx << " " << nedge << endl;

  // subdivide the unit cube into cells
  //double dx = (nx > 1 ? 1. / (nx - 1) : 0.);
  //double dy = (ny > 1 ? 1. / (ny - 1) : 0.);
  //double dz = (nz > 1 ? 1. / (nz - 1) : 0.);
  // miniFE has dx = 1, dy = 1, dx = -1
  double dx = 1.;
  double dy = 1.;
  double dz = -1.;

  int64_t cnt = 0;
  double x = 0., y = 0., z = 0.;
  for (int k = 0; k < nz; k++) {
    y = 0.;
    for (int j = 0; j < ny; j++) {
      x = 0.;
      for (int i = 0; i < nx; i++) {
        cnt++;
        // print neighbors
        compute_neighbors(gfp, mfp, cnt, i, j, k, nx, ny, nz);
        // print coordinates
        cfp << setw(15) << x << setw(15) << y << setw(15) << z << endl;
        x += dx;
      }
      y += dy;
    }
    z += dz;
  }

  cfp.close();
  mfp.close();
}
