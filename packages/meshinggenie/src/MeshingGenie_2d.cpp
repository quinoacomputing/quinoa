// 09/01/2011 5:41 pm

// R1.0: Non-convex domains (sampling + Voronoi meshing + output and testing planar faces + non-convex domains + internal faces)

// Author: Mohamed S. Ebeida

//@HEADER
// ************************************************************************
//
//               MeshingGenie: Fracture Meshing Services Package
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

// R5.0

#include "MeshingGenie_2d.h"

using namespace std;

MeshingGenie_2d::MeshingGenie_2d(double dm, std::vector<double> &ExternalBoundaries,
                                std::vector< std::vector<double> > &Holes,
                      std::vector< std::vector<double> > &Cracks,
                     int output, bool plot_results)
{
  _debug_code = false; _save_point_cloud = false;

  input_translator(ExternalBoundaries, Cracks);

    _dm = dm; _dm_squared = _dm * _dm;
  _ExternalBoundaries = ExternalBoundaries;
  _Holes = Holes; _Cracks = Cracks;
  _output = output; _plot_results = plot_results;

  _max_num_poly_pnts = 0;

  _num_neighbor_cells = 21;

  _MY_RAND_MAX = RAND_MAX; // 15 bits
    _MY_RAND_MAX = _MY_RAND_MAX << 15;
    _MY_RAND_MAX += RAND_MAX; // it's now 30 bits
}

int MeshingGenie_2d::execute()
{
  #pragma region Execution of Point Randomizer:

  size_t seed((unsigned) time(0));
  if (_fixed_seed > 0) seed = _fixed_seed;

  srand(seed);
  std::cout<< "* Seed of the random number generator = " << seed << std::endl;

  _num_polys = 10;
  _polys = new PolyPoint*[_num_polys];
  _polys_area = new double[_num_polys];
  for (size_t ipoly = 0; ipoly < _num_polys; ipoly++)
  {
    _polys[ipoly] = 0;
    _polys_area[ipoly] = 0.0;
  }

    clock_t start_time, end_time; double cpu_time;
    std::string ss; // string for printed messages

    start_time = clock();

    std::cout << "\n Process [ Random Point Distribution ]" << std::endl;

    apply_grid_based_method();

  if (_debug_code) return 0;

    end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    print_message("Total Execution Time = " + double_to_string(cpu_time) + " seconds.");

  if (_output == 1 || _output == 2)
  {
    std::cout << "\n Process [ Constrained Voronoi Tessellation ]" << std::endl;
    CDT(true);
  }

  if (_plot_results)
  {
    plot_vertices("sprinkled_vertices.ps", false, false);
    plot_vertices("sprinkled_vertices_grid.ps", false, true);
    plot_vertices("sprinkled_vertices_grid_circles.ps", true, true);
    plot_vertices("sprinkled_vertices_circles.ps", true, false);
  }
  //std::cout << "max poly points = " << _max_num_poly_pnts << std::endl;
    return 0;
  #pragma endregion
}

void MeshingGenie_2d::get_point_cloud(std::vector<double>& x, std::vector<double> &y)
{
  #pragma region Retrieve The generated Point Cloud:
  x.clear(); y.clear();
  size_t ii(0), num_points(0); Point* qi; Point* qj;
  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_points[icell] == 0) continue;

    if (ii == 0) qi = _cell_points[icell];

    x.push_back(qi->x); y.push_back(qi->y);

    if (qi->next == 0)
    {
      ii = 0;
      delete qi;
    }
    else
    {
      qj = qi;
      qi = qi->next; icell--; // stay in the same cell
      delete qj;
      ii++;
    }
  }
  delete [] _cell_points;
  #pragma endregion
}

void MeshingGenie_2d::get_CDT_Tessellation(std::vector<double> &x, std::vector<double> &y, std::vector< std::vector<size_t> > &elements)
{
  #pragma region Retrieve The generated Point Cloud:
  x.clear(); y.clear(); elements.clear();
  size_t ii(0), num_points(0); Point* qi; Point* qj; Point* qk;
  std::map<Point*, size_t> points_indices;
  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_points[icell] == 0) continue;

    if (ii == 0) qi = _cell_points[icell];

    points_indices[qi] = num_points; num_points++;

    if (qi->next == 0) ii = 0;
    else
    {
      qi = qi->next; icell--; // stay in the same cell
      ii++;
    }
  }

  x.reserve(num_points);
  y.reserve(num_points);
  elements.reserve(num_points * 3);

  ii = 0;
  PointEdge* edge; size_t qi_index, qj_index, qk_index;
  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_points[icell] == 0) continue;

    if (ii == 0) qi = _cell_points[icell];

    double xx = qi->x;
    double yy = qi->y;

    qi_index = points_indices[qi];
    x.push_back(xx); y.push_back(yy);

    edge = qi->edges;
    while (edge != 0)
    {
      if (edge->next == 0) break; // a closed loop

      qj = edge->edge_other_end;
      qk = edge->next->edge_other_end;

      if (connected_points(qj, qk))
      {
        qj_index = points_indices[qj];
        qk_index = points_indices[qk];

        if (qi_index < qj_index && qi_index < qk_index)
        {
          std::vector<size_t> triangle(3);
          triangle[0] = points_indices[qi];
          triangle[1] = points_indices[qj];
          triangle[2] = points_indices[qk];
          elements.push_back(triangle);
        }
      }

      edge = edge->next;
      if (edge == qi->edges) break; // a closed loop
    }

    if (qi->next == 0) ii = 0;
    else
    {
      qi = qi->next; icell--; // stay in the same cell
      ii++;
    }
  }
  // NEED TO CLEAN MEMORY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #pragma endregion
}

void MeshingGenie_2d::get_Voronoi_Tessellation(std::vector<double> &x, std::vector<double> &y,
                                             std::vector< std::vector<size_t> > &elements, double min_edge_length,
                         bool generate_side_sets, bool split_cracktipcells)
{
  #pragma region Retrieve The generated Point Cloud:

  double TOL = min_edge_length * min_edge_length;
  if (TOL < 1E-10) TOL = 1E-10;

  std::map<Point*, size_t> points_indices;
  std::map<Point*, Point*> point_target;
  std::map<Point*, Point*>::iterator target_iter;

  if (min_edge_length < 0) min_edge_length = 0;
  else min_edge_length *= min_edge_length;

  x.clear(); y.clear(); elements.clear();
  size_t num_points(0); Point* qi; Point* qj; Point* qk;

  bool first(true);
  _cell_nodes = new Point*[_num_cells];
  for (size_t icell = 0; icell < _num_cells; icell++) _cell_nodes[icell] = 0;

  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_points[icell] == 0) continue;

    if (first) qi = _cell_points[icell];

    num_points++;

    if (qi->next != 0)
    {
      qi = qi->next; icell--; // stay in the same cell
      first = false;
    }
    else first = true;
  }

  // Number of the points are almost the same
  elements.reserve(num_points);

  std::list<double> vxy;
  std::list<double>::iterator iter;
  std::vector<bool> vboundary;
  std::vector<bool> vboundary_(3);

  std::vector< std::vector<Point*> > elements_nodes;
  std::vector< std::vector<bool> > boundary_points;
  elements_nodes.reserve(num_points);
  boundary_points.reserve(num_points);

  std::vector<Point*> element;
  std::vector<size_t> element_nodes;

  std::vector<size_t> sidesets_element_index;
  std::vector<size_t> sidesets_local_index;
  std::vector<size_t> sidesets_edge_index;

  // internal edges
  int ked;
  double* xed_1 = new double[100];
  double* yed_1 = new double[100];
  double* xed_2 = new double[100];
  double* yed_2 = new double[100];

  first = true;
  PointEdge* edge; size_t edge_id_i, edge_id_j, edge_id_k;
  size_t edge_index; bool dir;
  PointEdge* edge_j;
  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_points[icell] == 0) continue;
    if (first) qi = _cell_points[icell];

    // collect internal edges from star
    edge = qi->edges; ked = 0;
    while (edge != 0)
    {
      #pragma region Collect edges of the surrounding triangles:
      if (edge->next == 0) break; // a closed loop

      qj = edge->edge_other_end;
      edge_j = qj->edges;
      while (edge_j != 0)
      {
        qk = edge_j->edge_other_end;
        if (boundary_edge_exists(qj, qk))
        {

          xed_1[ked] = qj->x;
          yed_1[ked] = qj->y;
          xed_2[ked] = qk->x;
          yed_2[ked] = qk->y;
          ked++;
        }
        edge_j = edge_j->next;
        if (edge_j ==  qj->edges) break;
      }
      edge = edge->next;
      if (edge == qi->edges) break;
      #pragma endregion
    }

    bool closed_loop(false);
    edge = qi->edges;
    vxy.clear();
    vboundary.clear();

    size_t ivtx(0);
    while (edge != 0)
    {
      #pragma region Collect the circum centers of the surrounding triangles:
      if (edge->next == 0) break; // a closed loop

      qj = edge->edge_other_end;
      qk = edge->next->edge_other_end;

      if (connected_points(qj, qk))
      {
        // get the two circles
        double cx, cy, cr;
        circumcircle(qi->x, qi->y, qj->x, qj->y, qk->x, qk->y, cx, cy, cr);

        bool bi = is_boundary_edge(qi, qj, edge_id_i);
        bool bj = is_boundary_edge(qj, qk, edge_id_j);
        bool bk = is_boundary_edge(qk, qi, edge_id_k);

        bool ci(false), cj(false), ck(false);
        if (bi)
        {
          get_edge_data(edge_id_i, edge_index, dir);
          if (edge_index >= _num_boundary_edges + _num_holes_edges) ci = true;
        }
        if (bj)
        {
          get_edge_data(edge_id_j, edge_index, dir);
          if (edge_index >= _num_boundary_edges + _num_holes_edges) cj = true;
        }
        if (bk)
        {
          get_edge_data(edge_id_k, edge_index, dir);
          if (edge_index >= _num_boundary_edges + _num_holes_edges) ck = true;
        }

        if (ked > 0)
        {
          double Li = distance_squared(qi->x - qj->x, qi->y - qj->y);
          double Lj = distance_squared(qj->x - qk->x, qj->y - qk->y);
          double Lk = distance_squared(qk->x - qi->x, qj->y - qk->y);

          bool skip(false);
          if (Li > Lj + Lk && !bi)
          {
            for (int kk = 0; kk < ked; kk++)
            {
              if (crossing_segments(xed_1[kk], yed_1[kk], xed_2[kk], yed_2[kk], cx, cy, qk->x, qk->y))
              {
                skip = true;
                break;
              }
            }
            if (skip)
            {
              if (bk)
              {
                vxy.push_back(0.5 * (qi->x + qk->x)); vxy.push_back(0.5 * (qi->y + qk->y));    vboundary.push_back(true);  ivtx++;
              }
              edge = edge->next;
              if (edge == qi->edges) break; // a closed loop
              continue;
            }
          }
          if (Lj > Lk + Li && !bj)
          {
            double xmid(0.0), ymid(0.0);
            for (int kk = 0; kk < ked; kk++)
            {
              if (crossing_segments(xed_1[kk], yed_1[kk], xed_2[kk], yed_2[kk], cx, cy, qi->x, qi->y))
              {
                xmid = 0.5 * (xed_1[kk] + xed_2[kk]);
                ymid = 0.5 * (yed_1[kk] + yed_2[kk]);
                skip = true;
                break;
              }
            }
            if (skip)
            {
              if (bi)
              {
                vxy.push_back(qi->x); vxy.push_back(qi->y);                                    vboundary.push_back(true);  ivtx++;
                vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));    vboundary.push_back(true);  ivtx++;
              }
              vxy.push_back(xmid); vxy.push_back(ymid);                                          vboundary.push_back(true);  ivtx++;
              if (bk)
              {
                vxy.push_back(0.5 * (qi->x + qk->x)); vxy.push_back(0.5 * (qi->y + qk->y));    vboundary.push_back(true);  ivtx++;
              }
              edge = edge->next;
              if (edge == qi->edges) break; // a closed loop
              continue;
            }
          }
          if (Lk > Li + Lj && !bk)
          {
            for (int kk = 0; kk < ked; kk++)
            {
              if (crossing_segments(xed_1[kk], yed_1[kk], xed_2[kk], yed_2[kk], cx, cy, qj->x, qj->y))
              {
                skip = true;
                break;
              }
            }
            if (skip)
            {
              if (bi)
              {
                vxy.push_back(qi->x); vxy.push_back(qi->y);                                    vboundary.push_back(true);  ivtx++;
                vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));    vboundary.push_back(true);  ivtx++;
              }
              edge = edge->next;
              if (edge == qi->edges) break; // a closed loop
              continue;
            }
          }
        }

        double px, py;

        if (bi && bk)
        {
          #pragma region a face with starting and ending boundary edges:
          bool in_tri = point_in_triangle(cx, cy, qi->x, qi->y, qj->x, qj->y, qk->x, qk->y);
          if (in_tri)
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                                 vboundary.push_back(true); ivtx++;
            vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y)); vboundary.push_back(true); ivtx++;
            vxy.push_back(cx); vxy.push_back(cy);                                       vboundary.push_back(false); ivtx++;

            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(0.5 * (qi->x + qk->x)); vxy.push_back(0.5 * (qi->y + qk->y)); vboundary.push_back(true); ivtx++;

            create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
            elements_nodes.push_back(element);
            boundary_points.push_back(vboundary); vboundary.clear();
          }
          else if (bk && crossing_segments(0.5 * (qi->x + qj->x), 0.5 * (qi->y + qj->y), cx, cy, qi->x, qi->y, qk->x, qk->y, px, py))
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                                   vboundary.push_back(true); ivtx++;
            vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));   vboundary.push_back(true); ivtx++;

            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            if (ck)
            {
              vxy.push_back(0.5 * (qi->x + qk->x)); vxy.push_back(0.5 * (qi->y + qk->y));   vboundary.push_back(true); ivtx++;
            }
            else
            {
              vxy.push_back(px); vxy.push_back(py);                             vboundary.push_back(true); ivtx++;
            }

            create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
            elements_nodes.push_back(element);
            boundary_points.push_back(vboundary); vboundary.clear();
          }
          else if (bi && crossing_segments(0.5 * (qi->x + qk->x), 0.5 * (qi->y + qk->y), cx, cy, qi->x, qi->y, qj->x, qj->y, px, py))
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                      vboundary.push_back(true); ivtx++;
            if (ci)
            {
              vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));    vboundary.push_back(true); ivtx++;
            }
            else
            {
              vxy.push_back(px); vxy.push_back(py);                                          vboundary.push_back(true); ivtx++;
            }
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(0.5 * (qi->x + qk->x)); vxy.push_back(0.5 * (qi->y + qk->y));    vboundary.push_back(true); ivtx++;
            create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
            elements_nodes.push_back(element);
            boundary_points.push_back(vboundary); vboundary.clear();
          }
          else if (bj && crossing_segments(0.5 * (qi->x + qj->x), 0.5 * (qi->y + qj->y), cx, cy, qj->x, qj->y, qk->x, qk->y, px, py))
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                                          vboundary.push_back(true); ivtx++;
            vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));          vboundary.push_back(true); ivtx++;
            if (cj)
            {
              vxy.push_back(0.5 * (qj->x + qk->x)); vxy.push_back(0.5 * (qj->y + qk->y));      vboundary.push_back(true); ivtx++;
            }
            else
            {
              if (generate_side_sets)
              {
                sidesets_element_index.push_back(elements_nodes.size());
                sidesets_local_index.push_back(ivtx);
                get_edge_data(edge_id_j, edge_index, dir);
                sidesets_edge_index.push_back(edge_index);
              }
              vxy.push_back(px); vxy.push_back(py);                                            vboundary.push_back(true); ivtx++;
              crossing_segments(0.5 * (qi->x + qk->x), 0.5 * (qi->y + qk->y), cx, cy, qj->x, qj->y, qk->x, qk->y, px, py);
              vxy.push_back(px); vxy.push_back(py);                                            vboundary.push_back(true); ivtx++;
            }
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(0.5 * (qi->x + qk->x)); vxy.push_back(0.5 * (qi->y + qk->y));          vboundary.push_back(true); ivtx++;
            create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
            elements_nodes.push_back(element);
            boundary_points.push_back(vboundary); vboundary.clear();
          }
          else
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                                      vboundary.push_back(true); ivtx++;
            vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));      vboundary.push_back(true); ivtx++;
            vxy.push_back(cx); vxy.push_back(cy);                                  vboundary.push_back(false); ivtx++;
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(0.5 * (qi->x + qk->x)); vxy.push_back(0.5 * (qi->y + qk->y));      vboundary.push_back(true); ivtx++;
            create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
            elements_nodes.push_back(element);
            boundary_points.push_back(vboundary); vboundary.clear();
          }
          #pragma endregion
        }
        else if (bi)
        {
          #pragma region a face with starting boundary edge:
          bool in_tri = point_in_triangle(cx, cy, qi->x, qi->y, qj->x, qj->y, qk->x, qk->y);
          if (in_tri)
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                                    vboundary.push_back(true);  ivtx++;
            vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));    vboundary.push_back(true);  ivtx++;
            vxy.push_back(cx); vxy.push_back(cy);                              vboundary.push_back(false); ivtx++;
          }
          else if (crossing_segments(0.5 * (qi->x + qk->x), 0.5 * (qi->y + qk->y), cx, cy, qi->x, qi->y, qj->x, qj->y, px, py))
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                                        vboundary.push_back(true); ivtx++;
            if (ci)
            {
              vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));    vboundary.push_back(true); ivtx++;
            }
            else
            {
              vxy.push_back(px); vxy.push_back(py);                              vboundary.push_back(true); ivtx++;
            }
          }
          else if (bj && crossing_segments(0.5 * (qi->x + qj->x), 0.5 * (qi->y + qj->y), cx, cy, qj->x, qj->y, qk->x, qk->y, px, py))
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                                        vboundary.push_back(true); ivtx++;
            vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));        vboundary.push_back(true); ivtx++;
            if (cj)
            {
              vxy.push_back(0.5 * (qj->x + qk->x)); vxy.push_back(0.5 * (qj->y + qk->y));    vboundary.push_back(true); ivtx++;
            }
            else
            {
              if (generate_side_sets)
              {
                sidesets_element_index.push_back(elements_nodes.size());
                sidesets_local_index.push_back(ivtx);
                get_edge_data(edge_id_j, edge_index, dir);
                sidesets_edge_index.push_back(edge_index);
              }
              vxy.push_back(px); vxy.push_back(py);                                          vboundary.push_back(true); ivtx++;
              crossing_segments(0.5 * (qi->x + qk->x), 0.5 * (qi->y + qk->y), cx, cy, qj->x, qj->y, qk->x, qk->y, px, py);
              vxy.push_back(px); vxy.push_back(py);                              vboundary.push_back(true); ivtx++;
            }
          }
          else
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_i, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(qi->x); vxy.push_back(qi->y);                                    vboundary.push_back(true); ivtx++;
            vxy.push_back(0.5 * (qi->x + qj->x)); vxy.push_back(0.5 * (qi->y + qj->y));    vboundary.push_back(true); ivtx++;
            vxy.push_back(cx); vxy.push_back(cy);                                vboundary.push_back(false); ivtx++;
          }
          #pragma endregion
        }
        else if (bk)
        {
          #pragma region a face with an ending boundary edge:
          bool in_tri = point_in_triangle(cx, cy, qi->x, qi->y, qj->x, qj->y, qk->x, qk->y);
          if (in_tri)
          {
            vxy.push_back(cx); vxy.push_back(cy);                                            vboundary.push_back(false); ivtx++;
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            double xen(0.5 * (qi->x + qk->x)), yen(0.5 * (qi->y + qk->y));
            vxy.push_back(0.5 * (qi->x + qk->x)); vxy.push_back(0.5 * (qi->y + qk->y));      vboundary.push_back(true);  ivtx++;

            // check if qi is an end of a crack tip
            iter = vxy.begin();
            double xo = *iter; iter++;
            double yo = *iter; iter++;
            double xe1 = *iter; iter++;
            double ye1 = *iter; iter++;
            bool split_cell(false);

            if (split_cracktipcells &&
              (invalid_triangle(xo, yo, xe1, ye1, xen, yen, 1E-10) ||
               ZeroAngle(xo, yo, xe1, ye1, xen, yen, 1E-10))) split_cell = true;

            if (split_cell)
            {
              Point* qq;
              Point* qm;
              Point* qp;
              iter = vxy.begin();
              double xx = *iter; iter++;
              double yy = *iter; iter++;
              size_t jj(0);
              get_closest_point(xx, yy, vboundary[jj], qq, TOL);
              vboundary_[0] = vboundary[jj]; jj++;
              xx = *iter; iter++;
              yy = *iter; iter++;
              get_closest_point(xx, yy, vboundary[jj], qp, TOL);
              vboundary_[1] = vboundary[jj]; jj++;
              while (iter != vxy.end())
              {
                xx = *iter; iter++;
                yy = *iter; iter++;
                get_closest_point(xx, yy, vboundary[jj], qm, TOL);
                if (qm != qp)
                {
                  vboundary_[2] = vboundary[jj]; jj++;
                  element.clear();
                  element.push_back(qq);
                  element.push_back(qp);
                  element.push_back(qm);
                  elements_nodes.push_back(element);
                  boundary_points.push_back(vboundary_);
                  qp = qm;
                  vboundary_[1] = vboundary_[2];
                }
              }
              element.clear();
              vboundary.clear();
              vxy.clear(); ivtx = 0;
            }
            else
            {
              create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
              elements_nodes.push_back(element);
              boundary_points.push_back(vboundary); vboundary.clear();
            }
          }
          else if (crossing_segments(0.5 * (qi->x + qj->x), 0.5 * (qi->y + qj->y), cx, cy, qi->x, qi->y, qk->x, qk->y, px, py))
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            double xen(0.5 * (qi->x + qk->x)), yen(0.5 * (qi->y + qk->y));
            if (ck)
            {
              vxy.push_back(xen); vxy.push_back(yen);   vboundary.push_back(true); ivtx++;
            }
            else
            {
              xen = px; yen = py;
              vxy.push_back(px); vxy.push_back(py);                             vboundary.push_back(true); ivtx++;
            }

            // check if qi is an end of a crack tip
            iter = vxy.begin();
            double xo = *iter; iter++;
            double yo = *iter; iter++;
            double xe1 = *iter; iter++;
            double ye1 = *iter; iter++;
            bool split_cell(false);


            if (split_cracktipcells &&
              (invalid_triangle(xo, yo, xe1, ye1, xen, yen, 1E-10) ||
               ZeroAngle(xo, yo, xe1, ye1, xen, yen, 1E-10))) split_cell = true;


            if (split_cell)
            {
              Point* qq;
              Point* qm;
              Point* qp;
              iter = vxy.begin();
              double xx = *iter; iter++;
              double yy = *iter; iter++;
              size_t jj(0);
              get_closest_point(xx, yy, vboundary[jj], qq, TOL);
              vboundary_[0] = vboundary[jj]; jj++;
              xx = *iter; iter++;
              yy = *iter; iter++;
              get_closest_point(xx, yy, vboundary[jj], qp, TOL);
              vboundary_[1] = vboundary[jj]; jj++;
              while (iter != vxy.end())
              {
                xx = *iter; iter++;
                yy = *iter; iter++;
                get_closest_point(xx, yy, vboundary[jj], qm, TOL);
                if (qm != qp)
                {
                  vboundary_[2] = vboundary[jj]; jj++;
                  element.clear();
                  element.push_back(qq);
                  element.push_back(qp);
                  element.push_back(qm);
                  elements_nodes.push_back(element);
                  boundary_points.push_back(vboundary_);
                  qp = qm;
                  vboundary_[1] = vboundary_[2];
                }
              }
              element.clear();
              vboundary.clear();
              vxy.clear(); ivtx = 0;
            }
            else
            {
              create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
              elements_nodes.push_back(element);
              boundary_points.push_back(vboundary); vboundary.clear();
            }
          }
          else if (bj && crossing_segments(0.5 * (qi->x + qj->x), 0.5 * (qi->y + qj->y), cx, cy, qj->x, qj->y, qk->x, qk->y, px, py))
          {
            if (cj)
            {
              vxy.push_back(0.5 * (qj->x + qk->x)); vxy.push_back(0.5 * (qj->y + qk->y));    vboundary.push_back(true); ivtx++;
            }
            else
            {
              if (generate_side_sets)
              {
                sidesets_element_index.push_back(elements_nodes.size());
                sidesets_local_index.push_back(ivtx);
                get_edge_data(edge_id_j, edge_index, dir);
                sidesets_edge_index.push_back(edge_index);
              }
              vxy.push_back(px); vxy.push_back(py);                                             vboundary.push_back(true); ivtx++;
              crossing_segments(0.5 * (qi->x + qk->x), 0.5 * (qi->y + qk->y), cx, cy, qj->x, qj->y, qk->x, qk->y, px, py);
              vxy.push_back(px); vxy.push_back(py);                                             vboundary.push_back(true); ivtx++;
            }
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            double xen(0.5 * (qi->x + qk->x)), yen(0.5 * (qi->y + qk->y));
            vxy.push_back(xen); vxy.push_back(yen);       vboundary.push_back(true); ivtx++;

            // check if qi is an end of a crack tip
            iter = vxy.begin();
            double xo = *iter; iter++;
            double yo = *iter; iter++;
            double xe1 = *iter; iter++;
            double ye1 = *iter; iter++;
            bool split_cell(false);

            if (split_cracktipcells &&
              (invalid_triangle(xo, yo, xe1, ye1, xen, yen, 1E-10) ||
               ZeroAngle(xo, yo, xe1, ye1, xen, yen, 1E-10))) split_cell = true;


            if (split_cell)
            {
              Point* qq;
              Point* qm;
              Point* qp;
              iter = vxy.begin();
              double xx = *iter; iter++;
              double yy = *iter; iter++;
              size_t jj(0);
              get_closest_point(xx, yy, vboundary[jj], qq, TOL);
              vboundary_[0] = vboundary[jj]; jj++;
              xx = *iter; iter++;
              yy = *iter; iter++;
              get_closest_point(xx, yy, vboundary[jj], qp, TOL);
              vboundary_[1] = vboundary[jj]; jj++;
              while (iter != vxy.end())
              {
                xx = *iter; iter++;
                yy = *iter; iter++;
                get_closest_point(xx, yy, vboundary[jj], qm, TOL);
                if (qm != qp)
                {
                  vboundary_[2] = vboundary[jj]; jj++;
                  element.clear();
                  element.push_back(qq);
                  element.push_back(qp);
                  element.push_back(qm);
                  elements_nodes.push_back(element);
                  boundary_points.push_back(vboundary_);
                  qp = qm;
                  vboundary_[1] = vboundary_[2];
                }
              }
              element.clear();
              vboundary.clear();
              vxy.clear(); ivtx = 0;
            }
            else
            {
              create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
              elements_nodes.push_back(element);
              boundary_points.push_back(vboundary); vboundary.clear();
            }
          }
          else
          {
            vxy.push_back(cx); vxy.push_back(cy);                                             vboundary.push_back(false); ivtx++;
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_k, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            double xen(0.5 * (qi->x + qk->x)), yen(0.5 * (qi->y + qk->y));
            vxy.push_back(xen); vxy.push_back(yen);       vboundary.push_back(true);  ivtx++;

            // check if qi is an end of a crack tip
            iter = vxy.begin();
            double xo = *iter; iter++;
            double yo = *iter; iter++;
            double xe1 = *iter; iter++;
            double ye1 = *iter; iter++;
            bool split_cell(false);

            if (split_cracktipcells &&
              (invalid_triangle(xo, yo, xe1, ye1, xen, yen, 1E-10) ||
               ZeroAngle(xo, yo, xe1, ye1, xen, yen, 1E-10))) split_cell = true;


            if (split_cell)
            {
              Point* qq;
              Point* qm;
              Point* qp;
              iter = vxy.begin();
              double xx = *iter; iter++;
              double yy = *iter; iter++;
              size_t jj(0);
              get_closest_point(xx, yy, vboundary[jj], qq, TOL);
              vboundary_[0] = vboundary[jj]; jj++;
              xx = *iter; iter++;
              yy = *iter; iter++;
              get_closest_point(xx, yy, vboundary[jj], qp, TOL);
              vboundary_[1] = vboundary[jj]; jj++;
              while (iter != vxy.end())
              {
                xx = *iter; iter++;
                yy = *iter; iter++;
                get_closest_point(xx, yy, vboundary[jj], qm, TOL);
                if (qm != qp)
                {
                  vboundary_[2] = vboundary[jj]; jj++;
                  element.clear();
                  element.push_back(qq);
                  element.push_back(qp);
                  element.push_back(qm);
                  elements_nodes.push_back(element);
                  boundary_points.push_back(vboundary_);
                  qp = qm;
                  vboundary_[1] = vboundary_[2];
                }
              }
              element.clear();
              vboundary.clear();
              vxy.clear(); ivtx = 0;
            }
            else
            {
              create_element(vxy, iter, vboundary, element, TOL); ivtx = 0;
              elements_nodes.push_back(element);
              boundary_points.push_back(vboundary); vboundary.clear();
            }
          }
          #pragma endregion
        }
        else if (bj && crossing_segments(0.5 * (qi->x + qj->x), 0.5 * (qi->y + qj->y), cx, cy, qj->x, qj->y, qk->x, qk->y, px, py))
        {
          if (cj)
          {
            vxy.push_back(0.5 * (qj->x + qk->x)); vxy.push_back(0.5 * (qj->y + qk->y));    vboundary.push_back(true); ivtx++;
          }
          else
          {
            if (generate_side_sets)
            {
              sidesets_element_index.push_back(elements_nodes.size());
              sidesets_local_index.push_back(ivtx);
              get_edge_data(edge_id_j, edge_index, dir);
              sidesets_edge_index.push_back(edge_index);
            }
            vxy.push_back(px); vxy.push_back(py);                      vboundary.push_back(true); ivtx++;
            crossing_segments(0.5 * (qi->x + qk->x), 0.5 * (qi->y + qk->y), cx, cy, qj->x, qj->y, qk->x, qk->y, px, py);
            vxy.push_back(px); vxy.push_back(py);            vboundary.push_back(true); ivtx++;
          }
        }
        else
        {
          vxy.push_back(cx); vxy.push_back(cy); vboundary.push_back(false); ivtx++;
        }
      }
      edge = edge->next;
      if (edge == qi->edges) break; // a closed loop
      #pragma endregion
    }

    if (vxy.size() > 0)
    {
      create_element(vxy, iter, vboundary, element, TOL);
      elements_nodes.push_back(element);
      boundary_points.push_back(vboundary); vboundary.clear();
    }

    if (qi->next != 0)
    {
      qi = qi->next; icell--; // stay in the same cell
      first = false;
    }
    else
    {
      first = true;
    }
  }

  size_t num_elements(elements_nodes.size());
  first = true; size_t node_index(0);
  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_nodes[icell] == 0) continue;
    if (first) qi = _cell_nodes[icell];

    points_indices[qi] = node_index; node_index++;

    if (qi->next != 0)
    {
      qi = qi->next; icell--; // stay in the same cell
      first = false;
    }
    else
    {
      first = true;
    }
  }

  if (generate_side_sets)
  {
    _side_sets.resize(_edges_x1.size());
  }

  size_t iss(0);
  for (size_t ie = 0; ie < num_elements; ie++)
  {
    #pragma region Storing elements using nodes indices instead of Points:
    size_t num_element_nodes(elements_nodes[ie].size());
    element_nodes.clear(); element_nodes.reserve(num_element_nodes);
    Point* pold(0); size_t iimax = num_element_nodes - 1;
    size_t num_collapsed(0);
    for (size_t ii = 0; ii < num_element_nodes; ii++)
    {
      if (elements_nodes[ie][ii] == pold)
      {
        num_collapsed++;
        if (generate_side_sets && iss < sidesets_element_index.size() && sidesets_element_index[iss] == ie &&
                                                                         sidesets_local_index[iss] == ii)
        {
          edge_index = sidesets_edge_index[iss];
          _side_sets[edge_index].push_back(ie);
          _side_sets[edge_index].push_back(ii - num_collapsed); iss++;
        }
        continue;
      }
      if (ii == iimax && elements_nodes[ie][ii] == elements_nodes[ie][0])
      {
        num_collapsed++;
        if (generate_side_sets && iss < sidesets_element_index.size() && sidesets_element_index[iss] == ie &&
                                                                         sidesets_local_index[iss] == ii)
        {
          edge_index = sidesets_edge_index[iss];
          _side_sets[edge_index].push_back(ie);
          _side_sets[edge_index].push_back(ii - num_collapsed); iss++;
        }
        break;
      }

      if (generate_side_sets && iss < sidesets_element_index.size() && sidesets_element_index[iss] == ie &&
        sidesets_local_index[iss] == ii)
      {
        edge_index = sidesets_edge_index[iss];
        _side_sets[edge_index].push_back(ie);
        _side_sets[edge_index].push_back(ii - num_collapsed); iss++;
      }

      element_nodes.push_back(points_indices[elements_nodes[ie][ii]]);
      pold = elements_nodes[ie][ii];
    }
    elements.push_back(element_nodes);
    #pragma endregion
  }

  x.resize(node_index); y.resize(node_index);
  std::map<Point*, size_t>::iterator  piter;
  for (piter = points_indices.begin(); piter!= points_indices.end(); piter++)
  {
    Point* q = piter->first;
    size_t index = piter->second;
    x[index] = q->x; y[index] = q->y;
  }

  // clear memory
  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_nodes[icell] == 0) continue;
    delete _cell_nodes[icell];
  }
  delete[] _cell_nodes;
  #pragma endregion
}


int MeshingGenie_2d::apply_grid_based_method()
{
    #pragma region Grid Based Method:

    clock_t start_time, end_time; double cpu_time;

    print_message("Extracting Boundary Edges Data");
    start_time = clock();
    extract_boudary_edges();
    end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    print_message_end("^ Execution Time = " + double_to_string(cpu_time) + " seconds.");

    print_message("Generating Background Grid ");
    start_time = clock();
    generate_backgorund_grid(0.70710678 * _dm);
    end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    print_message_end("^ Execution Time = " + double_to_string(cpu_time) + " seconds.");
    print_message_end("^ Number of utilized cells = " + double_to_string(1.0 * _num_cells));

    _num_bad_cells = 0; _num_sprinkled = 0;

    print_message("Marking Exterior Cells");
    start_time = clock();
    mark_cells_crossing_boundaries( );
  for (size_t i = 0; i < _num_cells; i++)
  {
    if (_bad_cells[i]) _blue_cells[i] = true;
  }

  _cell_points = new Point*[_num_cells];
  for (size_t ii = 0; ii < _num_cells; ii++)
  {
    _cell_points[ii] = 0;
  }

  fill_neighbors_data();

  for (size_t i = 0; i < _num_cells; i++)
  {
    _bad_cells[i] = false;
  }
  _bad_cells[0] = true; _num_bad_cells = 1;

  //plot_cells("blue_cells.ps", true);
  //plot_cells("valid_cells.ps", false);

  while (true)
  {
    bool nothing_new(true); int jcell;
    for (size_t icell = 0; icell < _num_cells; icell++)
    {
      if (_bad_cells[icell])
      {
        // four neighbors are external if not blue "boundary"
        for (size_t j = 1; j < 5; j++)
        {
          jcell = int(icell) + _neighbors[j];
          if (jcell >= 0 && jcell < int(_num_cells) && !_blue_cells[jcell] && !_bad_cells[jcell])
          {
            _bad_cells[jcell] = true; _num_bad_cells++;
            nothing_new = false;
          }
        }
      }
    }
    if (nothing_new) break;
  }

  //plot_cells("blue_cells.ps", true);
  //plot_cells("valid_cells.ps", false);

    //invalidate_external_cells();
    end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    print_message_end("^ Execution Time = " + double_to_string(cpu_time) + " seconds.");


  if (true)
  {
    print_message("Sprinkling points along domain boundaries");
    start_time = clock();
    //sprinkle_points_along_boundaries(0.5 * sqrt(3.0) * _dm);
    sprinkle_points_along_boundaries(_dm / sqrt(2.0));
    end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    print_message_end("^ Execution Time = " + double_to_string(cpu_time) + " seconds.");
    //plot_vertices("sprinkled_vertices_circles.ps", true, false);
  }

    print_message("Sprinkling points in the domain");
  apply_maximal_poisson_disk_sampler();
  //plot_vertices("sprinkled_vertices_circles.ps", true, false);
    return 0;
    #pragma endregion
}

int MeshingGenie_2d::extract_boudary_edges()
{
    #pragma region Extract the boundary edges data:
    _num_boundary_edges = 0;
    _num_holes_edges = 0;
    std::vector<double> polyline(_ExternalBoundaries);
    size_t num_points(polyline.size() - 1);
  _xmin = polyline[0]; _ymin = polyline[1];
  _xmax = polyline[0]; _ymax = polyline[1];
    for (size_t i = 0; i < num_points; i+=2)
    {
        size_t ip(i + 1);
        size_t j(i + 2), jp(i + 3);
        if (j >= num_points){j = 0; jp = 1;}

        double xx(polyline[i])  , yy(polyline[ip]);
        double xp(polyline[j])  , yp(polyline[jp]);

    if (xx < _xmin) _xmin = xx;
    if (xx > _xmax) _xmax = xx;
    if (yy < _ymin) _ymin = yy;
    if (yy > _ymax) _ymax = yy;
        _edges_x1.push_back(xx); _edges_y1.push_back(yy);
        _edges_x2.push_back(xp); _edges_y2.push_back(yp);
        _num_boundary_edges++;
    }

    // Add Internal holes to the edges list
  size_t num_holes(_Holes.size());
    for (size_t ii = 0; ii < num_holes; ii++)
    {
        polyline.clear();
        polyline = _Holes[ii];
        size_t num_points(polyline.size() - 1);
        for (size_t i = 0; i < num_points; i+=2)
        {
            size_t ip(i + 1);
            size_t j(i + 2), jp(i + 3);
            if (j >= num_points){j = 0; jp = 1;}

            double xx(polyline[i])  , yy(polyline[ip]);
            double xp(polyline[j])  , yp(polyline[jp]);

            _edges_x1.push_back(xx); _edges_y1.push_back(yy);
            _edges_x2.push_back(xp); _edges_y2.push_back(yp);
            _num_holes_edges++;
        }
    }

    // Add Internal holes to the edges list
  size_t num_cracks(_Cracks.size());
    for (size_t ii = 0; ii < num_cracks; ii++)
    {
        polyline.clear();
        polyline = _Cracks[ii];
        size_t num_points(polyline.size() - 3);
        for (size_t i = 0; i < num_points; i+=2)
        {
            size_t ip(i + 1);
            size_t j(i + 2), jp(i + 3);

            double xx(polyline[i])  , yy(polyline[ip]);
            double xp(polyline[j])  , yp(polyline[jp]);

            _edges_x1.push_back(xx);
            _edges_y1.push_back(yy);
            _edges_x2.push_back(xp);
            _edges_y2.push_back(yp);
        }
    }
    return 0;
    #pragma endregion
}

int MeshingGenie_2d::generate_backgorund_grid(double ds)
{
    #pragma region Generate Backgorund grid:
    _bg_s = ds;

  double xmin(_xmin), ymin(_ymin), xmax(_xmax), ymax(_ymax);

    double Lx(xmax - xmin); double Ly(ymax - ymin);
    // scale up the boundaing box
    xmin -= 4 * _bg_s; ymin -= 4 * _bg_s; // 4 extra cells
    Lx += 8 * _bg_s; Ly += 8 * _bg_s;

    _nc = int(ceil(Lx / ds)) + 1;
    _nr = int(ceil(Ly / ds)) + 1;

    _imax = _nc - 2;
    _jmax = _nr - 2;

    _bg_xo = xmin; _bg_yo = ymin;

    _num_cells = (_nr - 1) * (_nc - 1);
    _bad_cells = new bool[_num_cells];
  _blue_cells = new bool[_num_cells];

    for (size_t i = 0; i < _num_cells; i++)
  {
    _bad_cells[i] = false;
    _blue_cells[i] = false;
  }
    return 0;
    #pragma endregion
}

int MeshingGenie_2d::mark_cells_crossing_boundaries( )
{
    #pragma region Mark Cells crossing Boundaries:
  bool N, W, S, E;
    size_t num_edges(_num_boundary_edges + _num_holes_edges);
    for (size_t ied = 0; ied < num_edges; ied++)
    {
        double xo(_edges_x1[ied]), yo(_edges_y1[ied]);
        double xn(_edges_x2[ied]), yn(_edges_y2[ied]);
    N = false; W = false; S = false; E = false;

        // starting cell
        size_t io(int(floor((xo - _bg_xo) / _bg_s)));
        size_t jo(int(floor((yo - _bg_yo) / _bg_s)));

        // ending cell
        size_t in(int(floor((xn - _bg_xo) / _bg_s)));
        size_t jn(int(floor((yn - _bg_yo) / _bg_s)));

        int di(int(in - io)), dj(int(jn - jo));

        if (abs(di) > abs(dj))
        {
            bool reversed(false);
            if (di < 0)
            {
                #pragma region Start from left to right:
                double tmpd(xo); xo = xn; xn = tmpd;
                tmpd = yo; yo = yn; yn = tmpd;
                size_t tmp(io); io = in; in = tmp;
                tmp = jo; jo = jn; jn = tmp;
                di = -di; dj = -dj;
                reversed = true;
                #pragma endregion
            }

      if (reversed) N = true; // invalidate up
      else S = true; // invaldate down

            if (dj == 0)
            {
                #pragma region horizontal line:

                bool left(false), right(false), bottom(false), top(false);

                double xi = _bg_xo + io * _bg_s;
                if (fabs(xi - xo) < 1E-10) {io--; left = true;} // edge touches the left line
                if (fabs(xi - xo + _bg_s) < 1E-10) {left = true;} // edge touches the left line
                xi += _bg_s;
                if (fabs(xi - xn) < 1E-10) {in++; right = true;}   // edge touches the upper line
                if (fabs(xi - xn - _bg_s) < 1E-10) {right = true;} // edges touches the upper line

                double yj = _bg_yo + jo * _bg_s;
                if (fabs(yj - yo) < 1E-10) bottom = true; // edge touches the bottom line of the cell
                yj += _bg_s;
                if (fabs(yj - yo) < 1E-10) top = true; // edge touches the top line of the cell

                size_t jbottom(jo), jtop(jo);
                if (bottom && jbottom > 0) jbottom--;
                if (top && jtop < _nr) jtop++;

                for (size_t i = io; i <= in; i++)
                {
                    for (size_t j = jbottom; j <= jtop; j++)
                    {
                        invalidate_cell(i, j, N, W, S, E, ied);
                    }
                }
                #pragma endregion
            }
            else
            {
                #pragma region Slope less than one:

                double xi, yj;

                xi = _bg_xo + (in + 1) * _bg_s;
                if (fabs(xi - xn) < 1E-10) {in++; di++;}

                if (dj > 0)
                {
                    yj = _bg_yo + (jn + 1) * _bg_s;
                    if (fabs(yj - yn) < 1E-10) {jn++; dj++;}
                }
                else
                {
                    yj = _bg_yo + jn * _bg_s;
                    if (fabs(yj - yn) < 1E-10) {jn--; dj--;}
                }

                if (dj > 0)
                {
                    #pragma region line going up:
                    size_t i(io), j(jo);
                    double slope = (yn - yo) / (xn - xo);

                    #pragma region Working On First Cell:
                    invalidate_cell(i, j, N, W, S, E, ied);

                    // lower right corner of cell
                    xi = _bg_xo + i * _bg_s;
                    yj = _bg_yo + j * _bg_s;

                    if (fabs(xi - xo) < 1E-10 && fabs(yj - yo) < 1E-10)
                    {
                        // lower right corner touches edge start
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                        invalidate_cell(i - 1, j - 1, N, W, S, E, ied);
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                        i++; // move right
                    }
                    else if (fabs(xi - xo + _bg_s) < 1E-10 && fabs(yj - yo) < 1E-10)
                    {
                        // lower left corner touches edge start
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                        invalidate_cell(i + 1, j - 1, N, W, S, E, ied);
                        i++; // move right
                    }
                    else if (fabs(xi - xo + _bg_s) < 1E-10 && fabs(yj - yo + _bg_s) < 1E-10)
                    {
                        // upper right corner touches edge start
                        invalidate_cell(i, j + 1, N, W, S, E, ied);
                        invalidate_cell(i + 1, j, N, W, S, E, ied);
                        i++; j++; // move up and right
                    }
                    else if (fabs(xi - xo) < 1E-10 && fabs(yj - yo + _bg_s) < 1E-10)
                    {
                        // upper left corner touches edge start
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                        invalidate_cell(i - 1, j + 1, N, W, S, E, ied);
                        j++;
                    }
                    else if (fabs(xi - xo) < 1E-10)
                    {
                        // start of that sdge lies on the left border of that cell
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                    }
                    #pragma endregion

                    while (true)
                    {
                        // move from io, jo to in, jn to the right of the line
                        invalidate_cell(i, j, N, W, S, E, ied);

                        if (i >= in && j >= jn) break;

                        // lower right corner of cell
                        xi = _bg_xo + (i + 1) * _bg_s;
                        yj = _bg_yo + (j + 1) * _bg_s;

                        if (fabs((yj - yo) / (xi - xo)) > fabs(slope) - 1E-10)
                        {
                            if (fabs(fabs((yj - yo) / (xi - xo)) - fabs(slope)) < 1E-10)
                            {
                                // edge passes through the upper right corner
                                invalidate_cell(i + 1, j, N, W, S, E, ied);
                                invalidate_cell(i, j + 1, N, W, S, E, ied);
                                j++; // move up
                            }
                            i++; // move right
                        }
                        else j++; // move up

                    }
                    #pragma endregion
                }
                else
                {
                    #pragma region line going down:
                    size_t i(io), j(jo);
                    double slope = (yn - yo) / (xn - xo); // negative slope

                    #pragma region Working On First Cell:

                    invalidate_cell(i, j, N, W, S, E, ied);

                    // lower right corner of cell
                    xi = _bg_xo + i * _bg_s;
                    yj = _bg_yo + j * _bg_s;

                    if (fabs(xi - xo) < 1E-10 && fabs(yj - yo) < 1E-10)
                    {
                        // lower right corner touches edge start
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                        invalidate_cell(i - 1, j - 1, N, W, S, E, ied);
                        j--; // move down
                    }
                    else if (fabs(xi - xo + _bg_s) < 1E-10 && fabs(yj - yo) < 1E-10)
                    {
                        // lower left corner touches edge start
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                        invalidate_cell(i + 1, j, N, W, S, E, ied);
                        i++; j--;
                    }
                    else if (fabs(xi - xo + _bg_s) < 1E-10 && fabs(yj - yo + _bg_s) < 1E-10)
                    {
                        // upper right corner touches edge start
                        invalidate_cell(i, j + 1, N, W, S, E, ied);
                        invalidate_cell(i + 1, j + 1, N, W, S, E, ied);
                        i++;
                    }
                    else if (fabs(xi - xo) < 1E-10 && fabs(yj - yo + _bg_s) < 1E-10)
                    {
                        // upper left corner touches edge start
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                        invalidate_cell(i - 1, j + 1, N, W, S, E, ied);
                        invalidate_cell(i, j + 1, N, W, S, E, ied);
                        i++;
                    }
                    else if (fabs(xi - xo) < 1E-10)
                    {
                        // start of that sdge lies on the left border of that cell
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                    }
                    #pragma endregion

                    while (true)
                    {
                        // move from io, jo to in, jn to the right of the line
                        invalidate_cell(i, j, N, W, S, E, ied);

                        if (i >= in && j <= jn) break;

                        // lower right corner of cell
                        xi = _bg_xo + (i + 1) * _bg_s;
                        yj = _bg_yo + j * _bg_s;

                        if (fabs((yj - yo) / (xi - xo)) > fabs(slope) - 1E-10)
                        {
                            if (fabs(fabs((yj - yo) / (xi - xo)) - fabs(slope)) < 1E-10)
                            {
                                // edge passes through the upper left corner
                                invalidate_cell(i + 1, j, N, W, S, E, ied);
                                invalidate_cell(i, j - 1, N, W, S, E, ied);
                                j--; // move down
                            }
                            i++; // move right
                        }
                        else
                        {
                            j--; // move down
                        }
                    }
                    #pragma endregion
                }
                #pragma endregion
            }
        }
        else
        {
            bool reversed(false);

            if (dj < 0)
            {
                #pragma region start from bottom to top:
                double tmpd(xo); xo = xn; xn = tmpd;
                tmpd = yo; yo = yn; yn = tmpd;
                size_t tmp(io); io = in; in = tmp;
                tmp = jo; jo = jn; jn = tmp;
                di = -di; dj = -dj;
                reversed = true;
                #pragma endregion
            }


            if (reversed) W = true; // Invalidated Left
            else E = true; // Invalidate Right

            if (di == 0)
            {
                #pragma region vertical line:

                bool left(false), right(false), bottom(false), top(false);

                double yj = _bg_yo + jo * _bg_s;
                if (fabs(yj - yo) < 1E-10) {jo--; bottom = true;} // edge touches the bottom line
                if (fabs(yj - yo + _bg_s) < 1E-10) {bottom = true;} // edge touches the bottom line
                yj = _bg_yo + (jn + 1) * _bg_s;
                if (fabs(yj - yn) < 1E-10) {jn++; top = true;}   // edge touches the upper line
                if (fabs(yj - yn - _bg_s) < 1E-10) {top = true;} // edges touches the upper line

                double xi = _bg_xo + io * _bg_s;
                if (fabs(xi - xo) < 1E-10) left = true; // edge touches the left line of the cell
                xi += _bg_s;
                if (fabs(xi - xo) < 1E-10) right = true; // edge touches the right line of the cell

                size_t ileft(io), iright(io);
                if (left && ileft > 0) ileft--;
                if (right && iright < _nc) iright++;

                for (size_t j = jo; j <= jn; j++)
                {
                    for (size_t i = ileft; i <= iright; i++)
                    {
                        invalidate_cell(i, j, N, W, S, E, ied);
                    }
                }
                #pragma endregion
            }
            else
            {
                #pragma region Slope greater than one:

                double xi, yj;

                if (di > 0)
                {
                    xi = _bg_xo + (in + 1) * _bg_s;
                    if (fabs(xi - xn) < 1E-10) {in++; di++;}
                }
                else
                {
                    xi = _bg_xo + in * _bg_s;
                    if (fabs(xi - xn) < 1E-10) {in--; di--;}
                }

                yj = _bg_yo + (jn + 1) * _bg_s;
                if (fabs(yj - yn) < 1E-10) {jn++; dj++;}

                if (di > 0)
                {
                    #pragma region line going right:
                    size_t i(io), j(jo);
                    double slope = (yn - yo) / (xn - xo);

                    #pragma region Working On First Cell:
                    invalidate_cell(i, j, N, W, S, E, ied);
                    // lower left corner of cell
                    xi = _bg_xo + i * _bg_s;
                    yj = _bg_yo + j * _bg_s;

                    if (fabs(xi - xo) < 1E-10 && fabs(yj - yo) < 1E-10)
                    {
                        // lower left corner touches edge start
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                        invalidate_cell(i - 1, j - 1, N, W, S, E, ied);
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                        if (fabs(fabs((yj - yo + _bg_s) / (xi - xo + _bg_s)) - fabs(slope)) < 1E-10)
                        {
                            invalidate_cell(i + 1, j, N, W, S, E, ied);
                        }
                        j++; // move up
                    }
                    else if (fabs(xi - xo + _bg_s) < 1E-10 && fabs(yj - yo) < 1E-10)
                    {
                        // lower right corner touches edge start
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                        invalidate_cell(i + 1, j - 1, N, W, S, E, ied);
                        i++; // move right
                    }
                    else if (fabs(xi - xo + _bg_s) < 1E-10 && fabs(yj - yo + _bg_s) < 1E-10)
                    {
                        // upper right corner touches edge start
                        invalidate_cell(i, j + 1, N, W, S, E, ied);
                        invalidate_cell(i + 1, j, N, W, S, E, ied);
                        i++; j++; // move up and right
                    }
                    else if (fabs(xi - xo) < 1E-10 && fabs(yj - yo + _bg_s) < 1E-10)
                    {
                        // upper left corner touches edge start
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                        invalidate_cell(i - 1, j + 1, N, W, S, E, ied);
                        j++;
                    }
                    else if (fabs(yj - yo) < 1E-10)
                    {
                        // start of that sdge lies on the left border of that cell
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                    }
                    #pragma endregion

                    while (true)
                    {
                        // move from io, jo to in, jn to the right of the line
                        invalidate_cell(i, j, N, W, S, E, ied);
                        if (i >= in && j >= jn) break;

                        // upper right corner
                        xi = _bg_xo + (i + 1) * _bg_s;
                        yj = _bg_yo + (j + 1) * _bg_s;

                        if (fabs((yj - yo) / (xi - xo)) < fabs(slope) + 1E-10)
                        {
                            if (fabs(fabs((yj - yo) / (xi - xo)) - fabs(slope)) < 1E-10)
                            {
                                // edge passes through the upper right corner
                                invalidate_cell(i + 1, j, N, W, S, E, ied);
                                invalidate_cell(i, j + 1, N, W, S, E, ied);
                                i++; // move right
                            }
                            j++; // move up
                        }
                        else i++; // move right

                    }
                    #pragma endregion
                }
                else
                {
                    #pragma region line going left:
                    size_t i(io), j(jo);
                    double slope = (yn - yo) / (xn - xo); // negative slope

                    #pragma region Working On First Cell:
                    invalidate_cell(i, j, N, W, S, E, ied);

                    xi = _bg_xo + i * _bg_s;
                    yj = _bg_yo + j * _bg_s;

                    if (fabs(xi - xo) < 1E-10 && fabs(yj - yo) < 1E-10)
                    {
                        // lower left corner touches edge start
                        invalidate_cell(i - 1, j - 1, N, W, S, E, ied);
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                        i--; // move left
                    }
                    else if (fabs(xi - xo + _bg_s) < 1E-10 && fabs(yj - yo) < 1E-10)
                    {
                        // lower right corner touches edge start
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                        invalidate_cell(i + 1, j - 1, N, W, S, E, ied);
                        invalidate_cell(i + 1, j, N, W, S, E, ied);
                        if (fabs(fabs((yj - yo + _bg_s) / (xi - xo)) - fabs(slope)) < 1E-10)
                        {
                            // upper left corner touches edge
                            invalidate_cell(i - 1, j, N, W, S, E, ied);
                        }
                        j++; // move up
                    }
                    else if (fabs(xi - xo + _bg_s) < 1E-10 && fabs(yj - yo + _bg_s) < 1E-10)
                    {
                        // upper right corner touches edge start
                        invalidate_cell(i + 1, j, N, W, S, E, ied);
                        invalidate_cell(i + 1, j + 1, N, W, S, E, ied);
                        j++; // move up
                    }
                    else if (fabs(xi - xo) < 1E-10 && fabs(yj - yo + _bg_s) < 1E-10)
                    {
                        // upper left corner touches edge start
                        invalidate_cell(i - 1, j, N, W, S, E, ied);
                        invalidate_cell(i, j + 1, N, W, S, E, ied);
                        i--; j++; // move left and up
                    }
                    else if (fabs(yj - yo) < 1E-10)
                    {
                        // start of that sdge lies on the left border of that cell
                        invalidate_cell(i, j - 1, N, W, S, E, ied);
                    }
                    #pragma endregion

                    while (true)
                    {
                        // move from io, jo to in, jn to the right of the line
                        invalidate_cell(i, j, N, W, S, E, ied);
                        if (i <= in && j >= jn) break;

                        // upper left corner of cell
                        xi = _bg_xo + i * _bg_s;
                        yj = _bg_yo + (j + 1) * _bg_s;

                        if (fabs((yj - yo) / (xi - xo)) < fabs(slope) + 1E-10)
                        {
                            if (fabs(fabs((yj - yo) / (xi - xo)) - fabs(slope)) < 1E-10)
                            {
                                // edge passes through the upper left corner
                                invalidate_cell(i - 1, j, N, W, S, E, ied);
                                invalidate_cell(i, j + 1, N, W, S, E, ied);
                                i--; // move left
                            }
                            j++; // move up
                        }
                        else i--; // move left

                    }
                    #pragma endregion
                }
                #pragma endregion
            }
        }
    }
    return 0;
    #pragma endregion
}

int MeshingGenie_2d::sprinkle_points_along_boundaries(double d)
{
  #pragma region Sprinkle Points Along Boundaries:
  std::list<Point*> boundary_points;
  std::list<size_t> boundary_edges;

  double max_dist_sq = 4 * d * d;

  size_t num_directed_edges(_num_boundary_edges + _num_holes_edges);

  size_t num_edges(_edges_x1.size());
    for (size_t ied = 0; ied < num_edges; ied++)
    {
        double xo(_edges_x1[ied]), yo(_edges_y1[ied]);
        double xn(_edges_x2[ied]), yn(_edges_y2[ied]);

    // construct an open poly representing that edge:
    PolyPoint* pst = new PolyPoint();
    PolyPoint* pend = new PolyPoint();
    pst->x = xo; pst->y = yo; pst->next = pend; pst->prev = 0;
    pend->x = xn; pend->y = yn; pend->prev = pst; pend->next = 0;

    PolyPoint* pi = pst;
    PolyPoint* pj;
    while (true)
    {
      #pragma region Split edge randomly till a maximal distribution is obtained:

      double Lsq = distance_squared(pi->x - pi->next->x, pi->y - pi->next->y);

      if (Lsq > max_dist_sq)
      {
        double L = sqrt(Lsq);
        double ds = L - 2 * d;
        double u = generate_a_random_number();
        double tmpL = d + u * ds;
        double Linv = 1.0 / L;
        double xx = pi->x + tmpL * Linv * ( pi->next->x - pi->x);
        double yy = pi->y + tmpL * Linv * ( pi->next->y - pi->y);

        pj = new PolyPoint();
        pj->x = xx; pj->y = yy;
        pj->next = pi->next; pi->next->prev = pj;
        pj->prev = pi; pi->next = pj;
      }
      else if (pi->next == pend) break;
      else pi = pi->next;

      #pragma endregion
    }

    pi = pst; Point* newpoint; Point* oldpoint;
    double dst_sq; size_t icell;
    while (true)
    {
      newpoint = closest_point(pi->x, pi->y, 1E-10, dst_sq, icell);

      if (newpoint == 0)
      {
        newpoint = new Point();
        newpoint->x = pi->x; newpoint->y = pi->y;
        newpoint->next = 0; newpoint->edges = 0;
        newpoint->done = false; newpoint->num_edges = 0;
        newpoint->on_boundary = true;

        if (_cell_points[icell] == 0) _cell_points[icell] = newpoint;
        else
        {
          Point* q = _cell_points[icell];
          while (q->next != 0) q = q->next;
          q->next = newpoint;
        }
        _num_sprinkled++;
        _bad_cells[icell] = true; _num_bad_cells++;
      }

      if (pi == pst) { pi = pi->next; oldpoint = newpoint; continue;}
      if (_save_point_cloud)
      {
        boundary_points.push_back(oldpoint);
        boundary_points.push_back(newpoint);
        boundary_edges.push_back(ied);
        boundary_edges.push_back(ied);
      }
      add_constraint_edge(oldpoint, newpoint, ied);
      oldpoint = newpoint;

      if (pi == pend) break;

      pi = pi->next;
    }
    delete_poly(pst);
  }
  if (_save_point_cloud) save_boundary_points("boundary_points.dat", boundary_points, boundary_edges);
  return 0;
  #pragma endregion
}

int MeshingGenie_2d::apply_maximal_poisson_disk_sampler()
{
    #pragma region Maximal Poisson Disk Sampler:

    clock_t start_time, end_time; double cpu_time;

  print_message("Mode I: ");
  start_time = clock();

    size_t ngcells(0), ngbadcells(0), num_bad_cells(0);

    size_t icell, icell_max(_num_cells - 1);
    size_t num_misses(0);

  size_t i, j, num_invalidated(0);
    double xx, yy;
  size_t ipnt(0);

    //while (_num_bad_cells != _num_cells)

    size_t nit_max = 5 * (_num_cells - _num_bad_cells);

  //size_t mit(0);std::vector<double> percentage;

  for (size_t nit = 0; nit < nit_max; nit++)
  {
    if (num_misses > 300) break;

        double r = generate_a_random_number();
        icell = static_cast<size_t>(floor(r * icell_max + 0.5));

    if (_bad_cells[icell]) continue;

        // A cell with some available space - pick a random point
        double u = generate_a_random_number();
        double v = generate_a_random_number();

        // check if we can add that point to that cell
        get_cell_indices(icell, i, j);
    double xo = _bg_xo + i * _bg_s;
    double yo = _bg_yo + j * _bg_s;
        xx = xo + u * _bg_s;
        yy = yo + v * _bg_s;

        // check if adding this point is OK
    bool ok(false);
    if (_blue_cells[icell])
    {
      // Check if point is in the domain:
      size_t ied;
      double ed_x1, ed_y1, ed_x2, ed_y2, ed_x3, ed_y3, ed_x4, ed_y4;
      _cell_edge_iter = _cell_second_edge_map.find(icell);
      if (_cell_edge_iter != _cell_second_edge_map.end())
      {
        // cell is intersected using two edges
        ied = _cell_edge_iter->second;
        ed_x1 = _edges_x1[ied]; ed_y1 = _edges_y1[ied];
        ed_x2 = _edges_x2[ied]; ed_y2 = _edges_y2[ied];

        _cell_edge_iter = _cell_first_edge_map.find(icell);
        ied = _cell_edge_iter->second;
        ed_x3 = _edges_x1[ied]; ed_y3 = _edges_y1[ied];
        ed_x4 = _edges_x2[ied]; ed_y4 = _edges_y2[ied];

        bool non_manifold(false);
        if (distance_squared(ed_x3 - ed_x2, ed_y3 - ed_y2) < 1E-10)
        {
          ed_x3 = _edges_x2[ied]; ed_y3 = _edges_y2[ied];
        }
        else if (distance_squared(ed_x4 - ed_x1, ed_y4 - ed_y1) < 1E-10)
        {
          double tmp = ed_x1; ed_x1 = ed_x3; ed_x3 = ed_x2; ed_x2 = tmp;
          tmp = ed_y1; ed_y1 = ed_y3; ed_y3 = ed_y2; ed_y2 = tmp;
        }
        else
        {
          // Cutting edges do not intersect
          non_manifold = true;
        }

        if (non_manifold)
        {
          #pragma region Non_manifold Case:
          double area = area_triangle(ed_x1, ed_y1, ed_x2, ed_y2, ed_x3, ed_y3);
          bool valid_poly(false);
          if (area > 0) valid_poly = true;
          bool in_poly = point_in_triangle(xx, yy, ed_x1, ed_y1, ed_x2, ed_y2, ed_x3, ed_y3);
          if (!in_poly) in_poly = point_in_triangle(xx, yy, ed_x1, ed_y1, ed_x3, ed_y3, ed_x4, ed_y4);
          if (in_poly == valid_poly) ok = true;
          #pragma endregion
        }
        else
        {
          #pragma region Manifold case:
          double area_1 = area_triangle(xx, yy, ed_x1, ed_y1, ed_x2, ed_y2);
          double area_2 = area_triangle(xx, yy, ed_x2, ed_y2, ed_x3, ed_y3);

          if (area_1 > 0.0 && area_2 > 0.0) ok = true;
          else if (area_1 < 0.0 && area_2 < 0.0) ok = false;
          else
          {
            bool flat_angle(false);
            if (crossing_segments(xo, yo, xo + _bg_s, yo, ed_x1, ed_y1, ed_x3, ed_y3)) flat_angle = true;
            else if (crossing_segments(xo + _bg_s, yo, xo + _bg_s, yo + _bg_s, ed_x1, ed_y1, ed_x3, ed_y3)) flat_angle = true;
            else if (crossing_segments(xo + _bg_s, yo + _bg_s, xo, yo + _bg_s, ed_x1, ed_y1, ed_x3, ed_y3)) flat_angle = true;
            else if (crossing_segments(xo, yo + _bg_s, xo, yo, ed_x1, ed_y1, ed_x3, ed_y3)) flat_angle = true;

            if (flat_angle)
            {
              if (point_in_triangle(xx, yy, ed_x1, ed_y1, ed_x2, ed_y2, ed_x3, ed_y3))
              {
                double area = area_triangle(ed_x1, ed_y1, ed_x2, ed_y2, ed_x3, ed_y3);
                if (area > 0) ok = true;
              }
              else
              {
                double area = area_triangle(xx, yy, ed_x1, ed_y1, ed_x3, ed_y3);
                if (area > 0) ok = true;
              }
            }
            else
            {
              double area = area_triangle(ed_x1, ed_y1, ed_x2, ed_y2, ed_x3, ed_y3);
              bool valid_tri(false);
              if (area > 0) valid_tri = true;
              bool in_tri = point_in_triangle(xx, yy, ed_x1, ed_y1, ed_x2, ed_y2, ed_x3, ed_y3);
              if (in_tri == valid_tri) ok = true;
            }
          }
          #pragma endregion
        }
      }
      else
      {
        // cell is intersected using one edge
        _cell_edge_iter = _cell_first_edge_map.find(icell);
        ied = _cell_edge_iter->second;
        ed_x1 = _edges_x1[ied]; ed_y1 = _edges_y1[ied];
        ed_x2 = _edges_x2[ied]; ed_y2 = _edges_y2[ied];

        double area = area_triangle(xx, yy, ed_x1, ed_y1, ed_x2, ed_y2);
        if (area > 0.0) ok = true;
      }
    }
    else ok = true;

    // check if this point is not inside a circle of a neighbor
    if (ok) ok = is_valid_internal_point(icell, xx, yy);

        if (!ok)
        {
      num_misses++;

      //mit++;
      //percentage.push_back(_num_sprinkled * 1.0 / mit);
      continue;
        }

    num_misses = 0;
        Point* p = new Point();
        p->x = xx; p->y = yy;
        p->next = 0; p->edges = 0;
    p->done = false; p->on_boundary = false; p->num_edges = 0;
    Point* q = _cell_points[icell];
    if (q == 0) _cell_points[icell] = p;
    else
    {
      while (q->next != 0) q = q->next;
      q->next = p;
    }
    _num_sprinkled++;
        _bad_cells[icell] = true; _num_bad_cells++; num_invalidated++;

    //mit++;
    //percentage.push_back(_num_sprinkled * 1.0 / mit);
    }

  end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC; start_time = clock();
    print_message_end("^ Execution Time = " + double_to_string(cpu_time) + " seconds.");

    print_message_end("^ Number of sprinkled points  = "
                    + double_to_string(1.0 * _num_sprinkled));
    print_message_end("^ Number of remaining cells  = "
                    + double_to_string(1.0 * (_num_cells - _num_bad_cells)));


  //plot_cells("valid_cells_ModeI.ps", false);

  print_message("Mode II: ");

  _circles = new Point*[120];
  for (int ii = 0; ii < 120; ii++) _circles[ii] = 0;

  start_time = clock();

  eliminate_voids_via_polys();

  end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC; start_time = clock();

  print_message_end("^ Execution Time = " + double_to_string(cpu_time) + " seconds.");

    print_message_end("^ Number of sprinkled points  = "
                    + double_to_string(1.0 * _num_sprinkled));

    print_message_end("^ Number of remaining cells  = "
                    + double_to_string(1.0 * (_num_cells - _num_bad_cells)));


  delete [] _circles;
    return 0;
    #pragma endregion
}

int MeshingGenie_2d::fill_neighbors_data()
{
  #pragma region Filling neighbors data:

  _num_neighbor_points = 100;
  _neighbor_points = new Point*[100];
  _neighbor_points_cells = new size_t[_num_neighbor_points];
  for (size_t i = 0; i < _num_neighbor_points; i++) _neighbor_points[i] = 0;

  _num_neighbor_cells = 21;
  _neighbors = new int[_num_neighbor_cells];

  _neighbors[0] = 0;
  _neighbors[1] = int(_nr) - 1;          // E_1
  _neighbors[2] = 1;                // N_1
  _neighbors[3] = -(int(_nr) - 1);       // W_1
  _neighbors[4] = -1;               // S_1

  _neighbors[5] =  (int(_nr) - 1) + 1;       // NE_1
  _neighbors[6] = -(int(_nr) - 1) + 1;       // NW_1
  _neighbors[7] = -(int(_nr) - 1) - 1;       // SW_1
  _neighbors[8] =  (int(_nr) - 1) - 1;       // SE_1

  _neighbors[9]  = 2 * (int(_nr) - 1) - 1;   // ES_2
  _neighbors[10]  = 2 * (int(_nr) - 1);       // EC_2
  _neighbors[11] = 2 * (int(_nr) - 1) + 1;   // EN_2

  _neighbors[12] = 2 + (int(_nr) - 1);       // NE_2
  _neighbors[13] = 2;                   // NC_2
  _neighbors[14] = 2 - (int(_nr) - 1);       // NW_2

  _neighbors[15] = -2 * (int(_nr) - 1) + 1;   // WN_2
  _neighbors[16] = -2 * (int(_nr) - 1);       // WC_2
  _neighbors[17] = -2 * (int(_nr) - 1) - 1;   // WS_2

  _neighbors[18] = -2 - (int(_nr) - 1);       // SW_2
  _neighbors[19] = -2;                   // SC_2
  _neighbors[20] = -2 + (int(_nr) - 1);       // SE_2

  return 0;
    #pragma endregion
}

int MeshingGenie_2d::refill_neighbors_data()
{
  #pragma region Filling neighbors data:
  delete [] _neighbors;

  _num_neighbor_cells = 45;
  _neighbors = new int[45];

  int di = int(_nr) - 1; int dj = 1;

  _neighbors[0] = 0;
  _neighbors[1] = di;
  _neighbors[2] = di + dj;
  _neighbors[3] = dj;
  _neighbors[4] = - di + dj;
  _neighbors[5] = - di;
  _neighbors[6] = - di - dj;
  _neighbors[7] = - dj;
  _neighbors[8] = di - dj;

  _neighbors[9]  = 2 * di - dj;
  _neighbors[10]  = 2 * di;
  _neighbors[11] = 2 * di + dj;
  _neighbors[12] = 2 * di + 2 * dj;
  _neighbors[13] = di + 2 * dj;
  _neighbors[14] = 2 * dj;
  _neighbors[15] = - di + 2 * dj;
  _neighbors[16] = -2 * di + 2 * dj;
  _neighbors[17] = -2 * di + dj;
  _neighbors[18] = -2 * di;
  _neighbors[19] = -2 * di - dj;
  _neighbors[20] = -2 * di - 2 * dj;
  _neighbors[21] = - di - 2 * dj;
  _neighbors[22] = - 2 * dj;
  _neighbors[23] = di - 2 * dj;
  _neighbors[24] = 2 * di - 2 * dj;

  _neighbors[25] = 3 * di - 2 * dj;
  _neighbors[26] = 3 * di - dj;
  _neighbors[27] = 3 * di;
  _neighbors[28] = 3 * di + dj;
  _neighbors[29] = 3 * di + 2 * dj;

  _neighbors[30] = 2 * di + 3 * dj;
  _neighbors[31] = di + 3 * dj;
  _neighbors[32] = 3 * dj;
  _neighbors[33] = - di + 3 * dj;
  _neighbors[34] = -2 * di + 3 * dj;

  _neighbors[35] = - 3 * di + 2 * dj;
  _neighbors[36] = - 3 * di + dj;
  _neighbors[37] = - 3 * di;
  _neighbors[38] = - 3 * di - dj;
  _neighbors[39] = - 3 * di - 2 * dj;

  _neighbors[40] = - 2 * di - 3 * dj;
  _neighbors[41] = - di - 3 * dj;
  _neighbors[42] = - 3 * dj;
  _neighbors[43] = di - 3 * dj;
  _neighbors[44] = 2 * di - 3 * dj;

  return 0;
    #pragma endregion
}

// NEED TO MAKE SURE THAT PI-PJ DOES NOT INTERSECT WITH ANY BOUNDARY EDGES
bool MeshingGenie_2d::is_valid_internal_point(size_t icell, double xx, double yy)
{
  #pragma region check if this point is not inside a circle of a neighbor
  bool first(true); Point* p;
  for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
  {
    size_t jcell = icell + _neighbors[ii];

    if (_cell_points[jcell] == 0) continue; // _bad cells due to boundaries

    if (first) p = _cell_points[jcell];

    // retrieve the point associated with that cell
    double xj = p->x;
    double yj = p->y;
    double dx = xj - xx;
    double dy = yj - yy;

    if (dx < _dm && dy < _dm && dx * dx + dy * dy < _dm_squared)
    {
      // NEED TO MAKE SURE THAT PI-PJ DOES NOT INTERSECT WITH ANY BOUNDARY EDGES
      return false;
    }

    if (p->next != 0)
    {
      p = p->next;
      first = false;
      ii--;
    }
    else first = true;
  }
  return true;
  #pragma endregion
}

int MeshingGenie_2d::eliminate_voids_via_polys()
{
  #pragma region Utilize Polys to capture the remaining voids:
  _failed = false;

  // construct Polys
  PolyPoint* pp; size_t num_poly_pnts(0);
  std::list<Poly*> polys_list;
  std::list<Poly*>::iterator iter;
  double total_area(0.0); size_t num_bad_cells(0), num_blue_cells(0);
  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    #pragma region Create Polys:
    if (_bad_cells[icell]) {num_bad_cells++; continue;}
    //if (_blue_cells[icell]) {num_blue_cells++; continue;}

    generate_approximate_poly(icell, pp, num_poly_pnts);

    if (num_poly_pnts == 0)
    {
      _bad_cells[icell] = true; _num_bad_cells++; num_bad_cells++;
      continue;
    }

    bool ok(false);

    _icell = icell;
    //plot_cell("void.ps", icell, pp, 0);
    if (split_polygon(pp))
    {
      for (size_t ii = 0; ii < _num_polys; ii++)
      {
        if (_polys[ii] == 0) continue;

        pp = _polys[ii];
        Poly* ply = new Poly();
        ply->icell = icell; ply->area = 0.0;
        PolyPoint* po = pp;
        while (true)
        {
          ply->area += pp->x * (pp->next->y - pp->prev->y);
          pp = pp->next;
          if (pp == po) break;
        }
        if (ply->area < 0.0)
        {
          plot_cell("invalid_split_void.ps", icell, pp, 0);
        }
        ply->pp = pp;
        total_area += ply->area;
        polys_list.push_back(ply);
        ok = true;
      }
    }
    else
    {
      // plot_cell("void.ps", icell, pp, 0);
      Poly* ply = new Poly();
      ply->icell = icell; ply->area = 0.0;
      PolyPoint* po = pp;
      while (true)
      {
        ply->area += pp->x * (pp->next->y - pp->prev->y);
        pp = pp->next;
        if (pp == po) break;
      }
      ply->pp = pp;
      if (ply->area < 0.0)
      {
        plot_cell("invalid_one_void.ps", icell, pp, 0);
      }
      total_area += ply->area;
      polys_list.push_back(ply);
      ok = true;
    }
    if (!ok)
    {
      std::cout << "Error! Could not generate an approximate poly!" << std::endl;
      plot_cell("not_ok_void.ps", icell, pp, 0);
    }
    #pragma endregion
  }

  for (size_t ii = 0; ii < _num_polys; ii++) _polys[ii] = 0;

  if (polys_list.size() == 0) return 0;

  double area_inv(1.0 / total_area);
  double area_prev(0.0);

  size_t ipoly(0);
  Poly** polys = new Poly*[polys_list.size()];
  bool* invalid = new bool[polys_list.size()];
  for (iter = polys_list.begin(); iter!= polys_list.end(); iter++)
  {
    (*iter)->area *= area_inv;
    (*iter)->area += area_prev;
    area_prev = (*iter)->area;
    polys[ipoly] = *iter;
    invalid[ipoly] = false;
    ipoly++;
  }

  //std::cout<< " ... average area = " << total_area / ipoly / _bg_s / _bg_s << std::endl;

  //plot_cells("valid_cells_ModeII.ps", false);

  double u;
  size_t num_crowded(0);
  size_t num_polys(ipoly), num_invalidated(0);
  size_t lower, upper, ipoly_max(num_polys - 1);
  double invalidated_area(0.0);

  size_t nit_max = 2 * num_polys;
  if (nit_max < 100 ) nit_max = 100;

  size_t num_misses(0);
  for (size_t nit = 0; nit < nit_max; nit++)
    {
    if (num_invalidated == num_polys) break;

    if (invalidated_area > 0.7) break;

    if (num_misses > 25) break;

    u = generate_a_random_number();

    // binary search for the selected poly
    lower = 0; upper = ipoly_max;
    while (true)
    {
      #pragma region Binary Search:
      if (u < polys[0]->area)
      {
        ipoly = 0;
        break;
      }

      ipoly = (lower + upper) / 2;
      if (polys[ipoly]->area < u)
      {
        lower = ipoly;
      }
      else
      {
        upper = ipoly;
      }

      if (upper == lower + 1)
      {
        ipoly = upper; break;
      }
      #pragma endregion
    }

    if (invalid[ipoly]) continue;

    if (_bad_cells[polys[ipoly]->icell])
    {
      // delete that poly
      invalid[ipoly] = true; num_invalidated++;
      if (ipoly == 0) invalidated_area += polys[ipoly]->area;
      else invalidated_area += polys[ipoly]->area - polys[ipoly - 1]->area;
      continue;
    }

    //plot_cell("void.ps", polys[ipoly]->icell, polys[ipoly]->pp, 0);

    bool point_inserted = insert_point_in_poly(polys[ipoly]->pp, polys[ipoly]->icell);
    if (point_inserted)
    {
      // Invalidated poly and cell:
      invalid[ipoly] = true; num_invalidated++;
      if (ipoly == 0) invalidated_area += polys[ipoly]->area;
      else invalidated_area += polys[ipoly]->area - polys[ipoly - 1]->area;
      num_misses = 0;
    }
    else num_misses++;
  }

  // clear memory
  for (ipoly = 0; ipoly < num_polys; ipoly++)
  {
    delete_poly(polys[ipoly]->pp);
    delete polys[ipoly];
  }

  delete[] invalid;
  delete[] polys;
  polys_list.clear();

  print_message_end("^ Number of generated polys  = "
                  + double_to_string(1.0 * num_polys));

  print_message_end("^ Number of invalidated polys  = "
                  + double_to_string(1.0 * num_invalidated));

  if (num_invalidated != num_polys) eliminate_voids_via_polys();

  return 0;
  #pragma endregion
}

void MeshingGenie_2d::pick_a_random_point_from_a_convex_hull(PolyPoint* pp, double &xx, double &yy)
{
  #pragma region Pick a random point From a convex hull:
  double u_sum(0.0); PolyPoint* pst(pp);
  while (true)
  {
    pp->u = generate_a_random_number();
    u_sum += pp->u;
    pp = pp->next;
    if (pp == pst) break;
  }
  u_sum = 1.0 / u_sum;
  while (true)
  {
    pp->u *= u_sum;
    pp = pp->next;
    if (pp == pst) break;
  }
  xx = 0.0; yy = 0.0;
  while (true)
  {
    xx += pp->u * pp->x;
    yy += pp->u * pp->y;
    pp = pp->next;
    if (pp == pst) break;
  }
  #pragma endregion
}

bool MeshingGenie_2d::cut_polygon(PolyPoint* &pp, size_t &num_poly_pnts, Point** circles, size_t num_circles)
{
  #pragma region Cut a polygon using a number of circles:

  // The polygon is assumed to have vertical and horizontal edges only
  bool cut(false);

  double xc, yc, dx, dy, dr;
  PolyPoint* po;
  PolyPoint* pL;
  PolyPoint* pR;
  PolyPoint* e_st; PolyPoint* e_end; PolyPoint* newpoint;

  int nit(0);
  while (true)
  {
    nit++;

    if (nit > 100)
    {
      std::cout << "*** Too many iteration is cut polygon!!! ****" << std::endl;
      break;
    }

    bool nothing_new = true;

    for (size_t ic = 0; ic < num_circles; ic++)
    {
      if (circles[ic] == 0) continue; // invaid circle pointer

      xc = circles[ic]->x; yc = circles[ic]->y;

      // check if all polypoints lies inside the circle
      bool invalid_poly(true);
      for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
      {
        dx = pp->x - xc; dy = pp->y - yc;

        if (dx * dx + dy * dy > _dm_squared + 1E-10) {invalid_poly = false; break;}

        pp = pp->next;
      }

      if (invalid_poly)
      {
        // All the points of the polyline lies in this circle
        for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
        {
          po = pp;
          pp = pp->next;
          delete po;
        }
        num_poly_pnts = 0; pp = 0;
        return true;
      }

      // loop over the poly points
      for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
      {
        if (pp->center_L == circles[ic] || pp->center_R == circles[ic])
        {
          pp = pp->next; continue;
        }

        dx = pp->x - xc; dy = pp->y - yc;

        if (dx * dx + dy * dy > _dm_squared - 1E-10)
        {
          pp = pp->next; continue;
        }

        // pp point is inside circle:

        bool poly_changed(false);
        pL = 0; pR = 0;

        #pragma region moving backward:
        e_st = pp; e_end = pp->prev;
        while (true) // exit if an intersection is found or if all points lies in the circle
        {
          // check if the next point lies in the same circle
          dx = xc - e_end->x; dy = yc - e_end->y;
          dr = dx * dx + dy * dy;
          if (e_end->center_R != circles[ic] && dr < _dm_squared - 1E-10)
          {
            // proceed to the prev segment
            e_st = e_st->prev;
            e_end = e_end->prev;
          }
          else if (e_end->center_R == circles[ic] || dr < _dm_squared + 1E-10)
          {
            // The end point is the point of intersection
            e_end->next = pp; pp->prev = e_end;
            e_end->center_R = circles[ic];
            // delete all the points between new_point and pp
            while (e_st != pp)
            {
              po = e_st;
              e_st = e_st->next;
              delete po; num_poly_pnts--;
            }
            break;
          }
          else
          {
            nothing_new  = false; poly_changed = true;
            // calculate intersection point
            double tx(e_end->x - e_st->x);
            double ty(e_end->y - e_st->y);
            double t = sqrt(tx * tx + ty * ty);
            tx /= t; ty /= t;
            double ax(e_end->x - xc);
            double ay(e_end->y - yc);
            double aproj(ax * tx + ay * ty);
            double xst(e_end->x - aproj * tx);
            double yst(e_end->y - aproj * ty);
            double dxx(xc - xst), dyy(yc - yst);
            double dt = sqrt(_dm_squared - dxx * dxx - dyy * dyy);

            xst += dt * tx; yst += dt *ty; // Intersection point

            newpoint = new PolyPoint(); num_poly_pnts++;
            newpoint->x = xst; newpoint->y = yst;
            newpoint->next = pp; pp->prev = newpoint;
            newpoint->prev = e_end; e_end->next = newpoint;
            newpoint->center_R = circles[ic];
            newpoint->center_L = 0;

            // mark the new point as pL if e_st lies on a left circle
            if (e_st->center_L != 0)
            {
              pL = newpoint;
              pL->center_L = pp->center_L;
            }

            // delete all the points between new_point and pp
            while (e_st != pp)
            {
              po = e_st;
              e_st = e_st->next;
              delete po; num_poly_pnts--;
            }
            break;
          }
        }
        #pragma endregion

        #pragma region moving forward:
        e_st = pp; e_end = pp->next;
        while (true) // exist if an intersection is found or if all points lies in the circle
        {
          // check if the next point lies in the same circle
          dx = xc - e_end->x; dy = yc - e_end->y;
          dr = dx * dx + dy * dy;
          if (e_end->center_L != circles[ic] && dr < _dm_squared - 1E-10)
          {
            // proceed to the next segment
            e_st = e_st->next;
            e_end = e_end->next;
            if (e_end == pp)
            {
              // All the points of the polyline lies in this circle
              for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
              {
                po = e_end;
                e_end = e_end->next;
                delete po;
              }
              num_poly_pnts = 0;
              break;
            }
          }
          else if (e_end->center_L == circles[ic] || dr < _dm_squared + 1E-10)
          {
            // end point is the intersection point
            e_end->prev = pp; pp->next = e_end;
            e_end->center_L = circles[ic];

            // delete all the points between new_point and pp
            while (e_st != pp)
            {
              po = e_st;
              e_st = e_st->prev;
              delete po; num_poly_pnts--;
            }
            pp->prev->next = pp->next;
            pp->next->prev = pp->prev;
            delete pp; num_poly_pnts--;
            pp = e_end;
            break;
          }
          else
          {
            nothing_new = false; poly_changed = true;

            // get intersection point
            double tx(e_end->x - e_st->x);
            double ty(e_end->y - e_st->y);
            double t = sqrt(tx * tx + ty * ty);
            tx /= t; ty /= t;
            double ax(e_end->x - xc);
            double ay(e_end->y - yc);
            double aproj(ax * tx + ay * ty);
            double xst(e_end->x - aproj * tx);
            double yst(e_end->y - aproj * ty);
            double dxx(xc - xst), dyy(yc - yst);
            double dt = sqrt(_dm_squared - dxx * dxx - dyy * dyy);

            xst += dt * tx; yst += dt * ty; // Intersection point

            newpoint = new PolyPoint(); num_poly_pnts++;
            newpoint->x = xst; newpoint->y = yst;
            newpoint->next = e_end; e_end->prev = newpoint;
            newpoint->prev = pp; pp->next = newpoint;
            newpoint->center_L = circles[ic];
            newpoint->center_R = 0;

            // activate P_R if e_st lies on a right circle
            if (e_st->center_R != 0)
            {
              pR = newpoint;
              pR->center_R = e_st->center_R;
            }

            // delete all the points between new_point and pp
            while (e_st != pp)
            {
              po = e_st;
              e_st = e_st->prev;
              delete po; num_poly_pnts--;
            }
            pp->prev->next = pp->next;
            pp->next->prev = pp->prev;
            delete pp; num_poly_pnts--;
            pp = e_end;
            break;
          }
        }
        #pragma endregion

        if (poly_changed)
        {
          cut = true;

          bool jammed(false);
          if (pL != 0 && pR != 0)
          {
            double xi(pL->center_L->x), yi(pL->center_L->y);
            double xj(pR->center_R->x), yj(pR->center_R->y);

            // check if the remaining poly lies in one of these two circles
            for (int jcirc = 0; jcirc < 2; jcirc++)
            {
              invalid_poly = true;
              for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
              {
                if (jcirc == 0) {dx = pp->x - xi; dy = pp->y - yi;}
                else {dx = pp->x - xj; dy = pp->y - yj;}

                if (dx * dx + dy * dy > _dm_squared + 1E-10) {invalid_poly = false; break;}

                pp = pp->next;
              }
              if (invalid_poly) break;
            }
            if (invalid_poly)
            {
              // All the points of the polyline lies in this circle
              for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
              {
                po = pp;
                pp = pp->next;
                delete po;
              }
              num_poly_pnts = 0; pp = 0;
              return true;
            }


            for (size_t ipnt = 0; ipnt < num_poly_pnts; ipnt++)
            {
              dx = pp->x - xi; dy = pp->y - yi;

              if (dx * dx + dy * dy > _dm_squared + 1E-10) {invalid_poly = false; break;}

              pp = pp->next;
            }


            double dx = xj - xi; double dy = yj - yi;
            double d_squared = dx * dx + dy * dy;

            if ( d_squared < 4 * _dm_squared - 1E-10)
            {
              // two intersecting circles
              double h_squared = _dm_squared - 0.25 * d_squared;
              double h_d = sqrt( h_squared / d_squared);

              double x_R = 0.5 * (xi + xj) + h_d * dy;
              double y_R = 0.5 * (yi + yj) - h_d * dx;

              dx = x_R - xc; dy = y_R - yc;
              if (dx * dx + dy * dy < _dm_squared)
              {
                // Right intersection is inside circle
                dx = xj - xi; dy = yj - yi;
                double x_L = 0.5 * (xi + xj) - h_d * dy;
                double y_L = 0.5 * (yi + yj) + h_d * dx;

                dx = x_L - xc; dy = y_L - yc;
                if (dx * dx + dy * dy > _dm_squared && point_in_convex_hull(x_L, y_L, pp))
                {
                  // Left intersection is out of circle
                  jammed = true;
                  pR->x = x_L; pR->y = y_L;
                  // delete pL
                  pR->prev = pL->prev; pL->prev->next = pR;
                  delete pL; num_poly_pnts--;
                }
              }
            }
          }
          if (!jammed)
          {
            for (int ii = 0; ii < 2; ii++)
            {
              if (ii == 0 && pL == 0) continue;
              if (ii == 1 && pR == 0) continue;
              if (ii == 1) pL = pR;

              // Move pL to the Left intersection of the two circles
              double xi(pL->center_L->x), yi(pL->center_L->y);
              double xj(pL->center_R->x), yj(pL->center_R->y);

              double dx = xj - xi; double dy = yj - yi;

              double d_squared = dx * dx + dy * dy;
              if ( d_squared > 4 * _dm_squared - 1E-10)
              {
                plot_cell("adjusted_poly.ps", _icell, pp, 0);
                std::cout << "*** Circles do not intersect!!! ****" << std::endl;
                continue;
              }
              // two intersecting circles
              double h_squared = _dm_squared - 0.25 * d_squared;
              double h_d = sqrt( h_squared / d_squared);

              double x_L = 0.5 * (xi + xj) - h_d * dy;
              double y_L = 0.5 * (yi + yj) + h_d * dx;

              if (point_in_convex_hull(x_L, y_L, pp))
              {
                pL->x = x_L; pL->y = y_L;
              }
            }
          }
          break;
        }
      }
    }
    if (nothing_new) break;
  }
  if (_max_num_poly_pnts == 0 || num_poly_pnts > _max_num_poly_pnts) _max_num_poly_pnts = num_poly_pnts;
  return cut;
  #pragma endregion
}

bool MeshingGenie_2d::split_polygon(PolyPoint* pp)
{
  #pragma region split the polygon to capture all voids better:

  for (size_t i = 0; i < _num_polys; i++) _polys[i] = 0;
  bool changed(false);
  size_t active_poly(0), ipoly(1);
  while (true)
  {
    // split polyline pp
    bool split(false);

    PolyPoint* pi = pp->next; size_t num_poly_pnts(1);
    while (pi != pp)
    {
      pi = pi->next;  num_poly_pnts++;
    }
    PolyPoint* pj;
    for (size_t iedge = 0; iedge < num_poly_pnts; iedge++)
    {
      pi = pi->next;
      if (pi->center_R == 0) continue;

      double xi = pi->center_R->x;
      double yi = pi->center_R->y;
      pj = pi->next;
      while (pj->next->next != pi)
      {
        pj = pj->next;
        if (pj->center_R == 0) continue;

        if (pi->center_R == pj->center_R)
        {
          // a non manifold case
          bool valid(false);
          PolyPoint* tmp = pi->next;
          pi->next = pj->next; pj->next->prev = pi;
          if (is_valid_poly(pi))
          {
            _polys[active_poly] = pi;
            valid = true;
          }

          pj->next = tmp; tmp->prev = pj;
          if (is_valid_poly(pj))
          {
            if (valid)
            {
              _polys[ipoly] = pj; ipoly++;
            }
            else
            {
              _polys[active_poly] = pj;
            }
          }
          else
          {
            delete_poly(pj);
          }

          if (!valid)
          {
            delete_poly(pi);
          }

          if (ipoly > 4 || _failed)
          {
            plot_cell("silly_cell.ps", _icell, 0, 0);
            _failed = true;
            return true;
          }

          changed = true;
          split = true; break;
        }

        // check the intersection of circles i and j
        double xj = pj->center_R->x;
        double yj = pj->center_R->y;

        double dx = xj - xi; double dy = yj - yi;

        if (dx > 2 *_dm - 1E-10) continue;
        if (dy > 2 *_dm - 1E-10) continue;

        double d_squared = dx * dx + dy * dy;
        if (d_squared > 4 * _dm_squared - 1E-10) continue;

        // two intersecting circles
        double h_squared = _dm_squared - 0.25 * d_squared;
        double h_d = sqrt( h_squared / d_squared);

        double x_L = 0.5 * (xi + xj) - h_d * dy;
        double y_L = 0.5 * (yi + yj) + h_d * dx;

        double x_R = 0.5 * (xi + xj) + h_d * dy;
        double y_R = 0.5 * (yi + yj) - h_d * dx;

        if (point_in_convex_hull(x_L, y_L, pi) && point_in_convex_hull(x_R, y_R, pi))
        {
          bool valid(false); bool both_invalid(true);
          pp = pi->next;

          PolyPoint* p1 = new PolyPoint();
          p1->center_R = pj->center_R;
          p1->center_L = pi->center_R;
          p1->x = x_L; p1->y = y_L;
          p1->next = pj->next; pj->next->prev = p1;
          p1->prev = pi;       pi->next = p1;
          if (is_valid_poly(p1))
          {
            _polys[active_poly] = p1;
            valid = true; both_invalid = false;
          }

          PolyPoint* p2 = new PolyPoint();
          p2->center_R = pi->center_R;
          p2->center_L = pj->center_R;
          p2->x = x_R; p2->y = y_R;
          p2->next = pp; pp->prev = p2;
          p2->prev = pj; pj->next = p2;
          if (is_valid_poly(p2))
          {
            both_invalid = false;
            if (valid)
            {
              _polys[ipoly] = p2; ipoly++;
            }
            else
            {
              _polys[active_poly] = p2;
            }
          }
          else
          {
            //plot_cell("silly_cell.ps", _icell, p2, 0);
            delete_poly(p2);
          }

          if (!valid)
          {
            //plot_cell("silly_cell.ps", _icell, p1, 0);
            delete_poly(p1);
          }

          if (both_invalid)
          {
            _polys[active_poly] = 0;
            if (ipoly = 1) return true;
          }

          if (ipoly > 4 || _failed)
          {
            plot_cell("silly_cell.ps", _icell, 0, 0);
            _failed = true;
            return true;
          }

          changed = true;
          split = true;
          break;
        }
      }

      if (split)
      {
        pp = _polys[active_poly];
        break;
      }
    }

    if (!split)
    {
      active_poly++;
      if (active_poly == ipoly) return changed;
    }
  }
  return changed;
  #pragma endregion
}

int MeshingGenie_2d::plot_vertices(std::string file_name, bool plot_circles, bool plot_grid)
{
    #pragma region Plot Background grid:
  fstream file(file_name.c_str(), ios::out);
    file << "%!PS-Adobe-3.0" << endl;
    file << "72 72 scale     % one unit = one inch" << endl;

    #pragma region Retrieve bounding box, scale, and translate:
    double xmin(_bg_xo), ymin(_bg_yo), Lx((_nc - 1) * _bg_s), Ly((_nr - 1) * _bg_s);

    double scale_x, scale_y, scale;
    double shift_x, shift_y;

    scale_x = 6.5 / Lx;
    scale_y = 9.0 / Ly;

    if (scale_x < scale_y)
    {
        scale = scale_x;
        shift_x = 1.0 - xmin * scale;
        shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
    }
    else
    {
        scale = scale_y;
        shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
        shift_y = 1.0 - ymin * scale;
    }
    file << shift_x << " " << shift_y << " translate" << endl;
    #pragma endregion

    #pragma region Definitions of Shapes:

  file << "/circ    % stack: x y r" << endl;
    file << "{0 360 arc" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
  file << " stroke" << endl;
    file << "} def" << endl;

  file << "/fcirc    % stack: x y r" << endl;
    file << "{0 360 arc" << endl;
    file << " closepath" << endl;
    file << " gsave" << endl;
    file << " 0.6 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.01 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg2      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.008 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/quad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

  file << "/quad2      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
  file << " gsave" << endl;
    file << " 0.9 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;
    #pragma endregion

    double s(_dm * 0.05);

    std::list<Point*>::iterator iter;
  if (plot_circles)
  {
    for (size_t icell = 0; icell < _num_cells; icell++)
    {
      if (_cell_points[icell] == 0) continue;
      Point* q = _cell_points[icell];
      while (q != 0)
      {
        file << q->x * scale << "  " << q->y * scale << "  " << _dm * scale << "  ";
        file << "fcirc"     << endl;
        q = q->next;
      }
    }
  }

  if (plot_grid)
  {
    size_t i, j;
    for (size_t icell = 0; icell < _num_cells; icell++)
    {
      get_cell_indices(icell, i, j);
      double xo = _bg_xo + i * _bg_s;
      double yo = _bg_yo + j * _bg_s;
      double xn = xo + _bg_s;
      double yn = yo + _bg_s;

      // plot cell
      file << xo * scale << "  " << yo * scale << "  ";
      file << xn * scale << "  " << yo * scale << "  ";
      file << xn * scale << "  " << yn * scale << "  ";
      file << xo * scale << "  " << yn * scale << "  ";

      file << "quad"     << endl;
    }
  }

  if (plot_circles)
  {
    for (size_t icell = 0; icell < _num_cells; icell++)
    {
      if (_cell_points[icell] == 0) continue;
      Point* q = _cell_points[icell];
      while (q != 0)
      {
        file << q->x * scale << "  " << q->y * scale << "  " << _dm * scale << "  ";
        file << "circ"     << endl;
        q = q->next;
      }
    }
  }

  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_points[icell] == 0) continue;

    Point* q = _cell_points[icell];
    while (q != 0)
    {
      // plot vertex
      file << (q->x - s) * scale << "  " << (q->y - s) * scale << "  ";
      file << (q->x + s) * scale << "  " << (q->y + s) * scale << "  ";
      file << "seg2"     << endl;

      file << (q->x + s) * scale << "  " << (q->y - s) * scale << "  ";
      file << (q->x - s) * scale << "  " << (q->y + s) * scale << "  ";
      file << "seg2"     << endl;

      q = q->next;
    }
  }

    size_t num_boundary_edges(_edges_x1.size());
    for (size_t i = 0; i < num_boundary_edges; i++)
    {
        file << _edges_x1[i] * scale << "  " << _edges_y1[i] * scale << "  ";
        file << _edges_x2[i] * scale << "  " << _edges_y2[i] * scale << "  ";
        file << "seg"      << endl;
    }

    file << "showpage" << endl;

    return 0;
    #pragma endregion
}

int MeshingGenie_2d::plot_cells(std::string file_name, bool blue_cells)
{
    #pragma region Plot Background grid:
  fstream file(file_name.c_str(), ios::out);
    file << "%!PS-Adobe-3.0" << endl;
    file << "72 72 scale     % one unit = one inch" << endl;

    #pragma region Retrieve bounding box, scale, and translate:
    double xmin(_bg_xo), ymin(_bg_yo), Lx((_nc - 1) * _bg_s), Ly((_nr - 1) * _bg_s);

    double scale_x, scale_y, scale;
    double shift_x, shift_y;

    scale_x = 6.5 / Lx;
    scale_y = 9.0 / Ly;

    if (scale_x < scale_y)
    {
        scale = scale_x;
        shift_x = 1.0 - xmin * scale;
        shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
    }
    else
    {
        scale = scale_y;
        shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
        shift_y = 1.0 - ymin * scale;
    }
    file << shift_x << " " << shift_y << " translate" << endl;
    #pragma endregion

    #pragma region Definitions of Shapes:

  file << "/seg      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.01 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg2      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/quad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

  file << "/quad2      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
  file << " gsave" << endl;
    file << " 0.6 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;
    #pragma endregion

    double s(_dm * 0.05);

  if (true)
  {
    size_t i, j;
    for (size_t icell = 0; icell < _num_cells; icell++)
    {
      get_cell_indices(icell, i, j);
      double xo = _bg_xo + i * _bg_s;
      double yo = _bg_yo + j * _bg_s;
      double xn = xo + _bg_s;
      double yn = yo + _bg_s;

      // plot cell
      file << xo * scale << "  " << yo * scale << "  ";
      file << xn * scale << "  " << yo * scale << "  ";
      file << xn * scale << "  " << yn * scale << "  ";
      file << xo * scale << "  " << yn * scale << "  ";

      if (blue_cells)
      {
        if (_blue_cells[icell])
        {
          file << "quad2"     << endl;
        }
        else
        {
          file << "quad"     << endl;
        }
      }
      else
      {
        if (_bad_cells[icell])
        {
          file << "quad"     << endl;
        }
        else
        {
          file << "quad2"     << endl;
        }
      }
    }
  }

  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_points[icell] == 0) continue;
    // plot vertex
        file << (_cell_points[icell]->x - s) * scale << "  " << (_cell_points[icell]->y - s) * scale << "  ";
        file << (_cell_points[icell]->x + s) * scale << "  " << (_cell_points[icell]->y + s) * scale << "  ";
        file << "seg2"     << endl;

    file << (_cell_points[icell]->x + s) * scale << "  " << (_cell_points[icell]->y - s) * scale << "  ";
        file << (_cell_points[icell]->x - s) * scale << "  " << (_cell_points[icell]->y + s) * scale << "  ";
        file << "seg2"     << endl;
  }

    size_t num_boundary_edges(_edges_x1.size());
    for (size_t i = 0; i < num_boundary_edges; i++)
    {
        file << _edges_x1[i] * scale << "  " << _edges_y1[i] * scale << "  ";
        file << _edges_x2[i] * scale << "  " << _edges_y2[i] * scale << "  ";
        file << "seg"      << endl;
    }

    file << "showpage" << endl;

    return 0;
    #pragma endregion
}

int MeshingGenie_2d::plot_cell(std::string file_name, size_t icell, PolyPoint* pp, size_t kcell)
{
    #pragma region Plot Background grid:
  fstream file(file_name.c_str(), ios::out);
    file << "%!PS-Adobe-3.0" << endl;
    file << "72 72 scale     % one unit = one inch" << endl;

    #pragma region Retrieve bounding box, scale, and translate:
  size_t i, j;
  get_cell_indices(icell, i, j);

    double xmin(_bg_xo + i * _bg_s - 3 * _bg_s), ymin(_bg_yo + j * _bg_s - 3* _bg_s), Lx(6 * _bg_s), Ly(6 * _bg_s);

    double scale_x, scale_y, scale;
    double shift_x, shift_y;

    scale_x = 6.5 / Lx;
    scale_y = 9.0 / Ly;

    if (scale_x < scale_y)
    {
        scale = scale_x;
        shift_x = 1.0 - xmin * scale;
        shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
    }
    else
    {
        scale = scale_y;
        shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
        shift_y = 1.0 - ymin * scale;
    }
    file << shift_x << " " << shift_y << " translate" << endl;
    #pragma endregion

    #pragma region Definitions of Shapes:

  file << "/circ    % stack: x y r" << endl;
    file << "{0 360 arc" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
  file << " stroke" << endl;
    file << "} def" << endl;

  file << "/fcirc    % stack: x y r" << endl;
    file << "{0 360 arc" << endl;
    file << " closepath" << endl;
    file << " gsave" << endl;
    file << " 0.6 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.01 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg2      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/quad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

  file << "/quad2      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
  file << " gsave" << endl;
    file << " 0.9 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;
    #pragma endregion

    double s(_dm * 0.05);

  for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
  {
    size_t jcell = icell + _neighbors[ii];

    if (jcell >= _num_cells) continue;

    if (!_bad_cells[jcell]) continue;

    if (_cell_points[jcell] == 0) continue; // _bad cells due to boundaries

    // retrieve the point associated with that cell
    double xj = _cell_points[jcell]->x;
    double yj = _cell_points[jcell]->y;

    file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
    file << "fcirc"     << endl;
  }

  if (kcell > 0)
  {
    for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
    {
      if (int(kcell) < _neighbors[ii]) continue;
      size_t jcell = kcell + _neighbors[ii];

      if (jcell >= _num_cells) continue;

      if (!_bad_cells[jcell]) continue;

      if (_cell_points[jcell] == 0) continue; // _bad cells due to boundaries

      // retrieve the point associated with that cell
      double xj = _cell_points[jcell]->x;
      double yj = _cell_points[jcell]->y;

      file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
      file << "fcirc"     << endl;
    }
  }

  if (pp != 0)
  {
    PolyPoint* pst = pp;
    while (true)
    {
      // plot polyline
      file << pp->x * scale << "  " << pp->y * scale << "  ";
      file << pp->next->x * scale << "  " << pp->next->y * scale << "  ";
      file << "seg"     << endl;
      pp = pp->next;
      if (pp == pst) break;
    }
  }

  for (size_t ipoly = 0; ipoly < _num_polys; ipoly++)
  {
    if (_polys[ipoly] == 0) break;
    PolyPoint* pst = _polys[ipoly];
    while (true)
    {
      // plot polyline
      file << _polys[ipoly]->x * scale << "  " << _polys[ipoly]->y * scale << "  ";
      file << _polys[ipoly]->next->x * scale << "  " << _polys[ipoly]->next->y * scale << "  ";
      file << "seg"     << endl;
      _polys[ipoly] = _polys[ipoly]->next;
      if (_polys[ipoly] == pst) break;
    }
  }

  // plot grid
  for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
  {
    size_t jcell = icell + _neighbors[ii];

    if (jcell >= _num_cells) continue;

    get_cell_indices(jcell, i, j);

    double xo = _bg_xo + i * _bg_s;
    double yo = _bg_yo + j * _bg_s;
    double xn = xo + _bg_s;
    double yn = yo + _bg_s;

    // plot cell
    file << xo * scale << "  " << yo * scale << "  ";
    file << xn * scale << "  " << yo * scale << "  ";
    file << xn * scale << "  " << yn * scale << "  ";
    file << xo * scale << "  " << yn * scale << "  ";
    file << "quad"     << endl;
  }

  for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
  {
    size_t jcell = icell + _neighbors[ii];

    if (jcell >= _num_cells) continue;

    if (!_bad_cells[jcell]) continue;

    if (_cell_points[jcell] == 0) continue; // _bad cells due to boundaries

    // retrieve the point associated with that cell
    double xj = _cell_points[jcell]->x;
    double yj = _cell_points[jcell]->y;

    file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
    file << "circ"     << endl;

    // plot vertex
        file << (xj - s) * scale << "  " << (yj - s) * scale << "  ";
        file << (xj + s) * scale << "  " << (yj + s) * scale << "  ";
        file << "seg2"     << endl;

        file << (xj + s) * scale << "  " << (yj - s) * scale << "  ";
        file << (xj - s) * scale << "  " << (yj + s) * scale << "  ";
        file << "seg2"     << endl;

  }

  if (kcell > 0)
  {
    for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
    {
      if (int(kcell) < _neighbors[ii]) continue;
      size_t jcell = kcell + _neighbors[ii];

      if (jcell >= _num_cells) continue;

      if (!_bad_cells[jcell]) continue;

      if (_cell_points[jcell] == 0) continue; // _bad cells due to boundaries

      // retrieve the point associated with that cell
      double xj = _cell_points[jcell]->x;
      double yj = _cell_points[jcell]->y;

      file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
      file << "circ"     << endl;

      // plot vertex
      file << (xj - s) * scale << "  " << (yj - s) * scale << "  ";
      file << (xj + s) * scale << "  " << (yj + s) * scale << "  ";
      file << "seg2"     << endl;

      file << (xj + s) * scale << "  " << (yj - s) * scale << "  ";
      file << (xj - s) * scale << "  " << (yj + s) * scale << "  ";
      file << "seg2"     << endl;
    }
  }

  if (_blue_cells[icell])
  {
    if (_cell_first_edge_map.find(icell) != _cell_first_edge_map.end())
    {
      size_t ied = _cell_first_edge_map[icell];
      file << _edges_x1[ied] * scale << "  " << _edges_y1[ied] * scale << "  ";
      file << _edges_x2[ied] * scale << "  " << _edges_y2[ied] * scale << "  ";
      file << "seg"     << endl;
    }
    if (_cell_second_edge_map.find(icell) != _cell_second_edge_map.end())
    {
      size_t ied = _cell_second_edge_map[icell];
      file << _edges_x1[ied] * scale << "  " << _edges_y1[ied] * scale << "  ";
      file << _edges_x2[ied] * scale << "  " << _edges_y2[ied] * scale << "  ";
      file << "seg"     << endl;
    }
  }
    file << "showpage" << endl;
    return 0;
    #pragma endregion
}

int MeshingGenie_2d::plot_point(std::string file_name, size_t icell, LoopPoint* po)
{
    #pragma region Plot Background grid:
  fstream file(file_name.c_str(), ios::out);
    file << "%!PS-Adobe-3.0" << endl;
    file << "72 72 scale     % one unit = one inch" << endl;

    #pragma region Retrieve bounding box, scale, and translate:
  size_t i, j;
  get_cell_indices(icell, i, j);

    double xmin(_bg_xo + i * _bg_s - 4 * _bg_s), ymin(_bg_yo + j * _bg_s - 4 * _bg_s), Lx(8 * _bg_s), Ly(8 * _bg_s);

    double scale_x, scale_y, scale;
    double shift_x, shift_y;

    scale_x = 6.5 / Lx;
    scale_y = 9.0 / Ly;

    if (scale_x < scale_y)
    {
        scale = scale_x;
        shift_x = 1.0 - xmin * scale;
        shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
    }
    else
    {
        scale = scale_y;
        shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
        shift_y = 1.0 - ymin * scale;
    }
    file << shift_x << " " << shift_y << " translate" << endl;
    #pragma endregion

    #pragma region Definitions of Shapes:

  file << "/circ    % stack: x y r" << endl;
    file << "{0 360 arc" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
  file << " stroke" << endl;
    file << "} def" << endl;

  file << "/fcirc    % stack: x y r" << endl;
    file << "{0 360 arc" << endl;
    file << " closepath" << endl;
    file << " gsave" << endl;
    file << " 0.6 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.01 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg2      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/quad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

  file << "/quad2      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
  file << " gsave" << endl;
    file << " 0.9 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;
    #pragma endregion

    double s(_dm * 0.05);
  double xj, yj;

  std::set<size_t> edges;
  get_neighbor_points(icell);
  for (size_t ii = 0; ii < _num_neighbor_points_dynamic; ii++)
  {
    if (ii == 0)
    {
      xj = _cell_points[icell]->x;
      yj = _cell_points[icell]->y;

      file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
      file << "fcirc"     << endl;
    }

    if (_neighbor_points[ii] == 0) break;

    // retrieve the point associated with that cell
    xj = _neighbor_points[ii]->x;
    yj = _neighbor_points[ii]->y;

    file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
    file << "fcirc"     << endl;
  }

  // plot grid
  for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
  {
    size_t jcell = icell + _neighbors[ii];

    if (jcell >= _num_cells) continue;

    get_cell_indices(jcell, i, j);

    double xo = _bg_xo + i * _bg_s;
    double yo = _bg_yo + j * _bg_s;
    double xn = xo + _bg_s;
    double yn = yo + _bg_s;

    // plot cell
    file << xo * scale << "  " << yo * scale << "  ";
    file << xn * scale << "  " << yo * scale << "  ";
    file << xn * scale << "  " << yn * scale << "  ";
    file << xo * scale << "  " << yn * scale << "  ";
    file << "quad"     << endl;
  }

  LoopPoint* pst = po;
  if (pst != 0)
  {
    while (pst->next != 0 )
    {
      file << (pst->pp->x) * scale << "  " << (pst->pp->y) * scale << "  ";
      file << (pst->next->pp->x) * scale << "  " << (pst->next->pp->y) * scale << "  ";
      file << "seg"     << endl;
      pst = pst->next;
      if (po == pst) break;
    }
  }

  for (size_t ii = 0; ii < _num_neighbor_points; ii++)
  {
    if (ii == 0)
    {
      xj = _cell_points[icell]->x;
      yj = _cell_points[icell]->y;

      file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
      file << "circ"     << endl;

      // plot vertex
      file << (xj - s) * scale << "  " << (yj - s) * scale << "  ";
      file << (xj + s) * scale << "  " << (yj + s) * scale << "  ";
      file << "seg2"     << endl;

      file << (xj + s) * scale << "  " << (yj - s) * scale << "  ";
      file << (xj - s) * scale << "  " << (yj + s) * scale << "  ";
      file << "seg2"     << endl;
    }

    if (_neighbor_points[ii] == 0) break; // _bad cells due to boundaries

    // retrieve the point associated with that cell
    xj = _neighbor_points[ii]->x;
    yj = _neighbor_points[ii]->y;

    file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
    file << "circ"     << endl;

    // plot vertex
        file << (xj - s) * scale << "  " << (yj - s) * scale << "  ";
        file << (xj + s) * scale << "  " << (yj + s) * scale << "  ";
        file << "seg2"     << endl;

        file << (xj + s) * scale << "  " << (yj - s) * scale << "  ";
        file << (xj - s) * scale << "  " << (yj + s) * scale << "  ";
        file << "seg2"     << endl;
  }
    file << "showpage" << endl;
    return 0;
    #pragma endregion
}

int MeshingGenie_2d::plot_star(std::string file_name, double xo, double yo, size_t icell, Point** neighbors, bool* invalid, int num)
{
    #pragma region Plot Star:
  fstream file(file_name.c_str(), ios::out);
    file << "%!PS-Adobe-3.0" << endl;
    file << "72 72 scale     % one unit = one inch" << endl;

    #pragma region Retrieve bounding box, scale, and translate:
  size_t i, j;
  get_cell_indices(icell, i, j);

    double xmin(_bg_xo + i * _bg_s - 4 * _bg_s), ymin(_bg_yo + j * _bg_s - 4 * _bg_s), Lx(8 * _bg_s), Ly(8 * _bg_s);

    double scale_x, scale_y, scale;
    double shift_x, shift_y;

    scale_x = 6.5 / Lx;
    scale_y = 9.0 / Ly;

    if (scale_x < scale_y)
    {
        scale = scale_x;
        shift_x = 1.0 - xmin * scale;
        shift_y = 0.5 * (11.0 - Ly * scale) - ymin * scale;
    }
    else
    {
        scale = scale_y;
        shift_x = 0.5 * (8.5 - Lx * scale) - xmin * scale;
        shift_y = 1.0 - ymin * scale;
    }
    file << shift_x << " " << shift_y << " translate" << endl;
    #pragma endregion

    #pragma region Definitions of Shapes:

  file << "/circ    % stack: x y r" << endl;
    file << "{0 360 arc" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
  file << " stroke" << endl;
    file << "} def" << endl;

  file << "/fcirc    % stack: x y r" << endl;
    file << "{0 360 arc" << endl;
    file << " closepath" << endl;
    file << " gsave" << endl;
    file << " 0.6 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.01 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/seg2      % stack: x1 y1 x2 y2" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

    file << "/quad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;

  file << "/quad2      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl;
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
  file << " gsave" << endl;
    file << " 0.9 setgray fill" << endl;
  file << " grestore" << endl;
    file << " 0.002 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;
    #pragma endregion

    double s(_dm * 0.05);
  double xj, yj;

  std::set<size_t> edges;
  get_neighbor_points(icell);
  for (int ii = 0; ii < num; ii++)
  {
    if (invalid[ii]) continue;

    if (_neighbor_points[ii] == 0) break;

    // retrieve the point associated with that cell
    xj = neighbors[ii]->x;
    yj = neighbors[ii]->y;

    file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
    file << "fcirc"     << endl;
  }

  // plot grid
  for (size_t ii = 0; ii < _num_neighbor_cells; ii++)
  {
    size_t jcell = icell + _neighbors[ii];

    if (jcell >= _num_cells) continue;

    get_cell_indices(jcell, i, j);

    double xo = _bg_xo + i * _bg_s;
    double yo = _bg_yo + j * _bg_s;
    double xn = xo + _bg_s;
    double yn = yo + _bg_s;

    // plot cell
    file << xo * scale << "  " << yo * scale << "  ";
    file << xn * scale << "  " << yo * scale << "  ";
    file << xn * scale << "  " << yn * scale << "  ";
    file << xo * scale << "  " << yn * scale << "  ";
    file << "quad"     << endl;
  }

  for (int ii = 0; ii < num; ii++)
  {
    if (invalid[ii]) continue;

    if (_neighbor_points[ii] == 0) break;

    int iip(ii);
    iip++; if (iip == num) iip = 0;

    while (invalid[iip])
    {
      iip++; if (iip == num) iip = 0;
    }

    double A = area_triangle(xo, yo, neighbors[ii]->x, neighbors[ii]->y, neighbors[iip]->x, neighbors[iip]->y);
    if (A > 0.0)
    {
      file << (neighbors[ii]->x) * scale << "  " << (neighbors[ii]->y) * scale << "  ";
      file << (neighbors[iip]->x) * scale << "  " << (neighbors[iip]->y) * scale << "  ";
      file << "seg"     << endl;
    }

    file << neighbors[ii]->x * scale << "  " << neighbors[ii]->y * scale << "  " << _dm * scale << "  ";
    file << "circ"     << endl;
  }

  file << xo * scale << "  " << yo * scale << "  " << _dm * scale << "  ";
  file << "circ"     << endl;

  for (int ii = 0; ii < num; ii++)
  {
    if (_neighbor_points[ii] == 0) break;

    // retrieve the point associated with that cell
    xj = neighbors[ii]->x;
    yj = neighbors[ii]->y;

    file << xj * scale << "  " << yj * scale << "  " << _dm * scale << "  ";
    file << "circ"     << endl;
  }

  file << (xo - s) * scale << "  " << (yo - s) * scale << "  ";
  file << (xo + s) * scale << "  " << (yo + s) * scale << "  ";
    file << "seg2"     << endl;

  file << (xo + s) * scale << "  " << (yo - s) * scale << "  ";
    file << (xo - s) * scale << "  " << (yo + s) * scale << "  ";
  file << "seg2"     << endl;

  for (int ii = 0; ii < num; ii++)
  {
    if (_neighbor_points[ii] == 0) break;

    // retrieve the point associated with that cell
    xj = neighbors[ii]->x;
    yj = neighbors[ii]->y;

    file << (xj - s) * scale << "  " << (yj - s) * scale << "  ";
    file << (xj + s) * scale << "  " << (yj + s) * scale << "  ";
    file << "seg2"     << endl;

    file << (xj + s) * scale << "  " << (yj - s) * scale << "  ";
    file << (xj - s) * scale << "  " << (yj + s) * scale << "  ";
    file << "seg2"     << endl;
  }

    file << "showpage" << endl;
    return 0;
    #pragma endregion
}

///////////////////////////////////////////////////////////////////////////////////////////////
///
/// Connectivity: based on the new approach connect to all then clean up
///
///////////////////////////////////////////////////////////////////////////////////////////////

int MeshingGenie_2d::CDT(bool generate_Voronoi_Cells)
{
  #pragma region CDT:

  _num_oriented_edges = _num_boundary_edges + _num_holes_edges;

  //double xp(-1.0), yp(2.0); // for debuging

  #pragma region Angles data:
  /*
  double min_ang(180.0), max_ang(0.0);
  std::vector<size_t> angles(60);
  std::vector<size_t> bangles(60);
  std::vector<size_t> edlen(60);
  std::vector<size_t> bedlen(60);
  std::vector<double> edge_range(60);
  std::vector<double> ang_range(60);
  double edL(0.8), ang(25.0);
  for (int ia = 0; ia < 15; ia++)
  {
    angles[ia] = 0;
    bangles[ia] = 0;
    edlen[ia] = 0;
    bedlen[ia] = 0;
    edge_range[ia] = edL;
    ang_range[ia] = ang;

    if (ang < 40.0) ang += 5.0;
    else if (ang < 100.0) ang +=15.0;
    else ang += 5.0;

    edL += 0.1;
  }
  */
  #pragma endregion

  //if (_save_point_cloud) save_point_cloud("point.dat");

  print_message("Generating Tessellation");

  clock_t start_time, end_time; double cpu_time;
    start_time = clock();

  refill_neighbors_data();

  int k(0), k1(0), k2(0), k3(0), k4(0);
  Point** qi_neighbors = new Point*[400];
  bool* qi_invalid_neighbors = new bool[400];
  Point** qi_neighbors_1 = new Point*[100];
  Point** qi_neighbors_2 = new Point*[100];
  Point** qi_neighbors_3 = new Point*[100];
  Point** qi_neighbors_4 = new Point*[100];
  double* qi_slope_1 = new double[100]; // -45 to 45
  double* qi_slope_2 = new double[100]; // 45 to 135
  double* qi_slope_3 = new double[100]; // 135 to 225
  double* qi_slope_4 = new double[100]; // 225 to 315

  int ked;
  double* xed_1 = new double[100];
  double* yed_1 = new double[100];
  double* xed_2 = new double[100];
  double* yed_2 = new double[100];
  Point** ped   = new Point*[100];

  double max_dist_sq = 4.10 * _dm_squared;

  int num_rejected(0);
  double xo, yo;
  bool first(true);
  Point* qi; Point* qj; Point* tmpq;
  Point* q1;Point* q2; Point* q3;
  for (size_t icell = 0; icell < _num_cells; icell++)
  {
    if (_cell_points[icell] == 0) continue;

    if (first) qi = _cell_points[icell];

    xo = qi->x; yo = qi->y;

    bool plot(false);
    get_neighbor_points(icell);

    k = 0; k1 = 0; k2 = 0; k3 = 0; k4 = 0;

    for (size_t ii = 0; ii < _num_neighbor_points_dynamic; ii++)
    {
      #pragma region Retrieve neighbors in four zones:
      if (_neighbor_points[ii] == 0) break;
      if (_neighbor_points[ii] == qi) continue;

      qj = _neighbor_points[ii];
      if (notconnected_points(qi, qj)) continue;

      double dd = distance_squared(qi->x - qj->x, qi->y - qj->y);
      if (dd > max_dist_sq)
      {
        continue; // points are too far from each other
      }

      double xj(qj->x), yj(qj->y);

      double dx(xj - xo), dy(yj - yo);

      if (dx >= 0 && dx >= fabs(dy) - 1E-10)
      {
        qi_neighbors_1[k1] = qj;
        qi_slope_1[k1] = dy / dx;
        int jj(k1), jjm;
        while (jj > 0)
        {
          jjm = jj - 1;
          if (qi_slope_1[jj] < qi_slope_1[jj - 1])
          {
            double tmp = qi_slope_1[jjm]; qi_slope_1[jjm] = qi_slope_1[jj]; qi_slope_1[jj] = tmp;
            tmpq = qi_neighbors_1[jjm]; qi_neighbors_1[jjm] = qi_neighbors_1[jj]; qi_neighbors_1[jj] = tmpq;
            jj--;
          }
          else break;
        }
        k1++;
      }
      else if (dy >= 0 && dy >= fabs(dx) - 1E-10)
      {
        qi_neighbors_2[k2] = qj;
        qi_slope_2[k2] = dx / dy;
        int jj(k2), jjm;
        while (jj > 0)
        {
          jjm = jj - 1;
          if (qi_slope_2[jj] > qi_slope_2[jj - 1])
          {
            double tmp = qi_slope_2[jjm]; qi_slope_2[jjm] = qi_slope_2[jj]; qi_slope_2[jj] = tmp;
            tmpq = qi_neighbors_2[jjm]; qi_neighbors_2[jjm] = qi_neighbors_2[jj]; qi_neighbors_2[jj] = tmpq;
            jj--;
          }
          else break;
        }
        k2++;
      }
      else if (dx < 0 && -dx >= fabs(dy) - 1E-10)
      {
        qi_neighbors_3[k3] = qj;
        qi_slope_3[k3] = dy / dx;
        int jj(k3), jjm;
        while (jj > 0)
        {
          jjm = jj - 1;
          if (qi_slope_3[jj] < qi_slope_3[jj - 1])
          {
            double tmp = qi_slope_3[jjm]; qi_slope_3[jjm] = qi_slope_3[jj]; qi_slope_3[jj] = tmp;
            tmpq = qi_neighbors_3[jjm]; qi_neighbors_3[jjm] = qi_neighbors_3[jj]; qi_neighbors_3[jj] = tmpq;
            jj--;
          }
          else break;
        }
        k3++;
      }
      else if (dy < 0 && -dy >= fabs(dx) - 1E-10)
      {
        qi_neighbors_4[k4] = qj;
        qi_slope_4[k4] = dx / dy;
        int jj(k4), jjm;
        while (jj > 0)
        {
          jjm = jj - 1;
          if (qi_slope_4[jj] > qi_slope_4[jj - 1])
          {
            double tmp = qi_slope_4[jjm]; qi_slope_4[jjm] = qi_slope_4[jj]; qi_slope_4[jj] = tmp;
            tmpq = qi_neighbors_4[jjm]; qi_neighbors_4[jjm] = qi_neighbors_4[jj]; qi_neighbors_4[jj] = tmpq;
            jj--;
          }
          else break;
        }
        k4++;
      }
      #pragma endregion
    }

    #pragma region sort Points based on slope:
    for (int i = 0; i < k1; i++)
    {
      if (i == 0) {qi_neighbors[k] = qi_neighbors_1[i]; k++;}
      else
      {
        // check if same slope as previous point
        if (fabs(qi_slope_1[i] - qi_slope_1[i - 1]) < 1E-10)
        {
          // pick closest point
          int km = k - 1;
          double dqk = distance_squared(qi_neighbors[km]->x - xo, qi_neighbors[km]->y - yo);
          double dqi = distance_squared(qi_neighbors_1[i]->x - xo, qi_neighbors_1[i]->y - yo);
          if (dqi < dqk) qi_neighbors[k-1] = qi_neighbors_1[i];
        }
        else {qi_neighbors[k] = qi_neighbors_1[i]; k++;}
      }
    }

    for (int i = 0; i < k2; i++)
    {
      if (i == 0)
      {
        // check if same slope as previous point
        if (k1 > 0 && fabs(qi_slope_2[i] - 1) < 1E-10 && fabs(qi_slope_1[k1 - 1] - 1) < 1E-10)
        {
          // pick closest point
          int km = k - 1;
          double dqk = distance_squared(qi_neighbors[km]->x - xo, qi_neighbors[km]->y - yo);
          double dqi = distance_squared(qi_neighbors_2[i]->x - xo, qi_neighbors_2[i]->y - yo);
          if (dqi < dqk) qi_neighbors[k-1] = qi_neighbors_2[i];
        }
        else {qi_neighbors[k] = qi_neighbors_2[i]; k++;}
      }
      else
      {
        // check if same slope as previous point
        if (fabs(qi_slope_2[i] - qi_slope_2[i - 1]) < 1E-10)
        {
          // pick closest point
          int km = k - 1;
          double dqk = distance_squared(qi_neighbors[km]->x - xo, qi_neighbors[km]->y - yo);
          double dqi = distance_squared(qi_neighbors_2[i]->x - xo, qi_neighbors_2[i]->y - yo);
          if (dqi < dqk) qi_neighbors[k-1] = qi_neighbors_2[i];
        }
        else {qi_neighbors[k] = qi_neighbors_2[i]; k++;}
      }
    }

    for (int i = 0; i < k3; i++)
    {
      if (i == 0)
      {
        // check if same slope as previous point
        if (k2 > 0 && fabs(qi_slope_3[i] + 1) < 1E-10 && fabs(qi_slope_2[k2 - 1] + 1) < 1E-10)
        {
          // pick closest point
          int km = k - 1;
          double dqk = distance_squared(qi_neighbors[km]->x - xo, qi_neighbors[km]->y - yo);
          double dqi = distance_squared(qi_neighbors_3[i]->x - xo, qi_neighbors_3[i]->y - yo);
          if (dqi < dqk) qi_neighbors[k-1] = qi_neighbors_3[i];
        }
        else {qi_neighbors[k] = qi_neighbors_3[i]; k++;}
      }
      else
      {
        // check if same slope as previous point
        if (fabs(qi_slope_3[i] - qi_slope_3[i - 1]) < 1E-10)
        {
          // pick closest point
          int km = k - 1;
          double dqk = distance_squared(qi_neighbors[km]->x - xo, qi_neighbors[km]->y - yo);
          double dqi = distance_squared(qi_neighbors_3[i]->x - xo, qi_neighbors_3[i]->y - yo);
          if (dqi < dqk) qi_neighbors[k-1] = qi_neighbors_3[i];
        }
        else {qi_neighbors[k] = qi_neighbors_3[i]; k++;}
      }
    }

    for (int i = 0; i < k4; i++)
    {
      if (i == 0)
      {
        // check if same slope as previous point
        if (k3 > 0 && fabs(qi_slope_4[i] - 1) < 1E-10 && fabs(qi_slope_3[k3 - 1] - 1) < 1E-10)
        {
          // pick closest point
          int km = k - 1;
          double dqk = distance_squared(qi_neighbors[km]->x - xo, qi_neighbors[km]->y - yo);
          double dqi = distance_squared(qi_neighbors_4[i]->x - xo, qi_neighbors_4[i]->y - yo);
          if (dqi < dqk) qi_neighbors[k-1] = qi_neighbors_4[i];
        }
        else {qi_neighbors[k] = qi_neighbors_4[i]; k++;}
      }
      else
      {
        // check if same slope as previous point
        if (fabs(qi_slope_4[i] - qi_slope_4[i - 1]) < 1E-10)
        {
          // pick closest point
          int km = k - 1;
          double dqk = distance_squared(qi_neighbors[km]->x - xo, qi_neighbors[km]->y - yo);
          double dqi = distance_squared(qi_neighbors_4[i]->x - xo, qi_neighbors_4[i]->y - yo);
          if (dqi < dqk) qi_neighbors[k-1] = qi_neighbors_4[i];
        }
        else {qi_neighbors[k] = qi_neighbors_4[i]; k++;}
      }
    }

    for (int i = 0; i < k; i++) qi_invalid_neighbors[i] = false;

    #pragma endregion


    //plot_star("star.ps", xo, yo, icell, qi_neighbors, qi_invalid_neighbors, k);

    bool forward(false);
    if (qi->on_boundary)
    {
      #pragma region Remove edges in the exterior of the domain:
      for (int i = 0; i < k; i++)
      {
        if (qi_invalid_neighbors[i]) continue;

        qj = qi_neighbors[i];

        if (oriented_edge_exists(qi, qj, forward) && forward)
        {
          int im(i);
          while (true)
          {
            im--; if (im < 0) im = k - 1;
            if (qi_invalid_neighbors[im]) continue;
            qj = qi_neighbors[im];
            if (oriented_edge_exists(qi, qj, forward) && !forward) break;
            qi_invalid_neighbors[im] = true;
          }
        }
      }
      #pragma endregion
    }

    int io(0);
    for (int i = 0; i < k; i++)
    {
      #pragma region Clean up due to respect Constraint Edges:
      if (qi_invalid_neighbors[i]) continue;

      int im(i);
      im--; if (im < 0) im = k - 1;

      int j = i;
      j++; if (j == k) j = 0;
      j++; if (j == k) j = 0;
      qj = qi_neighbors[i];
      while (j != im)
      {
        if (qi_invalid_neighbors[j])
        {
          j++; if (j == k) j = 0;
          continue;
        }

        tmpq = qi_neighbors[j];
        if (connected_points(qj, tmpq))
        {
          double A = area_parallelogram(xo, yo, qj->x, qj->y, tmpq->x, tmpq->y);
          if (A > 0.0)
          {
            // delete points from i + 1 to j - 1 (move forward)
            int ip(i);
            ip++; if (ip == k) ip = 0;
            int jj(ip);
            while (jj != j)
            {
              A = area_parallelogram(qj->x, qj->y, qi_neighbors[jj]->x, qi_neighbors[jj]->y, tmpq->x, tmpq->y);
              if (A > 0.0)
              {
                qi_invalid_neighbors[jj] = true;
              }
              jj++; if (jj == k) jj = 0;
            }
          }
          else
          {
            // delete points from i - 1 to j + 1 (move backward)
            int im(i);
            im--; if (im < 0) im = k - 1;
            int jj(im);
            while (jj != j)
            {
              A = area_parallelogram(qj->x, qj->y, qi_neighbors[im]->x, qi_neighbors[im]->y, tmpq->x, tmpq->y);
              if (A < 0.0)
              {
                qi_invalid_neighbors[jj] = true;
              }
              jj--; if (jj < 0) jj = k - 1;
            }
          }
        }
        j++; if (j == k) j = 0;
        continue;
      }
      #pragma endregion
    }

    PointEdge* ed_o;
    PointEdge* ed_i;
    ked = 0;
    for (int i = 0; i < k; i++)
    {
      #pragma region get all boundary edges associated with qj such that the other end does not exist in loop:
      if (qi_invalid_neighbors[i]) continue;

      qj = qi_neighbors[i];

      if (qj->on_boundary)
      {
        ed_i = ed_o = qj->edges;
        while (ed_i != 0)
        {
          tmpq = ed_i->edge_other_end;
          if (tmpq == 0 || tmpq == qi)
          {
            ed_i = ed_i->next;
            if (ed_i == ed_o) break;
            continue;
          }

          bool found(false);
          int j(i);
          while (true)
          {
            j++; if (j == k) j = 0;
            if (j == i) break;
            if (qi_invalid_neighbors[j]) continue;
            if (qi_neighbors[j] == tmpq) {found = true; break;}
          }
          if (!found)
          {
            xed_1[ked] = qj->x; yed_1[ked] = qj->y;
            xed_2[ked] = tmpq->x; yed_2[ked] = tmpq->y;
            ped[ked] = qj;
            ked++;
          }

          ed_i = ed_i->next;
          if (ed_i == ed_o) break;
        }
      }
      #pragma endregion
    }


    for (int ied = 0; ied < ked; ied++)
    {
      #pragma region Clean up to respect Constrained boundary edges ouside star:
      for (int i = 0; i < k; i++)
      {
        if (qi_invalid_neighbors[i]) continue;
        qj = qi_neighbors[i];

        if (qj == ped[ied]) continue;

        double A1 = area_parallelogram(xed_1[ied], yed_1[ied], xed_2[ied], yed_2[ied], xo, yo);
        double A2 = area_parallelogram(xed_1[ied], yed_1[ied], xed_2[ied], yed_2[ied], qj->x, qj->y);
        if (A1 > 1E-10 && A2 > 1E-10) continue;
        if (A1 < -1E-10 && A2 < -1E-10) continue;
        if (fabs(A1) < 1E-10 && fabs(A2) < 1E-10) continue;
        A1 = area_parallelogram(xed_1[ied], yed_1[ied], xo, yo, qj->x, qj->y);
        A2 = area_parallelogram(xed_2[ied], yed_2[ied], xo, yo, qj->x, qj->y);
        if (A1 > 1E-10 && A2 > 1E-10) continue;
        if (A1 < -1E-10 && A2 < -1E-10) continue;
        if (fabs(A1) < 1E-10 && fabs(A2) < 1E-10) continue;
        qi_invalid_neighbors[i] = true;
      }
      #pragma endregion
    }


    int i(-1); bool second_loop(false);
    while (true)
    {
      #pragma region Clean up due to Delaunay Principal:
      i++;
      if (i == k)
      {
        i = 0; second_loop = true;
      }

      if (qi_invalid_neighbors[i]) continue;

      q1 = qi_neighbors[i];
      bool connected = connected_points(qi, q1);

      if (connected && second_loop) break;
      else if (connected) continue;

      int ip(i);
      while (true)
      {
        ip++; if (ip == k) ip = 0;
        if (!qi_invalid_neighbors[ip]) break;
      }
      q2 = qi_neighbors[ip];

      int im(i);
      while (true)
      {
        im--; if (im < 0) im = k - 1;
        if (!qi_invalid_neighbors[im]) break;
      }
      q3 = qi_neighbors[im];

      // check if q1 lies in xo, yo, q2, q3 circum circle
      if (!incircumcircle(xo, yo, q2->x, q2->y, q3->x, q3->y, q1->x, q1->y))
      {
        qi_invalid_neighbors[i] = true;
        if (im < i) {i = im - 1; continue;}
      }
      else if (second_loop) break;
      #pragma endregion
    }


    for (int i = 0; i < k; i++)
    {
      #pragma region Add star constarints:
      if (qi_invalid_neighbors[i]) continue;

      int j(i);
      while (true)
      {
        j++; if (j == k) j = 0;
        if (!qi_invalid_neighbors[j]) break;
      }

      q1 = qi_neighbors[i];

      q2 = qi_neighbors[j];

      if (q1->done || q2->done) continue;

      if (qi->on_boundary)
      {
        double A = area_parallelogram(xo, yo, q1->x, q1->y, q2->x, q2->y);
        if (A > 1E-10) add_constraint_edge(q1, q2);
      }
      else add_constraint_edge(q1, q2);
      #pragma endregion
    }


    size_t edge_id; size_t num_edges(0);
    PointEdge* edges(0); PointEdge* newedge(0);  PointEdge* firstedge(0); PointEdge* boundaryedge(0);
    for (int i = 0; i < k; i++)
    {
      #pragma region Stores ordered connectivity:
      if (qi_invalid_neighbors[i]) continue;
      q1 = qi_neighbors[i];

      newedge = new PointEdge();
      newedge->next = 0;
      newedge->edge_other_end = q1;
      newedge->edge_id = 0; num_edges++;
      if (is_boundary_edge(qi, q1, edge_id))
      {
        newedge->edge_id = edge_id;
        bool forward_direction = bool((edge_id - 1) & 1);
        if (forward_direction || boundaryedge == 0) boundaryedge = newedge;
      }
      else if (!q1->done)
      {
        // Add qi to q1 edges
        PointEdge* ed = q1->edges; bool found(false);
        while (true)
        {
          if (ed->edge_other_end == qi) {found = true; break;} // edge already exists
          if (ed->next == 0 || ed->next == q1->edges) break;
          ed = ed->next;
        }
        if (!found)
        {
          ed->next = new PointEdge();
          ed->next->edge_other_end = qi;
          ed->next->edge_id = 0; // unclassified edge
          ed->next->next = 0;
        }
      }

      if (edges == 0) {firstedge = newedge;}
      else {edges->next = newedge;}

      edges = newedge;
      #pragma endregion
    }

    newedge->next = firstedge;
    if (boundaryedge != 0 && firstedge != boundaryedge)
    {
      if (num_edges == 2) firstedge->next = 0;
      firstedge = boundaryedge;
    }
    else if (num_edges == 2 && firstedge == boundaryedge) newedge->next = 0;


    // delete old edges
    edges = qi->edges; PointEdge* nextedge;
    while (edges != 0)
    {
      nextedge = edges->next;
      delete edges;
      if (nextedge == 0 || nextedge ==  qi->edges) break;
      edges = nextedge;
    }
    qi->edges = firstedge;
    qi->done = true; qi->num_edges = num_edges;

    if (qi->next != 0)
    {
      qi = qi->next; icell--; // stay in the same cell
      first = false;
    }
    else first = true;

    #pragma region Angle Caluculations:
    /* calculate angles and edge lengths
    edges = qi->edges;
    while (edges != 0)
    {
      nextedge = edges->next;
      double x1 = edges->edge_other_end->x; double y1 = edges->edge_other_end->y;

      if (!edges->edge_other_end->done)
      {
        double tmpL = distance_squared(x1 - xo, y1 - yo);
        tmpL = sqrt(tmpL);
        for (int ia = 0; ia < 60; ia++)
        {
          if (tmpL < edge_range[ia] * _dm + 1E-10)
          {
            size_t edge_id;
            if (is_boundary_edge(qi, edges->edge_other_end, edge_id))
            {
              bedlen[ia]++;
            }
            else edlen[ia]++;

            break;
          }
        }
      }

      if (nextedge == 0 || nextedge ==  qi->edges) break;

      double x2 = nextedge->edge_other_end->x; double y2 = nextedge->edge_other_end->y;
      double dot = (x1 - xo) * (x2 - xo) + (y1 - yo) * (y2 - yo);
      double L1 = sqrt(distance_squared(x1-xo, y1-yo));
      double L2 = sqrt(distance_squared(x2-xo, y2-yo));
      double L3 = sqrt(distance_squared(x2-x1, y2-y1));
      dot /= (L1 * L2);
      double angle = acos(dot)*180/3.1415926;
      if (angle > max_ang) max_ang = angle;
      if (angle < min_ang) min_ang = angle;

      for (int ia = 0; ia < 60; ia++)
      {
        if (angle < ang_range[ia] + 0.001)
        {
          size_t edge_id;
          if (is_boundary_edge(qi, edges->edge_other_end, edge_id) ||
            is_boundary_edge(qi, nextedge->edge_other_end, edge_id) ||
            is_boundary_edge(edges->edge_other_end, nextedge->edge_other_end, edge_id))
          {
            bangles[ia]++;
          }
          else angles[ia]++;

          break;
        }
      }

      edges = nextedge;
    }
    */
    #pragma endregion
  }


  end_time = clock();
    cpu_time = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    print_message_end("^ Execution Time = " + double_to_string(cpu_time) + " seconds.");

  #pragma region Saving Angles Data:
  /*
  std::cout<< " angle bounds " << min_ang << " , " << max_ang << std::endl;

  std::fstream file("edges.dat", std::ios::out);
  for (size_t i = 1; i < 15; i++)
  {
    file << edlen[i] << std::endl;
  }

  std::fstream file_bed("bedges.dat", std::ios::out);
  for (size_t i = 1; i < 15; i++)
  {
    file_bed << bedlen[i] << std::endl;
  }

  std::fstream file_("angles.dat", std::ios::out);
  for (size_t i = 1; i < 15; i++)
  {
    file_ << angles[i] << std::endl;
  }

  std::fstream file_ba("bangles.dat", std::ios::out);
  for (size_t i = 1; i < 15; i++)
  {
    file_ba << bangles[i] << std::endl;
  }
  */
  #pragma endregion

  return 0;
  #pragma endregion
}
