/*!
 * \file CLookUpTable.cpp
 * \brief Definition of the look-up table class members.
 * \author D. Mayer, T. Economon, N.Beishuizen, G. Tamanampudi
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/CLookUpTable.hpp"

CLookUpTable::CLookUpTable(string file_name_lut) {

  rank = SU2_MPI::GetRank();

  LoadTableRaw(file_name_lut);

  FindTableLimits();

  if (rank == MASTER_NODE)
    cout << "Detecting all unique edges and setting edge to triangle connectivity "
            "..." << endl;

  IdentifyUniqueEdges();

  if (rank == MASTER_NODE) cout << " done." << endl;

  PrintTableInfo();

  if (rank == MASTER_NODE)
    cout << "Building a trapezoidal map for the (progress variable, enthalpy) "
            "space ..." << endl;

  trap_map_prog_enth = CTrapezoidalMap(GetData("ProgVar"), GetData("Enthalpy"),
                                       GetEdges(), GetEdgeToTriangle());

  if (rank == MASTER_NODE) cout << " done." << endl;

  if (rank == MASTER_NODE) {
    cout << "Precomputing interpolation coefficients..." << endl;
  }

  ComputeInterpCoeffs();
  
  if (rank == MASTER_NODE) cout << " done." << endl;

  if (rank == MASTER_NODE) cout << "LUT fluid model ready for use" << endl;
  
}

void CLookUpTable::LoadTableRaw(string file_name_lut) {

  CFileReaderLUT file_reader;
  
  if (rank == MASTER_NODE)
  cout << "Loading look-up-table-file " << file_name_lut << " ..." << endl;

  file_reader.ReadRawDRG(file_name_lut);

  n_points       = file_reader.GetNPoints();
  n_triangles    = file_reader.GetNTriangles();
  n_variables    = file_reader.GetNVariables();
  type_lut       = file_reader.GetTypeLUT();
  version_lut    = file_reader.GetVersionLUT();
  version_reader = file_reader.GetVersionReader();

  names_var  = file_reader.GetNamesVar();
  table_data = file_reader.GetTableData();
  triangles  = file_reader.GetTriangles();

  if (rank == MASTER_NODE) cout << " done." << endl;
  
}

void CLookUpTable::FindTableLimits() {

  int ixEnth = GetIndexOfVar("Enthalpy");
  int ixProg = GetIndexOfVar("ProgVar");

  limits_table_enth[0] =
      *min_element(table_data.at(ixEnth).begin(), table_data.at(ixEnth).end());
  limits_table_enth[1] =
      *max_element(table_data.at(ixEnth).begin(), table_data.at(ixEnth).end());
  limits_table_prog[0] =
      *min_element(table_data.at(ixProg).begin(), table_data.at(ixProg).end());
  limits_table_prog[1] =
      *max_element(table_data.at(ixProg).begin(), table_data.at(ixProg).end());

}

void CLookUpTable::PrintTableInfo() {

  if (rank == MASTER_NODE) {
    cout << setfill(' ');
    cout << endl;
    cout << "+----------------------------------------------+"                         << endl;
    cout << "|           Look-Up-Table (LUT) info           |"                         << endl;
    cout << "+----------------------------------------------+"                         << endl;
    cout << "| Table type:"           << setw(33) << type_lut                  << " |" << endl;
    cout << "| Table version:"        << setw(30) << version_lut               << " |" << endl;
    cout << "| Table reader version:" << setw(23) << version_reader            << " |" << endl;
    cout << "| Number of variables:"  << setw(24) << n_variables               << " |" << endl;
    cout << "| Number of points:"     << setw(27) << n_points                  << " |" << endl;
    cout << "| Number of triangles:"  << setw(24) << n_triangles               << " |" << endl;
    cout << "| Number of edges:"      << setw(28) << edges.size()              << " |" << endl;
    cout << "+----------------------------------------------+"                         << endl;
    cout << "| Minimum enthalpy:" << setw(27) << limits_table_enth[0]          << " |" << endl;
    cout << "| Maximum enthalpy:" << setw(27) << limits_table_enth[1]          << " |" << endl;
    cout << "| Minimum progress variable:" << setw(18) << limits_table_prog[0] << " |" << endl;
    cout << "| Maximum progress variable:" << setw(18) << limits_table_prog[1] << " |" << endl;
    cout << "+----------------------------------------------+" << endl;
    cout << "| Variable names:                              |" << endl;
    cout << "|                                              |" << endl;

    for (unsigned long i_var = 0; i_var < names_var.size(); i_var++)
      cout << "| " << right << setw(3) << i_var << ": " << left << setw(39) 
      << names_var.at(i_var) << " |" << endl;

    cout << "+----------------------------------------------+" << endl;

    cout << endl;
  }
}

void CLookUpTable::IdentifyUniqueEdges() {
  
  /* We will fill these two data members with the point pair and adjacent
   elements, respectively, for each unique edge in the grid. */
  edges.clear();
  edge_to_triangle.clear();
  
  /* Loop through elements and store the element ID as
   a neighbor for each point in the element. */
  vector<vector<unsigned long> > neighborElemsOfPoint;
  neighborElemsOfPoint.resize(n_points);
  for (unsigned long iElem = 0; iElem < n_triangles; iElem++) {
    
    /* loop over 3 points per triangle */
    for (unsigned long iPoint = 0; iPoint < N_POINTS_TRIANGLE; iPoint++) {
      
      /* Get the global ID of the current point. */
      const unsigned long GlobalIndex = triangles.at(iElem).at(iPoint);
      
      /* Add the current element ID to the neighbor list for this point. */
      neighborElemsOfPoint[GlobalIndex].push_back(iElem);
    }
  }
  
  /* Post-process the neighbor element lists for each point. */
  vector<unsigned long>::iterator vecIt;
  for (unsigned long iPoint = 0; iPoint < n_points; iPoint++) {
    sort(neighborElemsOfPoint[iPoint].begin(), neighborElemsOfPoint[iPoint].end());
    vecIt = unique(neighborElemsOfPoint[iPoint].begin(), neighborElemsOfPoint[iPoint].end());
    neighborElemsOfPoint[iPoint].resize(vecIt - neighborElemsOfPoint[iPoint].begin());
  }
  
  /* Loop through all neighbor elements for each point and store
   the point IDs that we find in each element as neighboring points. */
  vector<vector<unsigned long> > neighborPointsOfPoint;
  neighborPointsOfPoint.resize(n_points);
  for (unsigned long iPoint = 0; iPoint < n_points; iPoint++) {
    for (unsigned long iElem = 0; iElem < neighborElemsOfPoint[iPoint].size(); iElem++) {
      
      /* loop over 3 points per triangle */
      for (unsigned long jPoint = 0; jPoint < N_POINTS_TRIANGLE; jPoint++) {
        
        /* Get the global ID of the current point. */
        const unsigned long GlobalIndex = triangles.at(neighborElemsOfPoint[iPoint][iElem]).at(jPoint);
        
        /* Add the current element ID to the neighbor list for this point. */
        if (GlobalIndex != iPoint)
          neighborPointsOfPoint[iPoint].push_back(GlobalIndex);
        
      }
    }
  }
  
  /* Post-process the neighbor point lists for each point. */
  for (unsigned long iPoint = 0; iPoint < n_points; iPoint++) {
    sort(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());
    vecIt = unique(neighborPointsOfPoint[iPoint].begin(), neighborPointsOfPoint[iPoint].end());
    neighborPointsOfPoint[iPoint].resize(vecIt - neighborPointsOfPoint[iPoint].begin());
  }
  
  /* Loop through our point neighbors and fill the vector of the unique
   point pairs making up each edge in the grid. We impose a requirement
   that the smaller global index is in the first position and the larger
   is in the second to make for a unique set, so there's no need to
   remove duplicates. */
  for (unsigned long iPoint = 0; iPoint < n_points; iPoint++) {
    for (unsigned long jPoint = 0; jPoint < neighborPointsOfPoint[iPoint].size(); jPoint++) {
      
      /* Store the neighbor index more clearly. */
      const unsigned long GlobalIndex = neighborPointsOfPoint[iPoint][jPoint];
      
      /* Store the edge so that the lower index of the pair is always first. */
      if (iPoint < GlobalIndex) {
        vector<unsigned long> edge(2);
        edge[0] = iPoint;
        edge[1] = GlobalIndex;
        edges.push_back(edge);
      }
      
    }
  }
  
  /* Loop over our edges data structure. For the first point in each
   pair, loop through the neighboring elements and store the two
   elements that contain the second point in the edge. */
  edge_to_triangle.resize(edges.size());
  for (unsigned long iEdge = 0; iEdge < edges.size(); iEdge++) {
    
    /* Store the two points of the edge more clearly. */
    const unsigned long iPoint = edges[iEdge][0];
    const unsigned long jPoint = edges[iEdge][1];
    
    /* Loop over all neighobring elements to iPoint. */
    for (unsigned long iElem = 0; iElem < neighborElemsOfPoint[iPoint].size(); iElem++) {
      
      /* loop over 3 points per triangle */
      for (unsigned long kPoint = 0; kPoint < N_POINTS_TRIANGLE; kPoint++) {
        
        /* Get the global ID of the current point. */
        const unsigned long GlobalIndex = triangles.at(neighborElemsOfPoint[iPoint][iElem]).at(kPoint);
        
        /* Add the current element ID to the neighbor list for this point. */
        if (GlobalIndex == jPoint)
          edge_to_triangle[iEdge].push_back(neighborElemsOfPoint[iPoint][iElem]);
        
      }
    }
  }
  
}

void CLookUpTable::ComputeInterpCoeffs() {
  
  /* build KD tree for enthalpy, progress variable space */
  vector< unsigned long > points(n_points);
  vector< su2double > weights(2, 0);

  vector< su2double > prog_enth_pairs(2 * n_points);

  vector< unsigned long > next_triangle;
  
  const vector<su2double> &prog = GetData("ProgVar");
  const vector<su2double> &enth = GetData("Enthalpy");

  vector< unsigned long > result_ids;
  vector< int > result_ranks;
  vector< su2double > best_dist;
  
  for (unsigned long i_point = 0; i_point < n_points; i_point++) {

    points.push_back(i_point);

    prog_enth_pairs.push_back(prog.at(i_point));
    prog_enth_pairs.push_back(enth.at(i_point));
  }

  /* calculate weights for each triangle (basically a distance function) and
   * build inverse interpolation matrices */
  for (unsigned long i_triangle = 0; i_triangle < n_triangles; i_triangle++) {

    next_triangle = triangles.at(i_triangle);

    /* the query point is the weighted average of the vertexes of the triangle */
    weights.at(0) = 0;
    weights.at(1) = 0;
    
    /* enthalpy */
    weights.at(0) += enth.at(next_triangle.at(0));
    weights.at(0) += enth.at(next_triangle.at(1));
    weights.at(0) += enth.at(next_triangle.at(2));
    weights.at(0) /= 3;

    /* progress variable */
    weights.at(1) += prog.at(next_triangle.at(0));
    weights.at(1) += prog.at(next_triangle.at(1));
    weights.at(1) += prog.at(next_triangle.at(2));
    weights.at(1) /= 3;

    interp_points.push_back(next_triangle);

    // Now use the nearest 16 points to construct an interpolation function
    // for each search pair option
    vector< vector< su2double > > prog_interp_mat_inv(3, vector< su2double >(3,0));
    GetInterpMatInv(prog, enth, next_triangle, prog_interp_mat_inv);
    interp_mat_inv_prog_enth.push_back(prog_interp_mat_inv);

  }
}

void CLookUpTable::GetInterpMatInv(const vector<su2double>       &vec_x,
                                   const vector<su2double>       &vec_y,
                                   vector<unsigned long>         &point_ids,
                                   vector< vector< su2double > > &interp_mat_inv) {

  vector<vector<su2double> > interp_mat(3, vector<su2double>(3, 0));

  /* setup LHM matrix for the interpolation */
  for (int i_point = 0; i_point < 3; i_point++) {

    su2double x = vec_x.at(point_ids.at(i_point));
    su2double y = vec_y.at(point_ids.at(i_point));

    interp_mat.at(i_point).at(0) = 1;
    interp_mat.at(i_point).at(1) = x;
    interp_mat.at(i_point).at(2) = y;
  }

  /* invert the Interpolation matrix using Gaussian elimination with pivoting */
  GaussianInverse(interp_mat, interp_mat_inv);

  /* transpose the inverse */
  su2double swap_helper;
  for (int i = 0; i < 2; i++) {
    for (int j = i + 1; j < 3; j++) {
      swap_helper                = interp_mat_inv.at(i).at(j);
      interp_mat_inv.at(i).at(j) = interp_mat_inv.at(j).at(i);
      interp_mat_inv.at(j).at(i) = swap_helper;
    }
  }
  
}

void CLookUpTable::GaussianInverse(vector< vector< su2double > > &mat,
                                   vector< vector< su2double > > &mat_inv) {

  /* temp provides memory to invert mat */
  vector< vector< su2double > > temp;

  /* number of dimensions of mat */
  int n_dim = (int)mat.size();

  temp.resize(n_dim, vector< su2double >(2 * n_dim, 0));

  /* copy original matrix into inverse */
  for (int i = 0; i < n_dim; i++) {
    for (int j = 0; j < n_dim; j++) {
      temp.at(i).at(j)         = mat.at(i).at(j);
      temp.at(i).at(n_dim + j) = 0;
    }
    temp.at(i).at(n_dim + i) = 1;
  }

  su2double max_val;
  int max_idx;
  /* pivot each column such that the largest number possible divides the other
   * rows The goal is to avoid zeros or small numbers in division. */
  for (int k = 0; k < n_dim - 1; k++) {
    max_idx = k;
    max_val = abs(temp.at(k).at(k));

    /* find largest value (pivot) in column */
    for (int j = k; j < n_dim; j++) {
      if (abs(temp.at(j).at(k)) > max_val) {
        max_idx = j;
        max_val = abs(temp.at(j).at(k));
      }
    }

    /* move row with the highest value up */
    for (int j = 0; j < (n_dim * 2); j++) {
      su2double d = temp.at(k).at(j);
      temp.at(k).at(j)       = temp.at(max_idx).at(j);
      temp.at(max_idx).at(j) = d;
    }

    /* subtract the moved row from all other rows */
    for (int i = k + 1; i < n_dim; i++) {
      su2double c = temp.at(i).at(k) / temp.at(k).at(k);
      for (int j = 0; j < (n_dim * 2); j++) {
        temp.at(i).at(j) = temp.at(i).at(j) - temp.at(k).at(j) * c;
      }
    }
  }

  /* perform back-substitution */
  for (int k = n_dim - 1; k > 0; k--) {
    if (temp.at(k).at(k) != 0) {
      for (int i = k - 1; i > -1; i--) {
        su2double c = temp.at(i).at(k) / temp.at(k).at(k);
        for (int j = 0; j < (n_dim * 2); j++) {
          temp.at(i).at(j) = temp.at(i).at(j) - temp.at(k).at(j) * c;
        }
      }
    }
  }

  /* normalize the inverse */
  for (int i = 0; i < n_dim; i++) {
    su2double c = temp.at(i).at(i);
    for (int j = 0; j < n_dim; j++) {
      temp.at(i).at(j + n_dim) = temp.at(i).at(j + n_dim) / c;
    }
  }

  /* copy inverse part into mat_inv */
  for (int i = 0; i < n_dim; i++) {
    for (int j = 0; j < n_dim; j++) {
      mat_inv.at(i).at(j) = temp.at(i).at(j + n_dim);
    }
  }
}

su2double CLookUpTable::LookUp_ProgEnth(string    name_var,
                                        su2double val_prog,
                                        su2double val_enth){
  
  /* find the triangle that holds the (prog, enth) point */
  unsigned long id_triangle = trap_map_prog_enth.GetTriangle(val_prog, val_enth);
  
  /* get interpolation coefficients for point on triangle */
  vector<su2double> interp_coeffs(3);
  GetInterpCoeffs(val_prog, val_enth, interp_mat_inv_prog_enth.at(id_triangle), interp_coeffs);
  
  return Interpolate(GetData(name_var), (triangles.at(id_triangle)), interp_coeffs);
  
}

  void CLookUpTable::GetInterpCoeffs(su2double                     val_x,
                                     su2double                     val_y,
                                     vector< vector< su2double > > &interp_mat_inv,
                                     vector<su2double>             &interp_coeffs) {

    vector<su2double> query_vector;
    query_vector.push_back(1);
    query_vector.push_back(val_x);
    query_vector.push_back(val_y);

    su2double d;
    for (int i = 0; i < 3; i++) {
      d = 0;
      for (int j = 0; j < 3; j++) {
        d = d + interp_mat_inv.at(i).at(j) * query_vector.at(j);
      }
      interp_coeffs.at(i) = d;
    }

}

su2double CLookUpTable::Interpolate(const vector<su2double> &val_samples,
                                    vector<unsigned long>   &val_triangle,
                                    vector<su2double>       &val_interp_coeffs) {

  su2double result         = 0;
  su2double z              = 0;

  for (int i_point = 0; i_point < N_POINTS_TRIANGLE; i_point++) {
    z = val_samples.at(val_triangle.at(i_point));
    result += val_interp_coeffs.at(i_point) * z;
  }

  return result;

}
