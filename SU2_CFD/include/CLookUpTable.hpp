/*!
 * \file CLookUpTable.hpp
 * \brief Declaration and inlines for the look-up table class.
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

#pragma once

#include "../../Common/include/option_structure.hpp"
#include "CFileReaderLUT.hpp"
#include "CTrapezoidalMap.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

class CLookUpTable{

protected:

  int rank;  /*!< \brief MPI Rank. */

  string type_lut;
  string version_lut;
  string version_reader;
  unsigned long n_points;
  unsigned long n_triangles;
  unsigned long n_variables;

  su2double limits_table_enth[2];
  su2double limits_table_prog[2];

  /* !brief
   * Holds the variable names stored in the table file. 
   * Order is in sync with data
   */
  vector< string > names_var;

    /* !brief
   * Holds all data stored in the table. First index addresses the variable
   * while second index addresses the point.
   */
  vector< vector< su2double > > table_data;
    
  vector< vector< unsigned long > > triangles;
  vector< vector< unsigned long > > edges;
  vector< vector< unsigned long > > edge_to_triangle;

  CTrapezoidalMap trap_map_prog_enth;

  vector< vector< unsigned long > > interp_points;

  vector< vector< vector< su2double > > > interp_mat_inv_prog_enth;

  inline int GetIndexOfVar(string nameVar) {
    // FIXME dan: there should be an error handling when nameVar is not found in names_var
    return (int)(find(names_var.begin(), names_var.end(), nameVar) - names_var.begin());
  }
  
  inline const vector<su2double> &GetData(string name_var) {
    int ix_var = GetIndexOfVar(name_var);
    return table_data.at(ix_var);
  }
  
  inline const vector< vector< unsigned long > > &GetEdges() const {
    return edges;
  }
  
  inline const vector< vector< unsigned long > > &GetEdgeToTriangle() const {
    return edge_to_triangle;
  }
  
  void FindTableLimits();

  void IdentifyUniqueEdges();

  void LoadTableRaw(string file_name_lut);

  void ComputeInterpCoeffs();
  
  void GetInterpMatInv(const vector<su2double>       &vec_x,
                       const vector<su2double>       &vec_y,
                       vector<unsigned long>         &point_ids,
                       vector< vector< su2double > > &interp_mat_inv);

  void GetInterpCoeffs(su2double                     val_x,
                       su2double                     val_y,
                       vector< vector< su2double > > &interp_mat_inv,
                       vector<su2double>             &interp_coeffs);

  void GaussianInverse(vector<vector<su2double> > &mat,
                       vector<vector<su2double> > &mat_inv);

su2double Interpolate(const vector<su2double> &val_samples,
                      vector<unsigned long>   &val_triangle,
                      vector<su2double>       &val_interp_coeffs);

public:
CLookUpTable(string file_name_lut);

void PrintTableInfo();

su2double LookUp_ProgEnth(string name_var, su2double val_prog,
                          su2double val_enth);

inline pair<su2double, su2double> GetTableLimitsEnth() {
  return make_pair(limits_table_enth[0], limits_table_enth[1]);
}

inline pair<su2double, su2double> GetTableLimitsProg() {
  return make_pair(limits_table_prog[0], limits_table_prog[1]);
}

};
