/*!
 * \file CFileReaderLUT.hpp
 * \brief Declaration and inlines for the look-up table file reader class.
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

#include <string>
#include <vector>
#include <fstream>
#include "../../Common/include/mpi_structure.hpp"
/*#include "../../Common/include/datatypes/primitive_structure.hpp"*/

using namespace std;

class CFileReaderLUT {

 protected:

  int rank;

  string type_lut;
  string version_lut;
  string version_reader;
  unsigned long n_points;
  unsigned long n_triangles;
  unsigned long n_variables;

  /* !brief
   * Holds the variable names stored in the table file. Order is in sync with
   * tableFlamelet.
   */
  vector< string > names_var;
  
  /* !brief
   * Holds all data stored in the table. First index addresses the variable
   * while second index addresses the point.
   */
  vector< vector< su2double > > table_data;
  
  vector< vector< unsigned long > > triangles;

  string SkipToFlag(ifstream *file_stream, string flag);

  inline void SetTypeLUT(string value)    { type_lut    = value; }
  inline void SetVersionLUT(string value) { version_lut = value; }
  inline void SetNPoints(int value)       { n_points    = value; }
  inline void SetNTriangles(int value)    { n_triangles = value; }
  inline void SetNVariables(int value)    { n_variables = value; }
  
  inline void PushNameVar(string value) { names_var.push_back(value); }
  inline void PopNameVar() { names_var.pop_back(); }

  inline void AllocMemData() {
    table_data.resize(GetNVariables(), vector< su2double >(GetNPoints()));
  }
  
  inline void AllocMemTriangles() {
    triangles.resize(GetNTriangles(), vector< unsigned long >(3));
  }
  
 public:

  CFileReaderLUT();

  inline string GetTypeLUT() { return type_lut; }
  inline string GetVersionLUT() { return version_lut; }
  inline string GetVersionReader() { return version_reader; }
  inline int    GetNPoints() { return (int)n_points; }
  inline int    GetNTriangles() { return (int)n_triangles; }
  inline int    GetNVariables() { return (int)n_variables; }

  inline const vector< string >                  &GetNamesVar()  const { return names_var; }

  inline const vector< vector< su2double > >     &GetTableData() const { return table_data; }

  inline const vector< vector< unsigned long > > &GetTriangles() const { return triangles; };

  void ReadRawDRG(string file_name);

};
