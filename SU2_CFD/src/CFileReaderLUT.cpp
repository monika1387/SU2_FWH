/*!
 * \file CFileReaderLUT.cpp
 * \brief Definition of the look-up table file reader class members.
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

#include "../include/CFileReaderLUT.hpp"
#include "../../Common/include/mpi_structure.hpp"
#include "../../Common/include/option_structure.hpp"
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;

CFileReaderLUT::CFileReaderLUT() {}

void CFileReaderLUT::ReadRawDRG(string file_name) {

  version_reader = "1.0.0";

  /*--- Store MPI rank. ---*/
  
  rank = SU2_MPI::GetRank();

  string line;
  string word;

  istringstream stream_names_var;

  ifstream file_stream;

  int ixColon;

  bool eoHeader       = false;
  bool eoData         = false;
  bool eoConnectivity = false;

  file_stream.open(file_name.c_str(), ifstream::in);

  if (!file_stream.is_open()) {
    SU2_MPI::Error(string("There is no look-up-table file file called ") + file_name,
                   CURRENT_FUNCTION);
  }

  /* Read header */
  line = SkipToFlag(&file_stream, "<header>");

  while (getline(file_stream, line) && !eoHeader) {
    /* number of points in LUT */
    if (line.compare("[version]") == 0) {
      getline(file_stream, line);
      SetVersionLUT(line);
    }

    /* number of points in LUT */
    if (line.compare("[number of points]") == 0) {
      getline(file_stream, line);
      SetNPoints(stoi(line));
    }

    /* number of triangles in LUT */
    if (line.compare("[number of triangles]") == 0) {
      getline(file_stream, line);
      SetNTriangles(stoi(line));
    }

    /* number of variables in LUT */
    if (line.compare("[number of variables]") == 0) {
      getline(file_stream, line);
      SetNVariables(stoi(line));
    }

    /* variable names */
    if (line.compare("[variable names]") == 0) {

      getline(file_stream, line);
      stream_names_var.str(line);
      while (stream_names_var) {
        stream_names_var >> word;
        ixColon = (int)word.find(":");

        PushNameVar(word.substr(ixColon + 1, word.size() - 1)); 
      }
      PopNameVar();  // removes last redundant element
    }

    // check if end of header is reached
    if (line.compare("</header>") == 0) eoHeader = true;
  }

  // check version_lut
  if (version_lut.compare(version_reader) != 0)
    SU2_MPI::Error(
        "Version conflict between Dragon reader and Dragon library file.",
        CURRENT_FUNCTION);

  // check header quantities
  if (n_points == 0 || n_triangles == 0 || n_variables == 0)
    SU2_MPI::Error(
        "Number of points, triangles, or variables in Dragon library header is "
        "zero.",
        CURRENT_FUNCTION);

  // check if number of variables is consistent
  if (n_variables != names_var.size())
    SU2_MPI::Error(
        "Number of read variables does not match number of "
        "variables specified in Dragon "
        "library header.",
        CURRENT_FUNCTION);

  // now that n_variables, n_points, and n_variables is available, allocate memory
  if (rank == MASTER_NODE)
    cout << "allocating memory for the data" << endl;
  AllocMemData();
  if (rank == MASTER_NODE)
    cout << "allocating memory for the triangles" << endl;
  AllocMemTriangles();

  /* flush any cout */
  if (rank == MASTER_NODE)
    cout << endl;

  // read data block
  line = SkipToFlag(&file_stream, "<data>");

  unsigned long pointCounter = 0;
    while (getline(file_stream, line) && !eoData) {
    // check if end of data is reached
    if (line.compare("</data>") == 0) eoData = true;

    if (!eoData) {

      // if (rank == MASTER_NODE)
      //   cout << "\r Reading point " << setw(11) << pointCounter + 1 << " / "
      //        << GetNPoints() << flush;

      // one line contains values for one point for all variables
      istringstream streamDataLine(line);

      // add next line to table array
      for (unsigned long iVar = 0; iVar < n_variables; iVar++) {
        streamDataLine >> word;
        passivedouble tmp = stod(word);
        table_data.at(iVar).at(pointCounter) = (su2double) tmp;
      }
    }
    pointCounter++;
  }

  if (n_points != pointCounter - 1)
    SU2_MPI::Error(
        "Number of read points does not match number of points "
        "specified in Dragon "
        "library header.",
        CURRENT_FUNCTION);

  // if (rank == MASTER_NODE)
  //   cout << " done. " << flush << endl;

  // read connectivity
  line = SkipToFlag(&file_stream, "<connectivity>");

  unsigned long triCounter = 0;
  while (getline(file_stream, line) && !eoConnectivity) {
    // check if end of data is reached
    if (line.compare("</connectivity>") == 0) eoConnectivity = true;

    if (!eoConnectivity) {

      // if (rank == MASTER_NODE)
      //   cout << "\r Reading triangle " << setw(8) << triCounter + 1 << " / "
      //        << GetNTriangles() << flush;

      // one line contains values for one triangle (3 points)
      istringstream streamTriLine(line);

      // add next line to triangles
      for (int iPoint = 0; iPoint < 3; iPoint++) {
        streamTriLine >> word;
        // Dragon table index starts with 1, convert to c++ indexing starting
        // with 0:
        triangles.at(triCounter).at(iPoint) = stol(word) - 1;
      }
    }
    triCounter++;
  }

  // if (rank == MASTER_NODE)
  //   cout << " done. " << flush << endl;

  file_stream.close();

  type_lut = "DRG";

}

string CFileReaderLUT::SkipToFlag(ifstream *file_stream, string flag) {
  string line;
  getline(*file_stream, line);

  while (line.compare(flag) != 0 && !(*file_stream).eof()) {
    getline(*file_stream, line);
  }

  if ((*file_stream).eof())
    SU2_MPI::Error("Flag not found in file", CURRENT_FUNCTION);

  return line;
}

