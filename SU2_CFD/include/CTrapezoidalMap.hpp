/*!
 * \file CTrapezoidalMap.hpp
 * \brief Declaration and inlines for the trapezoidal map class.
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

#include <vector>

using namespace std;

class CTrapezoidalMap {
protected:

  /* The unique values of x which exist in the data */
  vector< su2double > unique_bands_x;

  vector< vector< su2double > > edge_limits_x;
  vector< vector< su2double > > edge_limits_y;

  vector< vector< unsigned long > > const *edge_to_triangle;

  /* The value that each edge which intersects the band takes within that
   * same band. Used to sort the edges */
  vector< vector< pair< su2double, unsigned long > > > y_edge_at_band_mid;

  
public:
  CTrapezoidalMap();

  CTrapezoidalMap(const vector< su2double >                &samples_x,
                  const vector< su2double >                &samples_y,
                  const vector< vector< unsigned long > >  &edges,
                  const vector< vector< unsigned long > >  &edge_to_triangle);

  ~CTrapezoidalMap(void);

  void Search_Band_For_Edge(su2double x, su2double y);

  unsigned long GetTriangle(su2double val_x, su2double val_y);

  pair<unsigned long, unsigned long> GetBand(su2double x);
  pair<unsigned long, unsigned long> GetEdges(
      pair<unsigned long, unsigned long> band, su2double val_x,
      su2double val_y);
};
