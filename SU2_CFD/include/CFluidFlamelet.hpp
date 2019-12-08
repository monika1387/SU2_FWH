/*!
 * \file CFluidFlamelet.hpp
 * \brief Declaration and inlines for the flamelet fluid class.
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

#include "../include/fluid_model.hpp"
#include "../include/CLookUpTable.hpp"


class CFluidFlamelet : public CFluidModel {

protected:

  int rank;

  CLookUpTable *look_up_table;
  
 public:
  CFluidFlamelet(CConfig *config, su2double value_pressure_operating);

  ~CFluidFlamelet();

  void SetTDState_T(su2double val_temperature, su2double *val_scalars);

  void SetFluidFlameletTableOutput(CConfig *config, su2double val_temperature, su2double *val_scalars);
  
  void SetTDState_prog_enth(su2double val_prog, su2double val_enth);

//  void SetTDState_prog_temp(su2double *val_prog, su2double TEMPERATURE);

  // void SetTDState_TY(su2double T, su2double *val_Scalar);

  // void SetTransTherm_TY(su2double T, su2double *val_Scalar);

  // void SetDiffusivity_TY(su2double T, su2double *val_Scalar);

  inline su2double GetMassDiffusivity() { return Diffusivity; }

  inline su2double GetThermalConductivity() { return Kt; }

  inline su2double GetLaminarViscosity() { return Mu; }

  su2double LookUp_ProgEnth(string val_name_var, su2double val_prog, su2double val_enth);

  su2double GetEnthFromTemp(su2double val_prog, su2double val_temp);

  inline pair<su2double, su2double> GetTableLimitsEnth() { return look_up_table->GetTableLimitsEnth(); }

  inline pair<su2double, su2double> GetTableLimitsProg() { return look_up_table->GetTableLimitsProg(); }

  inline su2double GetdDensitydPV() { return dDensitydPV; }

  inline su2double GetdSourcePVdPV() { return dSourcePVdPV; }

  inline su2double GetdDensitydEnth() { return dDensitydEnth; }

};
