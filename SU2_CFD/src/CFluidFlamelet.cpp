/*!
 * \file CFluidFlamelet.cpp
 * \brief Definition of the flamelet fluid class members.
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

#include "../include/CFluidFlamelet.hpp"
#include "../include/CLookUpTable.hpp"

CFluidFlamelet::CFluidFlamelet(CConfig *config, su2double value_pressure_operating) : CFluidModel() {

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    look_up_table = new CLookUpTable(config->GetFileNameLUT());
 
    Pressure = value_pressure_operating;
}

CFluidFlamelet::~CFluidFlamelet() {
  delete look_up_table;
  //if (scalar_table != NULL) delete [] scalar_table;
}

void CFluidFlamelet::SetTDState_T(su2double val_temperature, su2double *val_scalars){

  su2double val_prog        = max(look_up_table->GetTableLimitsProg().first,
                              min(look_up_table->GetTableLimitsProg().second, val_scalars[I_PROG_VAR]));
  su2double val_enth        = max(look_up_table->GetTableLimitsEnth().first,
                              min(look_up_table->GetTableLimitsEnth().second, val_scalars[I_ENTHALPY]));
  

  Temperature        =         look_up_table->LookUp_ProgEnth("Temperature"             , val_prog, val_enth ) ;
  Density            =         look_up_table->LookUp_ProgEnth("Density"                 , val_prog, val_enth ) ;
  Cp                 =         look_up_table->LookUp_ProgEnth("Cp"                      , val_prog, val_enth ) ;
  Mu                 =         look_up_table->LookUp_ProgEnth("ViscosityDyn"            , val_prog, val_enth ) ;
  Kt                 =         look_up_table->LookUp_ProgEnth("Conductivity"            , val_prog, val_enth ) ;
  Diffusivity        =         look_up_table->LookUp_ProgEnth("Diffusivity"             , val_prog, val_enth ) ;
  source_energy      = max(0., look_up_table->LookUp_ProgEnth("HeatRelease"             , val_prog, val_enth)) ;

  source_prog        = max(0., look_up_table->LookUp_ProgEnth("ProdRateTot-PV"          , val_prog, val_enth)) ;

  su2double delta = 1e-6;

  /* Use backward differences for maximum pv values */
  if ( (val_prog + delta*val_prog) >=  look_up_table->GetTableLimitsProg().second)
    delta *= -1;

  dSourcePVdPV  = (max(0., look_up_table->LookUp_ProgEnth("ProdRateTot-PV", val_prog + delta*val_prog, val_enth)) - source_prog)/(delta*val_prog+EPS);

  // FIXME TDE: these are for the preconditioner. commented out for now
  // bc we need to check this on paper and they are extra cost.
  //dDensitydPV   = (max(0., look_up_table->LookUp_ProgEnth("Density", val_prog + delta*val_prog, val_enth)) - Density)/(delta*val_prog+EPS);
  //dDensitydEnth = (max(0., look_up_table->LookUp_ProgEnth("Density", val_prog, val_enth + delta*val_enth)) - Density)/(delta*val_enth+EPS);
  
  Cv = Cp - UNIVERSAL_GAS_CONSTANT / (look_up_table->LookUp_ProgEnth("MolarWeightMix", val_prog, val_enth) / 1000.);

}

void CFluidFlamelet::SetFluidFlameletTableOutput(CConfig *config, su2double val_temperature, su2double *val_scalars){

su2double val_prog        = max(look_up_table->GetTableLimitsProg().first,
                              min(look_up_table->GetTableLimitsProg().second, val_scalars[I_PROG_VAR]));
su2double val_enth        = max(look_up_table->GetTableLimitsEnth().first,
                              min(look_up_table->GetTableLimitsEnth().second, val_scalars[I_ENTHALPY]));

scalar_table = new su2double[config->GetnFlameletTableOutput()];
for (int ivar = 0; ivar<config->GetnFlameletTableOutput();ivar++)
scalar_table[ivar] = look_up_table->LookUp_ProgEnth(config->GetFlameletTableOutput_Field(ivar)                    , val_prog, val_enth ) ;

}





su2double CFluidFlamelet::LookUp_ProgEnth(string val_name_var, su2double val_prog, su2double val_enth){

  return look_up_table->LookUp_ProgEnth(val_name_var, val_prog, val_enth);

}
  
su2double CFluidFlamelet::GetEnthFromTemp(su2double val_prog, su2double val_temp){

  su2double  delta_temp_final = 0.01 ; /* convergence criterion for temperature in [K] */
  su2double  enth_iter        = 0. ;   /* in CH4/Air flames, 0 is usually a good initial value for the iteration */
  su2double  cp_iter          = 73 ;
  su2double  temp_iter        = 73 ;
  su2double  delta_enth       = 73 ;
  su2double  delta_temp_iter  = 73 ;
  int        counter_limit    = 50 ;
  val_prog         = max(look_up_table->GetTableLimitsProg().first,
                     min(look_up_table->GetTableLimitsProg().second, val_prog));
  

  int counter = 0;
  while ( (abs(delta_temp_iter) > delta_temp_final) && (counter++ < counter_limit) ){
    temp_iter       = LookUp_ProgEnth("Temperature", val_prog, enth_iter);
    delta_temp_iter = val_temp - temp_iter;
    cp_iter         = LookUp_ProgEnth("Cp", val_prog, enth_iter);
    delta_enth      = cp_iter * delta_temp_iter;
    enth_iter      += delta_enth;
  }

  if (counter >= counter_limit) {
    if (rank == MASTER_NODE) {
      cout << " WARNING: Could not iterate enthalpy from progress variable and temperature " << endl;
      cout << endl;
      cout << "   T_target   = " << val_temp << endl;
      cout << "   T_iter     = " << temp_iter << endl;
      cout << "   delta_T    = " << delta_temp_iter << endl;
      cout << "   delta_enth = " << delta_enth << endl;
      cout << "   counter    = " << counter << endl;
    }
  }
  return enth_iter;


}
