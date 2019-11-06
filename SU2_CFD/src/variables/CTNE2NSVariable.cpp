/*!
 * \file CTNE2NSVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon, S.R. Copeland, W. Maier
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

#include "../../include/variables/CTNE2NSVariable.hpp"
#include <math.h>

CTNE2NSVariable::CTNE2NSVariable(void) : CTNE2EulerVariable() { }

CTNE2NSVariable::CTNE2NSVariable(unsigned short val_ndim,
                                 unsigned short val_nvar,
                                 unsigned short val_nprimvar,
                                 unsigned short val_nprimvargrad,
                                 CConfig *config, CFluidModel *FluidModel) :
                                                    CTNE2EulerVariable(val_ndim,
                                                                       val_nvar,
                                                                       val_nprimvar,
                                                                       val_nprimvargrad,
                                                                       config, FluidModel) {

  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  DiffusionCoeff  = new su2double[nSpecies];
  Dij = new su2double*[nSpecies];
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Dij[iSpecies] = new su2double[nSpecies];
}

CTNE2NSVariable::CTNE2NSVariable(su2double val_pressure, su2double *val_massfrac,
                                 su2double *val_mach, su2double val_temperature,
                                 su2double val_temperature_ve,
                                 unsigned short val_ndim,
                                 unsigned short val_nvar,
                                 unsigned short val_nvarprim,
                                 unsigned short val_nvarprimgrad,
                                 CConfig *config, CFluidModel *FluidModel) :
                                                    CTNE2EulerVariable(val_pressure,
                                                                       val_massfrac,
                                                                       val_mach,
                                                                       val_temperature,
                                                                       val_temperature_ve,
                                                                       val_ndim,
                                                                       val_nvar,
                                                                       val_nvarprim,
                                                                       val_nvarprimgrad,
                                                                       config, FluidModel) {

  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  DiffusionCoeff  = new su2double[nSpecies];
  Dij = new su2double*[nSpecies];
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Dij[iSpecies] = new su2double[nSpecies];
}

CTNE2NSVariable::CTNE2NSVariable(su2double *val_solution, unsigned short val_ndim,
                                 unsigned short val_nvar,
                                 unsigned short val_nprimvar,
                                 unsigned short val_nprimvargrad,
                                 CConfig *config, CFluidModel *FluidModel) :
                                                    CTNE2EulerVariable(val_solution,
                                                                       val_ndim,
                                                                       val_nvar,
                                                                       val_nprimvar,
                                                                       val_nprimvargrad,
                                                                       config, FluidModel) {
  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  DiffusionCoeff  = new su2double[nSpecies];
  Dij = new su2double*[nSpecies];
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Dij[iSpecies] = new su2double[nSpecies];
}

CTNE2NSVariable::~CTNE2NSVariable(void) {
  // This cause invalid free (delete me) delete [] DiffusionCoeff;
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] Dij[iSpecies];
  delete [] Dij;
}

bool CTNE2NSVariable::SetVorticity(void) {

  su2double u_y = Gradient_Primitive[VEL_INDEX][1];
  su2double v_x = Gradient_Primitive[VEL_INDEX+1][0];
  su2double u_z = 0.0;
  su2double v_z = 0.0;
  su2double w_x = 0.0;
  su2double w_y = 0.0;

  if (nDim == 3) {
    u_z = Gradient_Primitive[VEL_INDEX][2];
    v_z = Gradient_Primitive[VEL_INDEX+1][2];
    w_x = Gradient_Primitive[VEL_INDEX+2][0];
    w_y = Gradient_Primitive[VEL_INDEX+2][1];
  }

  Vorticity[0] = w_y-v_z;
  Vorticity[1] = -(w_x-u_z);
  Vorticity[2] = v_x-u_y;

  return false;
}

bool CTNE2NSVariable::SetStrainMag(void) {

  su2double Div;
  unsigned short iDim;

  Div = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Div += Gradient_Primitive[VEL_INDEX+iDim][iDim];
  }

  StrainMag = 0.0;

  /*--- Add diagonal part ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    StrainMag += pow(Gradient_Primitive[VEL_INDEX+iDim][iDim] - 1.0/3.0*Div, 2.0);
  }
  if (nDim == 2) {
    StrainMag += pow(1.0/3.0*Div, 2.0);
  }

  /*--- Add off diagonals ---*/

  StrainMag += 2.0*pow(0.5*(Gradient_Primitive[VEL_INDEX][1] + Gradient_Primitive[VEL_INDEX+1][0]), 2.0);

  if (nDim == 3) {
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[VEL_INDEX][2]   + Gradient_Primitive[VEL_INDEX+2][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[VEL_INDEX+1][2] + Gradient_Primitive[VEL_INDEX+2][1]), 2.0);
  }

  StrainMag = sqrt(2.0*StrainMag);

  return false;

}

bool CTNE2NSVariable::SetPrimVar_Compressible(CConfig *config, CFluidModel *FluidModel) {

  bool nonPhys, bkup;
  unsigned short iVar;

  nonPhys = Cons2PrimVar(config, FluidModel, Solution, Primitive, dPdU, dTdU, dTvedU, eves, Cvves);
  if (nonPhys) {
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    bkup = Cons2PrimVar(config, FluidModel, Solution, Primitive, dPdU, dTdU, dTvedU, eves, Cvves);
  }

  SetVelocity2();

  FluidModel->Calc_TransportCoeff(config, Primitive);

  return nonPhys;
}

void CTNE2NSVariable::CalcEddyVisc_BL(su2double ymax){

  su2double A, S, l, I, Tau;
  su2double rho, mu_t, muti, muto;
  su2double y, yplus;


  /*--- Set parameters for model ---*/
  A   = 26;
  S   = GetStrainMag();
  //y   = GetWallDistance(); //Ask Eddie
  Tau = GetTauWall();

  /*--- Rename for convenience ---*/
  rho = Primitive[RHO_INDEX];

  /*--- Compute ---*/
  //yplus =

}
