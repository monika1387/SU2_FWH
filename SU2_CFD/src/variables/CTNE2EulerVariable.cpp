/*!
 * \file CTNE2EulerVariable.cpp
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

#include "../../include/variables/CTNE2EulerVariable.hpp"
#include <math.h>

CTNE2EulerVariable::CTNE2EulerVariable(void) : CVariable() {

  /*--- Array initialization ---*/
  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;
  Limiter            = NULL;
  nPrimVar           = 0;
  nPrimVarGrad       = 0;
  //Undivided_Laplacian = NULL;

  dPdU   = NULL;   dTdU   = NULL;
  dTvedU = NULL;   eves   = NULL;
  Cvves  = NULL;

  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;
  LAM_VISC_INDEX = nSpecies+nDim+8;
  EDDY_VISC_INDEX  = nSpecies+nDim+9;
  DIFF_COEFF_INDEX = nSpecies+nDim+10;
  K_INDEX          = 2*nSpecies+nDim+10;
  KVE_INDEX        = 2*nSpecies+nDim+11;
}

CTNE2EulerVariable::CTNE2EulerVariable(unsigned short val_ndim,
                                       unsigned short val_nvar,
                                       unsigned short val_nprimvar,
                                       unsigned short val_nprimvargrad,
                                       CConfig *config,
                                       CFluidModel *FluidModel) : CVariable(val_ndim,
                                                                            val_nvar,
                                                                            config) {

  nDim         = val_ndim;
  nVar         = val_nvar;

  nPrimVar     = val_nprimvar;
  nPrimVarGrad = val_nprimvargrad;

  nSpecies     = config->GetnSpecies();
  ionization   = config->GetIonization();

  /*--- Array initialization ---*/
  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;
  Limiter            = NULL;

  dPdU   = NULL;   dTdU   = NULL;
  dTvedU = NULL;   eves   = NULL;
  Cvves  = NULL;

  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;
  LAM_VISC_INDEX = nSpecies+nDim+8;
  EDDY_VISC_INDEX  = nSpecies+nDim+9;
  DIFF_COEFF_INDEX = nSpecies+nDim+10;
  K_INDEX          = 2*nSpecies+nDim+10;
  KVE_INDEX        = 2*nSpecies+nDim+11;


}

CTNE2EulerVariable::CTNE2EulerVariable(su2double val_pressure,
                                       su2double *val_massfrac,
                                       su2double *val_mach,
                                       su2double val_temperature,
                                       su2double val_temperature_ve,
                                       unsigned short val_ndim,
                                       unsigned short val_nvar,
                                       unsigned short val_nvarprim,
                                       unsigned short val_nvarprimgrad,
                                       CConfig *config,
                                       CFluidModel *FluidModel) : CVariable(val_ndim,
                                                                            val_nvar,
                                                                            config ) {

  unsigned short iEl, iMesh, iDim, iSpecies, iVar, nDim, nEl, nHeavy, nMGSmooth;
  unsigned short *nElStates;
  su2double *xi, *Ms, *thetav, **thetae, **g, *hf, *Tref;
  su2double rhoE, rhoEve, Ev, Ee, Ef, T, Tve, rho, rhoCvtr, rhos;
  su2double RuSI, Ru, sqvel, num, denom, conc, soundspeed;

  /*--- Define var counters ---*/
  nSpecies     = config->GetnSpecies();
  nDim         = val_ndim;
  nPrimVar     = val_nvarprim;
  nPrimVarGrad = val_nvarprimgrad;
  nMGSmooth    = 0;

  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;
  LAM_VISC_INDEX = nSpecies+nDim+8;
  EDDY_VISC_INDEX  = nSpecies+nDim+9;
  DIFF_COEFF_INDEX = nSpecies+nDim+10;
  K_INDEX          = 2*nSpecies+nDim+10;
  KVE_INDEX        = 2*nSpecies+nDim+11;


  /*--- Array initialization ---*/
  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;
  Limiter            = NULL;

  /*--- Allocate & initialize residual vectors ---*/
  Res_TruncError = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }

  /*--- If using multigrid, allocate residual-smoothing vectors ---*/
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }

  /*--- Allocate & initialize primitive variable & gradient arrays ---*/
  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) Limiter[iVar] = 0.0;

  Limiter_Primitive = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) Limiter_Primitive[iVar] = 0.0;

  //Limiter_Secondary = new su2double [nSecondaryVarGrad];
  //for (iVar = 0; iVar < nSecondaryVarGrad; iVar++) Limiter_Secondary[iVar] = 0.0;

  Solution_Max      = new su2double [nPrimVarGrad];
  Solution_Min      = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nVar; iVar++){
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }

  /*--- Initialize primitive variables ---*/
  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

  /*--- Initialize solution vectors ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar]     = 0.0;
    Solution_Old[iVar] = 0.0;
  }

  /*--- Allocate and initialize solution for dual time strategy ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    Solution_time_n[iVar] = 0.0;
//    Solution_time_n1[iVar] = 0.0;
//  }

  //Secondary = new su2double [nSecondaryVar];
  //for (iVar = 0; iVar < nSecondaryVar; iVar++) Secondary[iVar] = 0.0;

  /*--- Compressible flow, gradients primitive variables ---*/
  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }

  /*--- Allocate partial derivative vectors ---*/
  dPdU      = new su2double [nVar];
  dTdU      = new su2double [nVar];
  dTvedU    = new su2double [nVar];

  /*--- Allocate primitive vibrational energy arrays ---*/
  eves      = new su2double [nSpecies];
  Cvves     = new su2double [nSpecies];

  /*--- Determine the number of heavy species ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Load variables from the config class --*/
  xi        = FluidModel->GetRotationModes();      // Rotational modes of energy storage
  Ms        = FluidModel->GetMolar_Mass();         // Species molar mass
  thetav    = FluidModel->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  thetae    = FluidModel->GetCharElTemp();         // Characteristic electron temperature [K]
  g         = FluidModel->GetElDegeneracy();       // Degeneracy of electron states
  nElStates = FluidModel->GetnElStates();          // Number of electron states
  Tref      = FluidModel->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = FluidModel->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]

  /*--- Rename & initialize for convenience ---*/
  RuSI      = UNIVERSAL_GAS_CONSTANT;          // Universal gas constant [J/(mol*K)]
  Ru        = 1000.0*RuSI;                     // Universal gas constant [J/(kmol*K)]
  Tve       = val_temperature_ve;              // Vibrational temperature [K]
  T         = val_temperature;                 // Translational-rotational temperature [K]
  sqvel     = 0.0;                             // Velocity^2 [m2/s2]
  rhoE      = 0.0;                             // Mixture total energy per mass [J/kg]
  rhoEve    = 0.0;                             // Mixture vib-el energy per mass [J/kg]
  denom     = 0.0;
  conc      = 0.0;
  rhoCvtr   = 0.0;

  /*--- Calculate mixture density from supplied primitive quantities ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    denom += val_massfrac[iSpecies] * (Ru/Ms[iSpecies]) * T;

  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    denom += val_massfrac[nSpecies-1] * (Ru/Ms[nSpecies-1]) * Tve;

  rho = val_pressure / denom;

  /*--- Calculate sound speed and extract velocities ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    conc += val_massfrac[iSpecies]*rho/Ms[iSpecies];
    rhoCvtr += rho*val_massfrac[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
  }

  soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * val_pressure/rho);

  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += val_mach[iDim]*soundspeed * val_mach[iDim]*soundspeed;

  /*--- Calculate energy (RRHO) from supplied primitive quanitites ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    // Species density
    rhos = val_massfrac[iSpecies]*rho;

    // Species formation energy
    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

    // Species vibrational energy
    if (thetav[iSpecies] != 0.0)
      Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;

    // Species electronic energy
    num = 0.0;
    denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);

    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
    }

    Ee = Ru/Ms[iSpecies] * (num/denom);

    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                    + Ev + Ee + Ef + 0.5*sqvel);

    // Mixture vibrational-electronic energy
    rhoEve += rhos * (Ev + Ee);
  }

  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
    // Species formation energy
    Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];

    // Electron t-r mode contributes to mixture vib-el energy
    rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
  }

  /*--- Initialize Solution & Solution_Old vectors ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Solution[iSpecies]     = rho*val_massfrac[iSpecies];
    Solution_Old[iSpecies] = rho*val_massfrac[iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    Solution[nSpecies+iDim]     = rho*val_mach[iDim]*soundspeed;
    Solution_Old[nSpecies+iDim] = rho*val_mach[iDim]*soundspeed;
  }
  Solution[nSpecies+nDim]       = rhoE;
  Solution_Old[nSpecies+nDim]   = rhoE;
  Solution[nSpecies+nDim+1]     = rhoEve;
  Solution_Old[nSpecies+nDim+1] = rhoEve;

  /*--- Assign primitive variables ---*/
  Primitive[T_INDEX]   = val_temperature;
  Primitive[TVE_INDEX] = val_temperature_ve;
  Primitive[P_INDEX]   = val_pressure;

}

CTNE2EulerVariable::CTNE2EulerVariable(su2double *val_solution,
                                       unsigned short val_ndim,
                                       unsigned short val_nvar,
                                       unsigned short val_nvarprim,
                                       unsigned short val_nvarprimgrad,
                                       CConfig *config,
                                       CFluidModel *FluidModel) : CVariable(val_ndim,
                                                                            val_nvar,
                                                                            config) {

  unsigned short iVar, iDim, iMesh, nMGSmooth;

  nSpecies     = config->GetnSpecies();
  nDim         = val_ndim;
  nPrimVar     = val_nvarprim;
  nPrimVarGrad = val_nvarprimgrad;
  nMGSmooth    = 0;

  /*--- Define structure of the primtive variable vector ---*/
  // Primitive: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T
  // GradPrim:  [rho1, ..., rhoNs, T, Tve, u, v, w, P]^T
  RHOS_INDEX    = 0;
  T_INDEX       = nSpecies;
  TVE_INDEX     = nSpecies+1;
  VEL_INDEX     = nSpecies+2;
  P_INDEX       = nSpecies+nDim+2;
  RHO_INDEX     = nSpecies+nDim+3;
  H_INDEX       = nSpecies+nDim+4;
  A_INDEX       = nSpecies+nDim+5;
  RHOCVTR_INDEX = nSpecies+nDim+6;
  RHOCVVE_INDEX = nSpecies+nDim+7;
  LAM_VISC_INDEX = nSpecies+nDim+8;
  EDDY_VISC_INDEX  = nSpecies+nDim+9;
  DIFF_COEFF_INDEX = nSpecies+nDim+10;
  K_INDEX          = 2*nSpecies+nDim+10;
  KVE_INDEX        = 2*nSpecies+nDim+11;
  /*--- Array initialization ---*/
  Primitive          = NULL;
  Gradient_Primitive = NULL;
  Limiter_Primitive  = NULL;
  Limiter            = NULL;

  /*--- Allocate & initialize residual vectors ---*/
  Res_TruncError = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Res_TruncError[iVar] = 0.0;
  }

  /*--- If using multigrid, allocate residual-smoothing vectors ---*/
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++)
    nMGSmooth += config->GetMG_CorrecSmooth(iMesh);

  if (nMGSmooth > 0) {
    Residual_Sum = new su2double [nVar];
    Residual_Old = new su2double [nVar];
  }

  /*--- If using limiters, allocate the arrays ---*/
  Limiter_Primitive = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Limiter_Primitive[iVar] = 0.0;

  Limiter = new su2double [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Limiter[iVar] = 0.0;

  Solution_Max = new su2double [nPrimVarGrad];
  Solution_Min = new su2double [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Solution_Max[iVar] = 0.0;
    Solution_Min[iVar] = 0.0;
  }

  /*--- Allocate & initialize primitive variable & gradient arrays ---*/
  Primitive = new su2double [nPrimVar];
  for (iVar = 0; iVar < nPrimVar; iVar++) Primitive[iVar] = 0.0;

  Gradient_Primitive = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Gradient_Primitive[iVar] = new su2double [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
  }

  /*--- Allocate partial derivative vectors ---*/
  dPdU   = new su2double [nVar];
  dTdU   = new su2double [nVar];
  dTvedU = new su2double [nVar];

  /*--- Allocate vibrational-electronic arrays ---*/
  eves  = new su2double[nSpecies];
  Cvves = new su2double[nSpecies];

  /*--- Determine the number of heavy species ---*/
  ionization = config->GetIonization();

  /*--- Initialize Solution & Solution_Old vectors ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Solution[iVar]     = val_solution[iVar];
    Solution_Old[iVar] = val_solution[iVar];
  }

  /*--- Initialize Tve to the free stream for Newton-Raphson method ---*/
  Primitive[TVE_INDEX] = config->GetTemperature_FreeStream();
  Primitive[T_INDEX]   = config->GetTemperature_FreeStream();
  Primitive[P_INDEX]   = config->GetPressure_FreeStream();
}

CTNE2EulerVariable::~CTNE2EulerVariable(void) {

  unsigned short iVar;

  if (Primitive          != NULL) delete [] Primitive;
  if (Limiter_Primitive  != NULL) delete [] Limiter_Primitive;

  if (Gradient_Primitive != NULL) {
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      if (Gradient_Primitive[iVar] != NULL) delete [] Gradient_Primitive[iVar];
    delete [] Gradient_Primitive;
  }

  if (dPdU   != NULL) delete [] dPdU;
  if (dTdU   != NULL) delete [] dTdU;
  if (dTvedU != NULL) delete [] dTvedU;
  if (eves   != NULL) delete [] eves;
  if (Cvves  != NULL) delete [] Cvves;
}

void CTNE2EulerVariable::SetGradient_PrimitiveZero(unsigned short val_primvar) {

  unsigned short iVar, iDim;

  for (iVar = 0; iVar < val_primvar; iVar++)
    for (iDim = 0; iDim < nDim; iDim++)
      Gradient_Primitive[iVar][iDim] = 0.0;
}

su2double CTNE2EulerVariable::GetProjVel(su2double *val_vector) {

  su2double ProjVel, density;
  unsigned short iDim, iSpecies;

  ProjVel = 0.0;
  density = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    density += Solution[iSpecies];
  for (iDim = 0; iDim < nDim; iDim++)
    ProjVel += Solution[nSpecies+iDim]*val_vector[iDim]/density;

  return ProjVel;
}

void CTNE2EulerVariable::SetVelocity2(void) {

  unsigned short iDim;

  Velocity2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Primitive[VEL_INDEX+iDim] = Solution[nSpecies+iDim] / Primitive[RHO_INDEX];
    Velocity2 +=  Solution[nSpecies+iDim]*Solution[nSpecies+iDim]
        / (Primitive[RHO_INDEX]*Primitive[RHO_INDEX]);
  }
}

bool CTNE2EulerVariable::SetPrimVar_Compressible(CConfig *config, CFluidModel *FluidModel) {

  bool nonPhys, bkup;
  unsigned short iVar;

  /*--- Convert conserved to primitive variables ---*/
  nonPhys = Cons2PrimVar(config, FluidModel, Solution, Primitive, dPdU, dTdU, dTvedU, eves, Cvves);
  if (nonPhys) {
    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];
    bkup = Cons2PrimVar(config, FluidModel, Solution, Primitive, dPdU, dTdU, dTvedU, eves, Cvves);
  }

  SetVelocity2();

  return nonPhys;
}

bool CTNE2EulerVariable::Cons2PrimVar(CConfig *config, CFluidModel *FluidModel,
                                      su2double *U, su2double *V,
                                      su2double *val_dPdU, su2double *val_dTdU,
                                      su2double *val_dTvedU, su2double *val_eves,
                                      su2double *val_Cvves) {

  bool ionization, errT, errTve, NRconvg, Bconvg, nonPhys;
  unsigned short iDim, iEl, iSpecies, nHeavy, nEl, iIter, maxBIter, maxNIter;
  su2double rho, rhoE, rhoEve, rhoE_f, rhoE_ref, rhoEve_min, rhoEve_max, rhoEve_t;
  su2double RuSI, Ru, sqvel, rhoCvtr, rhoCvve;
  su2double Tve, Tve2, Tve_o;
  su2double f, df, NRtol, Btol, scale;
  su2double Tmin, Tmax, Tvemin, Tvemax;
  su2double radical2;
  su2double *xi, *Ms, *hf, *Tref;
  unsigned short iVar;

  /*--- Conserved & primitive vector layout ---*/
  // U:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  // V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T

  /*--- Set booleans ---*/
  errT    = false;
  errTve  = false;
  nonPhys = false;

  /*--- Set temperature clipping values ---*/
  Tmin   = 50.0; Tmax   = 8E4;
  Tvemin = 50.0; Tvemax = 8E4;

  /*--- Set temperature algorithm paramters ---*/
  NRtol    = 1.0E-6;    // Tolerance for the Newton-Raphson method
  Btol     = 1.0E-4;    // Tolerance for the Bisection method
  maxNIter = 18;        // Maximum Newton-Raphson iterations
  maxBIter = 32;        // Maximum Bisection method iterations
  scale    = 0.5;       // Scaling factor for Newton-Raphson step

  /*--- Initialize primitive vector ---*/
  for (iVar = 0; iVar < nPrimVar; iVar++)
    V[iVar] = 0.0;

  /*--- Read parameters from config ---*/
  xi         = FluidModel->GetRotationModes();      // Rotational modes of energy storage
  Ms         = FluidModel->GetMolar_Mass();         // Species molar mass
  Tref       = FluidModel->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf         = FluidModel->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]
  ionization = FluidModel->GetIonization();         // Molecule Ionization

  /*--- Rename variables for convenience ---*/
  RuSI   = UNIVERSAL_GAS_CONSTANT;    // Universal gas constant [J/(mol*K)]
  Ru     = 1000.0*RuSI;               // Universal gas constant [J/(kmol*K)]
  rhoE   = U[nSpecies+nDim];          // Density * energy [J/m3]
  rhoEve = U[nSpecies+nDim+1];        // Density * energy_ve [J/m3]

  /*--- Determine the number of heavy species ---*/
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Assign species & mixture density ---*/
  // Note: if any species densities are < 0, these values are re-assigned
  //       in the primitive AND conserved vectors to ensure positive density
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    if (U[iSpecies] < 0.0) {
      V[RHOS_INDEX+iSpecies] = 1E-20;
      U[iSpecies]            = 1E-20;
      V[RHO_INDEX]           = 1E-20;
      nonPhys                = true;
    } else
      V[RHOS_INDEX+iSpecies] = U[iSpecies];
      V[RHO_INDEX]          += U[iSpecies];
  }

  /*--- Assign mixture velocity ---*/
  sqvel = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    V[VEL_INDEX+iDim] = U[nSpecies+iDim]/V[RHO_INDEX];
    sqvel            += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];
  }

  /*--- Translational-Rotational Temperature ---*/

  // Rename for convenience
  rho = V[RHO_INDEX];

  // Determine properties of the mixture at the current state
  rhoE_f   = 0.0;
  rhoE_ref = 0.0;
  rhoCvtr  = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhoCvtr  += U[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies];
    rhoE_ref += U[iSpecies] * (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * Tref[iSpecies];
    rhoE_f   += U[iSpecies] * (hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies]);
  }

  // Calculate T-R temperature
  V[T_INDEX] = (rhoE - rhoEve - rhoE_f + rhoE_ref - 0.5*rho*sqvel) / rhoCvtr;
  V[RHOCVTR_INDEX] = rhoCvtr;

  // Determine if the temperature lies within the acceptable range
  if (V[T_INDEX] < Tmin) {
    V[T_INDEX] = Tmin;
    nonPhys = true;
    errT    = true;
  } else if (V[T_INDEX] > Tmax){
    V[T_INDEX] = Tmax;
    nonPhys = true;
    errT    = true;
  }

  /*--- Vibrational-Electronic Temperature ---*/

  // Check for non-physical solutions
  rhoEve_min = 0.0;
  rhoEve_max = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhoEve_min += U[iSpecies]*FluidModel->CalcEve(Tvemin, iSpecies);
    rhoEve_max += U[iSpecies]*FluidModel->CalcEve(Tvemax, iSpecies);
  }
  if (rhoEve < rhoEve_min) {
    errTve       = true;
    nonPhys      = true;
    V[TVE_INDEX] = Tvemin;
  } else if (rhoEve > rhoEve_max) {
    errTve       = true;
    nonPhys      = true;
    V[TVE_INDEX] = Tvemax;
  } else {

    /*--- Execute a Newton-Raphson root-finding method to find Tve ---*/
    // Initialize to the translational-rotational temperature
    Tve   = V[T_INDEX];

    // Execute the root-finding method
    NRconvg = false;
    //    for (iIter = 0; iIter < maxNIter; iIter++) {
    //      rhoEve_t = 0.0;
    //      rhoCvve  = 0.0;
    //      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    //        val_eves[iSpecies]  = CalcEve( Tve, iSpecies);
    //        val_Cvves[iSpecies] = CalcCvve(Tve, iSpecies);
    //        rhoEve_t += U[iSpecies]*val_eves[iSpecies];
    //        rhoCvve  += U[iSpecies]*val_Cvves[iSpecies];
    //      }
    //
    //      // Find the root
    //      f  = U[nSpecies+nDim+1] - rhoEve_t;
    //      df = -rhoCvve;
    //      Tve2 = Tve - (f/df)*scale;
    //
    //      // Check for nonphysical steps
    //      if ((Tve2 < Tvemin) || (Tve2 > Tvemax))
    //        break;
    ////      if (Tve2 < Tvemin)
    ////        Tve2 = Tvemin;
    ////      else if (Tve2 > Tvemax)
    ////        Tve2 = Tvemax;
    //
    //      // Check for convergence
    //      if (fabs(Tve2 - Tve) < NRtol) {
    //        NRconvg = true;
    //        break;
    //      } else {
    //        Tve = Tve2;
    //      }
    //    }

    // If the Newton-Raphson method has converged, assign the value of Tve.
    // Otherwise, execute a bisection root-finding method
    if (NRconvg){
      V[TVE_INDEX] = Tve;
     } else {

      // Assign the bounds
      Tve_o = Tvemin;
      Tve2  = Tvemax;

      // Execute the root-finding method
      Bconvg = false;
      for (iIter = 0; iIter < maxBIter; iIter++) {
        Tve      = (Tve_o+Tve2)/2.0;
        rhoEve_t = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          val_eves[iSpecies] = FluidModel->CalcEve(Tve, iSpecies);
          rhoEve_t          += U[iSpecies] * val_eves[iSpecies];
        }

        if (fabs(rhoEve_t - U[nSpecies+nDim+1]) < Btol) {
          V[TVE_INDEX] = Tve;
          for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
            val_Cvves[iSpecies] = FluidModel -> Calc_CvVibElSpecies(Tve, iSpecies);
          Bconvg = true;
          break;
        } else {
          if (rhoEve_t > rhoEve) Tve2 = Tve;
          else                  Tve_o = Tve;
        }
      }

      // If absolutely no convergence, then assign to the TR temperature
      if (!Bconvg) {
        V[TVE_INDEX] = V[T_INDEX];
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          val_eves[iSpecies]  = FluidModel->CalcEve(V[TVE_INDEX], iSpecies);
          val_Cvves[iSpecies] = FluidModel->Calc_CvVibElSpecies(V[TVE_INDEX], iSpecies);
        }
      }
    }
  }

  /*--- Set mixture rhoCvve ---*/
  rhoCvve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhoCvve += U[iSpecies]*val_Cvves[iSpecies];
  V[RHOCVVE_INDEX] = rhoCvve;

  /*--- If there are clipped temperatures, correct the energy terms ---*/
  //  if (errT) {
  //    U[nSpecies+nDim]   = rhoCvtr*V[T_INDEX] + rhoCvve*V[TVE_INDEX] + rhoE_f
  //                       - rhoE_ref + 0.5*rho*sqvel;
  //  }
  //  if (errTve) {
  //    U[nSpecies+nDim]   = rhoCvtr*V[T_INDEX] + rhoCvve*V[TVE_INDEX] + rhoE_f
  //                       - rhoE_ref + 0.5*rho*sqvel;
  //    U[nSpecies+nDim+1] = rhoCvve*V[TVE_INDEX];
  //  }

  /*--- Pressure ---*/
  V[P_INDEX] = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    V[P_INDEX] += U[iSpecies] * Ru/Ms[iSpecies] * V[T_INDEX];
  for (iEl = 0; iEl < nEl; iEl++)
    V[P_INDEX] += U[nSpecies-1] * Ru/Ms[nSpecies-1] * V[TVE_INDEX];

  if (V[P_INDEX] < 0.0) {
    V[P_INDEX] = 1E-20;
    nonPhys = true;
  }

  /*--- Partial derivatives of pressure and temperature ---*/
  FluidModel -> CalcdPdU(  V, val_eves, val_dPdU  );
  FluidModel -> CalcdTdU(  V, val_dTdU  );
  FluidModel -> CalcdTvedU(V, val_eves, val_dTvedU);

  /*--- Sound speed ---*/
  radical2 = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    radical2 += V[RHOS_INDEX+iSpecies]/V[RHO_INDEX] * val_dPdU[iSpecies];
  for (iDim = 0; iDim < nDim; iDim++)
    radical2 += V[VEL_INDEX+iDim]*val_dPdU[nSpecies+iDim];
  radical2 += (U[nSpecies+nDim]+V[P_INDEX])/V[RHO_INDEX] * val_dPdU[nSpecies+nDim];
  radical2 += U[nSpecies+nDim+1]/V[RHO_INDEX] * val_dPdU[nSpecies+nDim+1];
  V[A_INDEX] = sqrt(radical2);

  if (radical2 < 0.0) {
    nonPhys = true;
    V[A_INDEX] = EPS;
  }

  /*--- Enthalpy ---*/
  V[H_INDEX] = (U[nSpecies+nDim] + V[P_INDEX])/V[RHO_INDEX];

  return nonPhys;
}

void CTNE2EulerVariable::Prim2ConsVar(CFluidModel *FluidModel, su2double *V, su2double *U) {
  unsigned short iDim, iEl, iSpecies, nEl, nHeavy;
  unsigned short *nElStates;
  su2double Ru, RuSI, Tve, T, sqvel, rhoE, rhoEve, Ef, Ev, Ee, rhos, rhoCvtr, num, denom;
  su2double *thetav, *Ms, *xi, *hf, *Tref;
  su2double **thetae, **g;

  /*--- Determine the number of heavy species ---*/
  ionization = FluidModel->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Load variables from the config class --*/
  xi        = FluidModel->GetRotationModes();      // Rotational modes of energy storage
  Ms        = FluidModel->GetMolar_Mass();         // Species molar mass
  thetav    = FluidModel->GetCharVibTemp();        // Species characteristic vib. temperature [K]
  thetae    = FluidModel->GetCharElTemp();         // Characteristic electron temperature [K]
  g         = FluidModel->GetElDegeneracy();       // Degeneracy of electron states
  nElStates = FluidModel->GetnElStates();          // Number of electron states
  Tref      = FluidModel->GetRefTemperature();     // Thermodynamic reference temperature [K]
  hf        = FluidModel->GetEnthalpy_Formation(); // Formation enthalpy [J/kg]

  /*--- Rename & initialize for convenience ---*/
  RuSI    = UNIVERSAL_GAS_CONSTANT;         // Universal gas constant [J/(mol*K)] (SI units)
  Ru      = 1000.0*RuSI;                    // Universal gas constant [J/(kmol*K)]
  Tve     = V[TVE_INDEX];                   // Vibrational temperature [K]
  T       = V[T_INDEX];                     // Translational-rotational temperature [K]
  sqvel   = 0.0;                            // Velocity^2 [m2/s2]
  rhoE    = 0.0;                            // Mixture total energy per mass [J/kg]
  rhoEve  = 0.0;                            // Mixture vib-el energy per mass [J/kg]
  denom   = 0.0;
  rhoCvtr = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];

  /*--- Set species density ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    U[iSpecies] = V[RHOS_INDEX+iSpecies];

  /*--- Set momentum ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    U[nSpecies+iDim] = V[RHO_INDEX]*V[VEL_INDEX+iDim];

  /*--- Set the total energy ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhos = U[iSpecies];

    // Species formation energy
    Ef = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];

    // Species vibrational energy
    if (thetav[iSpecies] != 0.0)
      Ev = Ru/Ms[iSpecies] * thetav[iSpecies] / (exp(thetav[iSpecies]/Tve)-1.0);
    else
      Ev = 0.0;

    // Species electronic energy
    num = 0.0;
    denom = g[iSpecies][0] * exp(thetae[iSpecies][0]/Tve);
    for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
    }
    Ee = Ru/Ms[iSpecies] * (num/denom);

    // Mixture total energy
    rhoE += rhos * ((3.0/2.0+xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * (T-Tref[iSpecies])
                    + Ev + Ee + Ef + 0.5*sqvel);

    // Mixture vibrational-electronic energy
    rhoEve += rhos * (Ev + Ee);
  }
  for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
    // Species formation energy
    Ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1] * Tref[nSpecies-1];

    // Electron t-r mode contributes to mixture vib-el energy
    rhoEve += (3.0/2.0) * Ru/Ms[nSpecies-1] * (Tve - Tref[nSpecies-1]);
  }

  /*--- Set energies ---*/
  U[nSpecies+nDim]   = rhoE;
  U[nSpecies+nDim+1] = rhoEve;

  return;
}

bool CTNE2EulerVariable::GradCons2GradPrimVar(CConfig *config, CFluidModel *FluidModel, su2double *U,
                                              su2double *V, su2double **GradU,
                                              su2double **GradV) {

  unsigned short iSpecies, iEl, iDim, jDim, iVar, nHeavy, *nElStates;
  su2double rho, rhoCvtr, rhoCvve, T, Tve, eve, ef, eref, RuSI, Ru;
  su2double Cvvs, Cves, dCvvs, dCves;
  su2double thoTve, exptv;
  su2double An1, Bd1, Bn1, Bn2, Bn3, Bn4, A, B;
  su2double *rhou, *Grhou2;
  su2double *xi, *Ms, *Tref, *hf, *thetav;
  su2double **thetae, **g;

  /*--- Conserved & primitive vector layout ---*/
  // U:  [rho1, ..., rhoNs, rhou, rhov, rhow, rhoe, rhoeve]^T
  // V: [rho1, ..., rhoNs, T, Tve, u, v, w, P, rho, h, a, rhoCvtr, rhoCvve]^T

  /*--- Allocate arrays ---*/
  Grhou2 = new su2double[nDim];
  rhou   = new su2double[nDim];

  /*--- Determine number of heavy-particle species ---*/
  if (config->GetIonization()) { nHeavy = nSpecies-1; }
  else                         { nHeavy = nSpecies;   }

  /*--- Rename for convenience ---*/
  rho       = V[RHO_INDEX];
  rhoCvtr   = V[RHOCVTR_INDEX];
  rhoCvve   = V[RHOCVVE_INDEX];
  T         = V[T_INDEX];
  Tve       = V[TVE_INDEX];
  xi        = FluidModel->GetRotationModes();
  Ms        = FluidModel->GetMolar_Mass();
  Tref      = FluidModel->GetRefTemperature();
  hf        = FluidModel->GetEnthalpy_Formation();
  RuSI      = UNIVERSAL_GAS_CONSTANT;
  Ru        = 1000.0*RuSI;
  thetav    = FluidModel->GetCharVibTemp();
  g         = FluidModel->GetElDegeneracy();
  thetae    = FluidModel->GetCharElTemp();
  nElStates = FluidModel->GetnElStates();

  for (iDim = 0; iDim < nDim; iDim++)
    for (iVar = 0; iVar < nPrimVarGrad; iVar++)
      GradV[iVar][iDim] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++) {

    Grhou2[iDim] = 0.0;

    /*--- Species density ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      GradV[RHOS_INDEX+iSpecies][iDim] = GradU[iSpecies][iDim];
      GradV[RHO_INDEX][iDim]          += GradU[iSpecies][iDim];
    }

    /*--- Velocity ---*/
    for (jDim = 0; jDim < nDim; jDim++) {
      GradV[VEL_INDEX+jDim][iDim] = (rho*GradU[nSpecies+jDim][iDim] -
          rhou[jDim]*GradV[RHO_INDEX][iDim]) / (rho*rho);
    }

    /*--- Specific Heat (T-R) ---*/
    GradV[RHOCVTR_INDEX][iDim] = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      GradV[RHOCVTR_INDEX][iDim] += GradU[iSpecies][iDim]*(3.0+xi[iSpecies])/2.0 *Ru/Ms[iSpecies];

    /*--- Temperature ---*/
    // Calculate the gradient of rho*u^2
    for (jDim = 0; jDim < nDim; jDim++) {
      Grhou2[iDim] += 2.0/rho*(rhou[jDim]*GradU[nSpecies+jDim][iDim]) -
          GradV[RHO_INDEX][iDim]/(rho*rho)*(GradU[nSpecies+jDim][iDim]*
          GradU[nSpecies+jDim][iDim]);
    }
    // Calculate baseline GradT
    GradV[T_INDEX][iDim] = 1.0/rhoCvtr*(GradU[nSpecies+nDim][iDim] -
        GradU[nSpecies+nDim+1][iDim] -
        0.5*Grhou2[iDim] -
        T*GradV[RHOCVTR_INDEX][iDim]);
    // Subtract formation/reference energies
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      eref = (3.0/2.0 + xi[iSpecies]/2.0) * Ru/Ms[iSpecies] * Tref[iSpecies];
      ef   = (hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies]);
      GradV[T_INDEX][iDim] -= 1.0/rhoCvtr*(GradU[iSpecies][iDim]*(ef+eref));
    }

    /*--- Vibrational-electronic temperature ---*/
    GradV[TVE_INDEX][iDim] = GradU[nSpecies+nDim+1][iDim]/rhoCvve;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      eve = CalcEve(config, V[TVE_INDEX], iSpecies);
      GradV[TVE_INDEX][iDim] -= GradU[iSpecies][iDim]*eve;
    }

    /*--- Pressure ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      GradV[P_INDEX][iDim] += GradU[iSpecies][iDim]*Ru/Ms[iSpecies]*T +
          U[iSpecies]*Ru/Ms[iSpecies]*GradV[T_INDEX][iDim];
    }

    /*--- Enthalpy ---*/
    GradV[H_INDEX][iDim] = rho*(GradU[nSpecies+nDim][iDim] + GradV[P_INDEX][iDim]) -
        (U[nSpecies+nDim]+V[P_INDEX])*GradV[RHO_INDEX][iDim] / (rho*rho);

    /*--- Specific Heat (V-E) ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      // Vibrational energy specific heat
      if (thetav[iSpecies] != 0) {
        Cvvs = Ru/Ms[iSpecies]*(thetav[iSpecies]*thetav[iSpecies]/(Tve*Tve))*
            exp(thetav[iSpecies]/Tve)/((exp(thetav[iSpecies]/Tve)-1)*
                                       (exp(thetav[iSpecies]/Tve)-1));
        dCvvs = (-2/Tve - thetav[iSpecies]/(Tve*Tve) +
                 2*thetav[iSpecies]*exp(thetav[iSpecies]/Tve)/
                 (Tve*Tve*(exp(thetav[iSpecies]/Tve)-1)))*Cvvs;
      } else {
        Cvvs = 0.0;
        dCvvs = 0.0;
      }


      // Electronic energy specific heat
      An1 = 0.0;
      Bn1 = 0.0;
      Bn2 = g[iSpecies][0]*thetae[iSpecies][0]/(Tve*Tve)*exp(-thetae[iSpecies][0]/Tve);
      Bn3 = g[iSpecies][0]*(thetae[iSpecies][0]*thetae[iSpecies][0]/
          (Tve*Tve*Tve*Tve))*exp(-thetae[iSpecies][0]/Tve);
      Bn4 = 0.0;
      Bd1 = g[iSpecies][0]*exp(-thetae[iSpecies][0]/Tve);
      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
        thoTve = thetae[iSpecies][iEl]/Tve;
        exptv = exp(-thetae[iSpecies][iEl]/Tve);

        An1 += g[iSpecies][iEl]*thoTve*thoTve*exptv;
        Bn1 += g[iSpecies][iEl]*thetae[iSpecies][iEl]*exptv;
        Bn2 += g[iSpecies][iEl]*thoTve/Tve*exptv;
        Bn3 += g[iSpecies][iEl]*thoTve*thoTve/(Tve*Tve)*exptv;
        Bn4 += g[iSpecies][iEl]*thoTve*thoTve*thoTve/Tve*exptv;
        Bd1 += g[iSpecies][iEl]*exptv;
      }
      A = An1/Bd1;
      B = Bn1*Bn2/(Bd1*Bd1);
      Cves = Ru/Ms[iSpecies]*(A-B);

      dCves = Ru/Ms[iSpecies]*(-2.0/Tve*(A-B) - 2*Bn2/Bd1*(A-B) -
                               Bn1*Bn3/(Bd1*Bd1) + Bn4/Bd1);

      GradV[RHOCVVE_INDEX][iDim] += V[RHOS_INDEX+iSpecies]*(dCvvs+dCves)*GradV[TVE_INDEX][iDim] +
          GradV[RHOS_INDEX+iSpecies][iDim]*(Cvvs+Cves);

    }
  }

  delete [] Grhou2;
  delete [] rhou;
  return false;
}

void CTNE2EulerVariable::SetPrimVar_Gradient(CConfig *config, CFluidModel *FluidModel) {

  /*--- Use built-in method on TNE2 variable global types ---*/
  GradCons2GradPrimVar(config, FluidModel, Solution, Primitive,
                       Gradient, Gradient_Primitive);

}
