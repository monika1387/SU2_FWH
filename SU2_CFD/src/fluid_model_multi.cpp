/*!
* fluid_model_multi.cpp
* \brief Source of the two-temperature model
* \author W. Maier, C. Garbacz
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

#include "../include/fluid_model.hpp"

CMultiSpeciesGas::CMultiSpeciesGas(unsigned short val_nSpecies,
                                   unsigned short val_nDim,
                                   bool val_ionization)  : CFluidModel() {

  nSpecies      = val_nSpecies;
  nDim          = val_nDim;
  ionization    = val_ionization;
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
}

CMultiSpeciesGas::~CMultiSpeciesGas(void) { }


CMutationGas::CMutationGas(string Mixfile, string Transport, unsigned short val_nSpecies,
                           unsigned short val_nDim, bool val_ionization) :
                           CMultiSpeciesGas(val_nSpecies, val_nDim, val_ionization)//,
                           /*opt(Mixfile)*/{

  //opt.setMechanism(Mixfile);                //DELETE ME,THESE DONT WORK YET
  //opt.setStateModel("ChemNonEqTTv");
  //opt.setViscosityAlgorithm(Transport);
  //opt.setThermalConductivityAlgorithm(Transport);
}

CMutationGas::~CMutationGas(void) {

   delete [] mutation;
   delete [] Ds;
   delete [] hs;

 }

void CMutationGas::InitializeMixture(CConfig *config) {

//  comp = new su2double[nSpecies];
//  hs   = new su2double[nSpecies];
//  Ds   = new su2double[nSpecies];

//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    comp[iSpecies]=config->GetInitial_Gas_Composition(iSpecies);

//  mutation = new CMutation(comp, nSpecies/*, opt*/); //DELETE ME, Opt Not working
//  Ms       = mutation->Mutation_MolarMass();

}

vector<su2double> CMutationGas::Calc_CvTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {

//  Cv_ks = mutation->Mutation_Get_CvModeSpecies(cs, rho, T, Tve);
//  Cv_trs.resize(nSpecies);

//  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    Cv_trs[iSpecies] = Cv_ks[iSpecies];

//  return Cv_trs;

}

vector<su2double> CMutationGas::Calc_CvVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {

//  Cv_ks = mutation->Mutation_Get_CvModeSpecies(cs, rho, T, Tve);
//  Cv_ves.resize(nSpecies);

//  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    Cv_ves[iSpecies] = Cv_ks[nSpecies+iSpecies];

//  return Cv_ves;

}

su2double* CMutationGas::Calc_CpTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {

//  Cp_ks  = mutation->Mutation_Get_CpModeSpecies(cs, rho, T, Tve);
//  Cp_trs = new su2double[nSpecies];

//  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    Cp_trs[iSpecies] = Cp_ks[iSpecies];

//  return Cp_trs;

}

su2double* CMutationGas::Calc_CpVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {

//  Cp_ks = mutation->Mutation_Get_CpModeSpecies(cs, rho, T, Tve);
//  Cp_ves = new su2double[nSpecies];

//  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    Cp_ves[iSpecies] = Cp_ks[nSpecies+iSpecies];

//  return Cp_ves;

}

su2double CMutationGas::Calc_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve){

//  su2double Cp, Cv, *Cp_trs, *Cp_ves, Cvtr, Cvve, Cptr, Cpve;
//  vector<su2double> Cv_trs, Cv_ves;

//  Cp_trs = Calc_CpTraRotSpecies(cs, rho, T, Tve);
//  Cp_ves = Calc_CpVibElSpecies(cs, rho, T, Tve);
//  Cv_trs = Calc_CvTraRotSpecies(cs, rho, T, Tve);
//  Cv_ves = Calc_CvVibElSpecies(cs, rho, T, Tve);

//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    Cptr += cs[iSpecies] * Cp_trs[iSpecies];
//    Cpve += cs[iSpecies] * Cp_ves[iSpecies];
//    Cvtr += cs[iSpecies] * Cv_trs[iSpecies];
//    Cvve += cs[iSpecies] * Cv_ves[iSpecies];
//  }

//  Cp = Cptr + Cpve;
//  Cv = Cvtr + Cvve;

//  gamma = 1.4; // WHY IS THIS SET TO 1.4 -> ASK CATARINA, DELETE ME

//  return gamma;

}

su2double CMutationGas::Calc_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){

//  gammaFrozen = mutation->Mutation_Get_GammaFrozen(cs, rho, T, Tve);
//  return gammaFrozen;

}

su2double CMutationGas::Calc_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve){

//  gammaEquilibrium = mutation->Mutation_Get_GammaEquilibrium(cs, rho, T, Tve);
//  return gammaEquilibrium;

}

su2double CMutationGas::Calc_MixtureEnergy(su2double* cs, su2double rho, su2double T, su2double Tve) {

//  E = mutation->Mutation_Get_MixtureEnergy(cs, rho, T, Tve);
//  return E;

}

vector<su2double> CMutationGas::Calc_MixtureEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {

//  Energies = mutation->Mutation_Get_MixtureEnergies(cs, rho, T, Tve);
//  return Energies;

}

vector<su2double> CMutationGas::Calc_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {

//  Energies_Species = mutation->Mutation_Get_SpeciesEnergies(cs, rho, T, Tve);
//  return Energies_Species;

}

su2double* CMutationGas::Calc_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve) {

//  Ws = new su2double[nSpecies];
//  Ws = mutation->Mutation_Get_NetProductionRates(cs, rho, T, Tve);
//  return Ws;

}

vector<su2double> CMutationGas::Calc_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve) {

//  OmegaVT = mutation->Mutation_Get_VTEnergysourceTerm(cs, rho, T, Tve);
//  return OmegaVT;

}

su2double CMutationGas::Calc_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve) {

//  Tref = mutation->Mutation_Get_ReferenceTemperature(cs, rho, T, Tve);
//  return Tref;

}

vector<su2double> CMutationGas::Calc_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve) {

//  hf = mutation->Mutation_Get_EnthalpiesFormation(cs, rho, T, Tve);
//  return hf;

}

su2double* CMutationGas::Calc_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve) {

//  hs = mutation->Mutation_Get_Enthalpies(cs, rho, T, Tve);
//  return hs;

}

su2double* CMutationGas::Calc_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve) {

//   Ds = mutation->Mutation_Get_DiffusionCoeff(cs, rho, T, Tve);
//   return Ds;

}

su2double  CMutationGas::Calc_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve) {

//   mu = mutation->Mutation_Get_Viscosity(cs, rho, T, Tve);
//   return mu;

}

vector<su2double> CMutationGas::Calc_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve) {

//   lambda = mutation->Mutation_Get_ThermalConductivity(cs, rho, T, Tve);
//   return lambda;

}

vector<su2double> CMutationGas::Calc_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve){

//  Temp = mutation->Mutation_Get_Temperatures(cs, rho, rhoE, rhoEve);
//  return Temp;

}

su2double CMutationGas::Calc_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){

//  a = mutation->Mutation_Get_SoundSpeedFrozen(cs, rho, T, Tve);
//  return a;

}

su2double CMutationGas::Calc_Density(su2double T, su2double *Xs, su2double P){

//  Density = mutation->Mutation_Get_Density(T, Xs, P);
//  return Density;

}

void CMultiSpeciesGas::CalcdPdU(su2double *V, su2double *val_eves,
                             CConfig *config, su2double *val_dPdU) {

  // Note: Requires SetDensity(), SetTemperature(), SetPressure(), & SetGasProperties()
  // Note: Electron energy not included properly.

  unsigned short iDim, iSpecies, iEl, nHeavy, nEl, *nElStates;
  su2double *Ms, *Tref, *hf, *xi, *thetav, **thetae, **g;
  su2double RuSI, Ru, RuBAR, CvtrBAR, rhoCvtr, rhoCvve, Cvtrs, rho_el, sqvel, conc;
  su2double rho, rhos, T, Tve, ef;
  su2double num, denom;

  /*--- Determine the number of heavy species ---*/
  if (ionization) {
    nHeavy = nSpecies-1;
    nEl    = 1;
    rho_el = V[RHOS_INDEX+nSpecies-1];
  } else {
      nHeavy = nSpecies;
      nEl    = 0;
      rho_el = 0.0;
    }

  /*--- Read gas mixture properties from config ---*/
  Ms        = config->GetMolar_Mass();
  Tref      = config->GetRefTemperature();
  hf        = config->GetEnthalpy_Formation();
  xi        = config->GetRotationModes();
  thetav    = config->GetCharVibTemp();
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  nElStates = config->GetnElStates();

  /*--- Rename for convenience ---*/
  RuSI    = UNIVERSAL_GAS_CONSTANT;
  Ru      = 1000.0*RuSI;
  T       = V[T_INDEX];
  Tve     = V[TVE_INDEX];
  rho     = V[RHO_INDEX];
  rhoCvtr = V[RHOCVTR_INDEX];
  rhoCvve = V[RHOCVVE_INDEX];

  /*--- Pre-compute useful quantities ---*/
  RuBAR   = 0.0;
  CvtrBAR = 0.0;
  sqvel   = 0.0;
  conc    = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += V[VEL_INDEX+iDim] * V[VEL_INDEX+iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    CvtrBAR += V[RHOS_INDEX+iSpecies]*(3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
    conc    += V[RHOS_INDEX+iSpecies]/Ms[iSpecies];
  }

  // Species density
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    rhos  = V[RHOS_INDEX+iSpecies];
    ef    = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
    Cvtrs = (3.0/2.0+xi[iSpecies]/2.0)*Ru/Ms[iSpecies];

    val_dPdU[iSpecies] =  T*Ru/Ms[iSpecies] + Ru*conc/rhoCvtr *
    (-Cvtrs*(T-Tref[iSpecies]) -
    ef + 0.5*sqvel);
  }
  if (ionization) {
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      //      evibs = Ru/Ms[iSpecies] * thetav[iSpecies]/(exp(thetav[iSpecies]/Tve)-1.0);
      //      num = 0.0;
      //      denom = g[iSpecies][0] * exp(-thetae[iSpecies][0]/Tve);
      //      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      //        num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      //        denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      //      }
      //      eels = Ru/Ms[iSpecies] * (num/denom);
      val_dPdU[iSpecies] -= rho_el * Ru/Ms[nSpecies-1] * (val_eves[iSpecies])/rhoCvve;
    }
    ef = hf[nSpecies-1] - Ru/Ms[nSpecies-1]*Tref[nSpecies-1];
    val_dPdU[nSpecies-1] = Ru*conc/rhoCvtr * (-ef + 0.5*sqvel)
      + Ru/Ms[nSpecies-1]*Tve
      - rho_el*Ru/Ms[nSpecies-1] * (-3.0/2.0*Ru/Ms[nSpecies-1]*Tve)/rhoCvve;
  }
  // Momentum
  for (iDim = 0; iDim < nDim; iDim++)
    val_dPdU[nSpecies+iDim] = -conc*Ru*V[VEL_INDEX+iDim]/rhoCvtr;

  // Total energy
  val_dPdU[nSpecies+nDim]   = conc*Ru / rhoCvtr;

  // Vib.-el energy
  val_dPdU[nSpecies+nDim+1] = -val_dPdU[nSpecies+nDim]
    + rho_el*Ru/Ms[nSpecies-1]*1.0/rhoCvve;

}

su2double CMultiSpeciesGas::CalcEve(CConfig *config, su2double val_Tve,
                            unsigned short val_Species) {

  unsigned short iEl, *nElStates;
  su2double *Ms, *thetav, **thetae, **g, *hf, *Tref, RuSI, Ru;
  su2double Tve, Ev, Eel, Ef;
  su2double num, denom;

  /*--- Read gas mixture properties from config ---*/
  Ms        = config->GetMolar_Mass();

  /*--- Rename for convenience ---*/
  RuSI  = UNIVERSAL_GAS_CONSTANT;
  Ru    = 1000.0*RuSI;
  Tve   = val_Tve;

  /*--- Electron species energy ---*/
  if ((ionization) && (val_Species == nSpecies-1)) {

    /*--- Get quantities from CConfig ---*/
    Tref = config->GetRefTemperature();
    hf   = config->GetEnthalpy_Formation();

    /*--- Calculate formation energy ---*/
    Ef = hf[val_Species] - Ru/Ms[val_Species] * Tref[val_Species];

    /*--- Electron t-r mode contributes to mixture vib-el energy ---*/
    Eel = (3.0/2.0) * Ru/Ms[val_Species] * (Tve - Tref[val_Species]) + Ef;
    Ev  = 0.0;

  }

  /*--- Heavy particle energy ---*/
  else {

    /*--- Read from CConfig ---*/
    thetav    = config->GetCharVibTemp();
    thetae    = config->GetCharElTemp();
    g         = config->GetElDegeneracy();
    nElStates = config->GetnElStates();

    /*--- Calculate vibrational energy (harmonic-oscillator model) ---*/
    if (thetav[val_Species] != 0.0)
    Ev = Ru/Ms[val_Species] * thetav[val_Species] / (exp(thetav[val_Species]/Tve)-1.0);
    else
    Ev = 0.0;

    /*--- Calculate electronic energy ---*/
    num = 0.0;
    denom = g[val_Species][0] * exp(-thetae[val_Species][0]/Tve);
    for (iEl = 1; iEl < nElStates[val_Species]; iEl++) {
      num   += g[val_Species][iEl] * thetae[val_Species][iEl] * exp(-thetae[val_Species][iEl]/Tve);
      denom += g[val_Species][iEl] * exp(-thetae[val_Species][iEl]/Tve);
    }
    Eel = Ru/Ms[val_Species] * (num/denom);
  }

  return Ev + Eel;
}

su2double CMultiSpeciesGas::CalcHs(CConfig *config, su2double val_T,
                           su2double val_eves, unsigned short val_Species) {

  su2double RuSI, Ru, *xi, *Ms, *hf, *Tref, T, eve, ef, hs;

  /*--- Read from config ---*/
  xi   = config->GetRotationModes();
  Ms   = config->GetMolar_Mass();
  hf   = config->GetEnthalpy_Formation();
  Tref = config->GetRefTemperature();

  /*--- Rename for convenience ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  T    = val_T;

  /*--- Calculate vibrational-electronic energy per unit mass ---*/
  eve = val_eves;

  /*--- Calculate formation energy ---*/
  ef = hf[val_Species] - Ru/Ms[val_Species]*Tref[val_Species];

  hs = Ru/Ms[val_Species]*T
      + (3.0/2.0+xi[val_Species]/2.0)*Ru/Ms[val_Species]*T
      + hf[val_Species] + eve;

  return hs;
}

su2double CMultiSpeciesGas::CalcCvve(su2double val_Tve, CConfig *config, unsigned short val_Species) {

  unsigned short iEl, *nElStates;
  su2double *Ms, *thetav, **thetae, **g, RuSI, Ru;
  su2double thoTve, exptv, thsqr, Cvvs, Cves;
  su2double Tve;
  su2double num, num2, num3, denom;

  /*--- Read from config ---*/
  thetav    = config->GetCharVibTemp();
  thetae    = config->GetCharElTemp();
  g         = config->GetElDegeneracy();
  Ms        = config->GetMolar_Mass();
  nElStates = config->GetnElStates();

  /*--- Rename for convenience ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  Tve  = val_Tve;

  /*--- If requesting electron specific heat ---*/
  if (ionization && val_Species == nSpecies-1) {
    Cvvs = 0.0;
    Cves = 3.0/2.0 * Ru/Ms[nSpecies-1];
  }

  /*--- Heavy particle specific heat ---*/
  else {

    /*--- Vibrational energy ---*/
    if (thetav[val_Species] != 0.0) {
      thoTve = thetav[val_Species]/Tve;
      exptv = exp(thetav[val_Species]/Tve);
      thsqr = thetav[val_Species]*thetav[val_Species];
      Cvvs  = Ru/Ms[val_Species] * thoTve*thoTve * exptv / ((exptv-1.0)*(exptv-1.0));
    } else {
      Cvvs = 0.0;
    }

    /*--- Electronic energy ---*/
    if (nElStates[val_Species] != 0) {
      num = 0.0; num2 = 0.0;
      denom = g[val_Species][0] * exp(-thetae[val_Species][0]/Tve);
      num3  = g[val_Species][0] * (thetae[val_Species][0]/(Tve*Tve))*exp(-thetae[val_Species][0]/Tve);
      for (iEl = 1; iEl < nElStates[val_Species]; iEl++) {
        thoTve = thetae[val_Species][iEl]/Tve;
        exptv = exp(-thetae[val_Species][iEl]/Tve);

        num   += g[val_Species][iEl] * thetae[val_Species][iEl] * exptv;
        denom += g[val_Species][iEl] * exptv;
        num2  += g[val_Species][iEl] * (thoTve*thoTve) * exptv;
        num3  += g[val_Species][iEl] * thoTve/Tve * exptv;
      }
      Cves = Ru/Ms[val_Species] * (num2/denom - num*num3/(denom*denom));
    } else {
      Cves = 0.0;
    }

  }
  return Cvvs + Cves;
}

void CMultiSpeciesGas::CalcdTdU(su2double *V, CConfig *config,
                             su2double *val_dTdU) {

  unsigned short iDim, iSpecies, nHeavy, nEl;
  su2double *Ms, *xi, *Tref, *hf;
  su2double v2, ef, T, Cvtrs, rhoCvtr, RuSI, Ru;

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Get gas properties from config settings ---*/
  Ms   = config->GetMolar_Mass();
  xi   = config->GetRotationModes();
  hf   = config->GetEnthalpy_Formation();
  Tref = config->GetRefTemperature();

  /*--- Rename for convenience ---*/
  rhoCvtr = V[RHOCVTR_INDEX];
  T       = V[T_INDEX];
  RuSI    = UNIVERSAL_GAS_CONSTANT;
  Ru      = 1000.0*RuSI;

  /*--- Calculate supporting quantities ---*/
  v2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    v2 += V[VEL_INDEX+iDim]*V[VEL_INDEX+iDim];

  /*--- Species density derivatives ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    ef    = hf[iSpecies] - Ru/Ms[iSpecies]*Tref[iSpecies];
    Cvtrs = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
    val_dTdU[iSpecies]   = (-ef + 0.5*v2 + Cvtrs*(Tref[iSpecies]-T)) / rhoCvtr;
  }
  if (ionization) {
    cout << "CTNE2Variable: NEED TO IMPLEMENT dTdU for IONIZED MIX" << endl;
    exit(1);
  }

  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dTdU[nSpecies+iDim] = -V[VEL_INDEX+iDim] / V[RHOCVTR_INDEX];

  /*--- Energy derivatives ---*/
  val_dTdU[nSpecies+nDim]   =  1.0 / V[RHOCVTR_INDEX];
  val_dTdU[nSpecies+nDim+1] = -1.0 / V[RHOCVTR_INDEX];

}

void CMultiSpeciesGas::CalcdTvedU(su2double *V, su2double *val_eves, CConfig *config,
                               su2double *val_dTvedU) {

  unsigned short iDim, iSpecies;
  su2double rhoCvve;

  /*--- Rename for convenience ---*/
  rhoCvve = V[RHOCVVE_INDEX];

  /*--- Species density derivatives ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    val_dTvedU[iSpecies] = -val_eves[iSpecies]/rhoCvve;
  }
  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dTvedU[nSpecies+iDim] = 0.0;

  /*--- Energy derivatives ---*/
  val_dTvedU[nSpecies+nDim]   = 0.0;
  val_dTvedU[nSpecies+nDim+1] = 1.0 / rhoCvve;
}
