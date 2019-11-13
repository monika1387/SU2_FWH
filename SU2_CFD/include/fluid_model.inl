/*!
 * \file fluid_model.inl
 * \brief In-Line subroutines of the <i>solver_structure.hpp</i> file.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
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

inline void CFluidModel::SetCpModel(CConfig *config) { }

inline su2double CFluidModel::GetPressure ()     { return Pressure; }
inline su2double CFluidModel::GetSoundSpeed ()   { return sqrt(SoundSpeed2); }
inline su2double CFluidModel::GetSoundSpeed2 ()  { return SoundSpeed2; }
inline su2double CFluidModel::GetDensity ()      { return Density; }
inline su2double CFluidModel::GetEntropy ()      { return Entropy; }
inline su2double CFluidModel::GetStaticEnergy () { return StaticEnergy; }
inline su2double CFluidModel::GetTemperature ()  { return Temperature; }
inline su2double CFluidModel::GetCp ()           { return Cp; }
inline su2double CFluidModel::GetCv ()           { return Cv; }
inline su2double CFluidModel::GetdPdrho_e () { return dPdrho_e; }
inline su2double CFluidModel::GetdPde_rho () { return dPde_rho; }
inline su2double CFluidModel::GetdTdrho_e () { return dTdrho_e; }
inline su2double CFluidModel::GetdTde_rho () { return dTde_rho; }
inline su2double CFluidModel::Getdhdrho_P () { return dhdrho_P; }
inline su2double CFluidModel::GetdhdP_rho () { return dhdP_rho; }
inline su2double CFluidModel::Getdsdrho_P () { return dsdrho_P; }
inline su2double CFluidModel::GetdsdP_rho () { return dsdP_rho; }

inline su2double CFluidModel::GetLaminarViscosity () {
        LaminarViscosity->SetViscosity(Temperature, Density);
        Mu = LaminarViscosity->GetViscosity();
        LaminarViscosity->SetDerViscosity(Temperature, Density);
        dmudrho_T= LaminarViscosity->Getdmudrho_T();
        dmudT_rho= LaminarViscosity->GetdmudT_rho();
        return Mu;
}

inline void CFluidModel::Calc_TransportCoeff(CConfig *config, su2double *V) {}

inline su2double CFluidModel::Getdmudrho_T () {
        return LaminarViscosity->Getdmudrho_T();
}

inline su2double CFluidModel::GetdmudT_rho () {
        return LaminarViscosity->GetdmudT_rho();
}

inline su2double CFluidModel::GetThermalConductivity () {
        ThermalConductivity->SetConductivity(Temperature, Density, Mu, Mu_Turb, Cp);
        Kt = ThermalConductivity->GetConductivity();
        ThermalConductivity->SetDerConductivity(Temperature, Density, dmudrho_T, dmudT_rho, Cp);
        dktdrho_T= ThermalConductivity->Getdktdrho_T();
        dktdT_rho= ThermalConductivity->GetdktdT_rho();
        return Kt;
}

inline su2double CFluidModel::Getdktdrho_T () {
        return dktdrho_T;
}

inline su2double CFluidModel::GetdktdT_rho () {
        return dktdT_rho;
}

inline void CFluidModel::SetTDState_rhoe (su2double rho, su2double e ) { }
inline void CFluidModel::SetTDState_PT (su2double P, su2double T ) { }
inline void CFluidModel::SetTDState_Prho (su2double P, su2double rho ) { }
inline void CFluidModel::SetTDState_hs (su2double h, su2double s ) { }
inline void CFluidModel::SetTDState_rhoT (su2double rho, su2double T ) { }
inline void CFluidModel::SetEnergy_Prho (su2double P, su2double rho ) { }
inline void CFluidModel::SetTDState_Ps (su2double P, su2double s ) { }
inline void CFluidModel::ComputeDerivativeNRBC_Prho (su2double P, su2double rho ){ }
inline void CFluidModel::SetTDState_T (su2double val_Temperature) { }
inline void CFluidModel::SetEddyViscosity (su2double val_Mu_Turb) { Mu_Turb = val_Mu_Turb; }
inline void CFluidModel::InitializeMixture(CConfig *config) { cout << "You aint initializing nothing pal" << endl;}
inline su2double* CFluidModel::GetMassFrac_FreeStream(void) {}

inline unsigned short CFluidModel::GetnReactions(void) { }
inline int ***CFluidModel::GetReaction_Map(void) { }
inline su2double ***CFluidModel::GetCollisionIntegral00(void) {  }
inline su2double ***CFluidModel::GetCollisionIntegral11(void) {  }
inline su2double **CFluidModel::GetBlottnerCoeff(void) {  }
inline bool CFluidModel::GetIonization (){}
inline su2double CFluidModel::GetArrheniusCoeff(unsigned short iReaction) { }
inline su2double CFluidModel::GetArrheniusEta(unsigned short iReaction) { }
inline su2double CFluidModel::GetArrheniusTheta(unsigned short iReaction) {}
inline su2double* CFluidModel::GetRxnTcf_a(void) { }
inline su2double* CFluidModel::GetRxnTcf_b(void) { }
inline su2double* CFluidModel::GetRxnTcb_a(void) { }
inline su2double* CFluidModel::GetRxnTcb_b(void) { }
inline su2double* CFluidModel::GetDissociationPot(void) { }
inline su2double CFluidModel::GetCharVibTemp(unsigned short iSpecies) { }
inline su2double* CFluidModel::GetCharVibTemp() { }
inline su2double** CFluidModel::GetCharElTemp() { }
inline unsigned short* CFluidModel::GetnElStates() { }
inline su2double** CFluidModel::GetElDegeneracy() { }
inline su2double* CFluidModel::GetRotationModes() {  }
inline su2double* CFluidModel::GetRefTemperature() { }
inline su2double* CFluidModel::GetWall_Catalycity() { }
inline su2double* CFluidModel::GetMolar_Mass() { }
inline su2double CFluidModel::GetMolar_Mass(unsigned short iSpecies) { }
inline su2double* CFluidModel::GetEnthalpy_Formation(void) {  }
inline su2double CFluidModel::GetEnthalpy_Formation(unsigned short iSpecies) { }
inline void CFluidModel::GetChemistryEquilConstants(CConfig *config, su2double **RxnConstantTable, unsigned short iReaction){}
inline su2double CFluidModel::CalcEve(su2double val_Tve, unsigned short val_Species) { return 0.0; }

inline su2double CFluidModel::Calc_CvVibElSpecies(su2double val_Tve, unsigned short val_Species){return 0.0;}

inline su2double CFluidModel::Calc_CvTraRotSpecies(su2double Ru, unsigned short val_Species) { return 0.0; }

inline su2double CFluidModel::Calc_MixtureEnergy(su2double* cs, su2double sqvel,
                                    su2double rho,
                                    su2double T, su2double Tve) { return 0.0; }

inline su2double CFluidModel::Calc_SoundSpeed(su2double *cs, su2double rhoCvtr,
                                              su2double rho, su2double P) { return 0.0; }

inline su2double CFluidModel::Calc_Enthalpies(su2double val_T, su2double val_eves,
                                              unsigned short val_Species) { return 0.0; }

inline su2double CFluidModel::Calc_Density(su2double *MassFrac,
                                           su2double T, su2double Tve, su2double P) { return 0.0; }

inline void CFluidModel::CalcdPdU(su2double *V, su2double *val_eves, su2double *val_dPdU) { }
inline void CFluidModel::CalcdTdU(su2double *V, su2double *val_dTdU) { }
inline void CFluidModel::CalcdTvedU(su2double *V, su2double *val_eves, su2double *val_dTvedU) {}

//DELETE ME, ADDRESS THS
//inline void CFluidModel::InitializeMixture(CConfig *config) {}
//inline bool CFluidModel::Get_Ionization() { return ionization; }
//inline vector<su2double> CFluidModel::Get_MolarMass(){ return Ms; }
//inline vector<su2double> CFluidModel::Get_CvTraRotSpecies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline vector<su2double> CFluidModel::Get_CvVibElSpecies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline su2double  CFluidModel::Get_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve){}
//inline su2double  CFluidModel::Get_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve){}
//inline su2double  CFluidModel::Get_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve){}
//inline su2double  CFluidModel::Get_MixtureEnergy(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline vector<su2double> CFluidModel::Get_MixtureEnergies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline vector<su2double> CFluidModel::Get_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {}
//inline su2double* CFluidModel::Get_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline vector<su2double> CFluidModel::Get_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline su2double  CFluidModel::Get_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline vector<su2double> CFluidModel::Get_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline su2double* CFluidModel::Get_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline su2double* CFluidModel::Get_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline su2double  CFluidModel::Get_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline vector<su2double> CFluidModel::Get_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline vector<su2double> CFluidModel::Get_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve) {}
//inline su2double  CFluidModel::Get_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve) {}
//inline su2double  CFluidModel::Get_Density(su2double T, su2double *Xs, su2double P){}
