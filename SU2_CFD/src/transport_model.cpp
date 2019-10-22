/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna, T. Economon
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

#include "../include/transport_model.hpp"

/*-------------------------------------------------*/
/*--------------- Viscosity Models ----------------*/
/*-------------------------------------------------*/

CViscosityModel::CViscosityModel(void) {

  /*--- Attributes initialization ---*/

  Mu = 0.0;
  dmudrho_T = 0.0;
  dmudT_rho = 0.0;

}

CViscosityModel::~CViscosityModel(void) { }

CConstantViscosity::CConstantViscosity(void) : CViscosityModel() { }

CConstantViscosity::CConstantViscosity(su2double mu_const) : CViscosityModel() {

  /*--- Attributes initialization ---*/

  Mu = mu_const;
  dmudrho_T = 0.0;
  dmudT_rho = 0.0;

}

CConstantViscosity::~CConstantViscosity(void) { }

CSutherland::CSutherland(void) : CViscosityModel() {
  Mu_ref = 0.0;
  T_ref = 0.0;
  S = 0.0;

}

CSutherland::CSutherland(su2double mu_ref, su2double t_ref, su2double s) : CViscosityModel() {

  Mu_ref = mu_ref;
  T_ref = t_ref;
  S = s;
}

CSutherland::~CSutherland(void) { }

void CSutherland::SetViscosity(su2double T, su2double rho) {

  const su2double TnonDim = T/T_ref;
  Mu = Mu_ref*TnonDim*sqrt(TnonDim)*((T_ref + S)/(T + S));

}

void CSutherland::SetDerViscosity(su2double T, su2double rho) {

  dmudrho_T = 0.0;

  const su2double T_refInv = 1.0/T_ref;
  const su2double TnonDim  = T_refInv*T;
  const su2double TSInv    = 1.0/(T + S);
 
  dmudT_rho = Mu_ref*(T_ref + S)*TSInv*sqrt(TnonDim)
            * (1.5*T_refInv - TnonDim*TSInv);

}

CPolynomialViscosity::CPolynomialViscosity(void) : CViscosityModel() {
  nPolyCoeffs = 0;
  b           = NULL;
}

CPolynomialViscosity::CPolynomialViscosity(unsigned short val_nCoeffs, su2double* val_b) : CViscosityModel() {
  
  /*--- Attributes initialization ---*/
  
  nPolyCoeffs = val_nCoeffs;
  b = new su2double[nPolyCoeffs];
  
  for (unsigned short iVar = 0; iVar < nPolyCoeffs; iVar++)
    b[iVar] = val_b[iVar];
  
}

CPolynomialViscosity::~CPolynomialViscosity(void) {
  if (b != NULL) delete [] b;
}

void CPolynomialViscosity::SetViscosity(su2double T, su2double rho) {
  
  /*--- Evaluate the new Mu from the coefficients and temperature. ---*/
  
  Mu = b[0];
  for (unsigned short iVar = 1; iVar < nPolyCoeffs; iVar++)
    Mu += b[iVar]*pow(T,iVar);
  
}

/*-------------------------------------------------*/
/*---------- Thermal Conductivity Models ----------*/
/*-------------------------------------------------*/

CConductivityModel::CConductivityModel(void) {

  /*--- Attributes initialization ---*/

  Kt = 0.0;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;

}

CConductivityModel::~CConductivityModel(void) { }

CConstantConductivity::CConstantConductivity(void) : CConductivityModel() { }

CConstantConductivity::CConstantConductivity(su2double kt_const) : CConductivityModel() {

  /*--- Attributes initialization ---*/

  Kt = kt_const;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;

}

CConstantConductivity::~CConstantConductivity(void) { }

CConstantConductivityRANS::CConstantConductivityRANS(void) : CConductivityModel() { }

CConstantConductivityRANS::CConstantConductivityRANS(su2double kt_const, su2double pr_turb) : CConductivityModel() {
  
  /*--- Attributes initialization ---*/
  
  Kt_Lam    = kt_const;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;
  
  Prandtl_Turb = pr_turb;

}

void CConstantConductivityRANS::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {
  
  Kt = Kt_Lam + cp*mu_turb/Prandtl_Turb;
  
}

CConstantConductivityRANS::~CConstantConductivityRANS(void) { }

CConstantPrandtl::CConstantPrandtl(void) : CConductivityModel() { }

CConstantPrandtl::CConstantPrandtl(su2double pr_const) : CConductivityModel() {

  /*--- Attributes initialization ---*/

  Pr_const = pr_const;

}

void CConstantPrandtl::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {

  Kt = mu_lam*cp/Pr_const;

}

void CConstantPrandtl::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {

  dktdrho_T = dmudrho_T*cp/Pr_const;
  dktdT_rho = dmudT_rho*cp/Pr_const;

}

CConstantPrandtl::~CConstantPrandtl(void) { }

CConstantPrandtlRANS::CConstantPrandtlRANS(void) : CConductivityModel() { }

CConstantPrandtlRANS::CConstantPrandtlRANS(su2double pr_lam, su2double pr_turb) : CConductivityModel() {

  /*--- Attributes initialization ---*/

  Prandtl_Lam  = pr_lam;
  Prandtl_Turb = pr_turb;
}

void CConstantPrandtlRANS::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {

  Kt = cp * ((mu_lam/Prandtl_Lam) + (mu_turb/Prandtl_Turb));
  
}

CConstantPrandtlRANS::~CConstantPrandtlRANS(void) { }

CPolynomialConductivity::CPolynomialConductivity(void) : CConductivityModel() {
  nPolyCoeffs = 0;
  b           = NULL;
}

CPolynomialConductivity::CPolynomialConductivity(unsigned short val_nCoeffs, su2double* val_b) : CConductivityModel() {
  
  /*--- Attributes initialization ---*/
  
  nPolyCoeffs = val_nCoeffs;
  b = new su2double[nPolyCoeffs];
  
  for (unsigned short iVar = 0; iVar < nPolyCoeffs; iVar++)
    b[iVar] = val_b[iVar];
  
}

CPolynomialConductivity::~CPolynomialConductivity(void) {
  if (b != NULL) delete [] b;
}

void CPolynomialConductivity::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {
  
  /*--- Evaluate the new Kt from the coefficients and temperature. ---*/
  
  Kt = b[0];
  for (unsigned short iVar = 1; iVar < nPolyCoeffs; iVar++)
    Kt += b[iVar]*pow(T,iVar);
  
}

CPolynomialConductivityRANS::CPolynomialConductivityRANS(unsigned short val_nCoeffs, su2double* val_b, su2double pr_turb) : CConductivityModel() {
  
  /*--- Attributes initialization ---*/
  
  nPolyCoeffs = val_nCoeffs;
  b = new su2double[nPolyCoeffs];
  
  for (unsigned short iVar = 0; iVar < nPolyCoeffs; iVar++)
    b[iVar] = val_b[iVar];
  
  Prandtl_Turb = pr_turb;
  
}

CPolynomialConductivityRANS::~CPolynomialConductivityRANS(void) {
  if (b != NULL) delete [] b;
}

void CPolynomialConductivityRANS::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {
  
  /*--- Evaluate the new Kt from the coefficients and temperature. ---*/
  
  Kt = b[0];
  for (unsigned short iVar = 1; iVar < nPolyCoeffs; iVar++)
    Kt += b[iVar]*pow(T,iVar);
  
  /*--- Add a component due to turbulence to compute effective conductivity. ---*/
  
  Kt += cp*mu_turb/Prandtl_Turb;
  
}

/*-------------------------------------------------*/
/*--------- MultiSpecies Transport Models ---------*/
/*-------------------------------------------------*/
//CTransportModel::CTransportModel(){}

//CTransportModel::~CTransportModel(void){}

//CWilkeBlottEucken::CWilkeBlottEucken(): CTransportModel() {

//  unsigned short iSpecies, jSpecies;
//  su2double *Ms, Mi, Mj, M;
//  su2double rho, T, Tve, RuSI, Ru, *xi;
//  su2double Xs[nSpecies], conc;
//  su2double Cves;
//  su2double phis[nSpecies], mus[nSpecies], ks[nSpecies], kves[nSpecies];
//  su2double denom, tmp1, tmp2;
//  su2double **Blottner;
//  su2double ***Omega00, Omega_ij;

//  /*--- Rename for convenience ---*/
//  rho  = Primitive[RHO_INDEX];
//  T    = Primitive[T_INDEX];
//  Tve  = Primitive[TVE_INDEX];
//  Ms   = FluidModel->GetMolar_Mass();
//  xi   = FluidModel->GetRotationModes();
//  RuSI = UNIVERSAL_GAS_CONSTANT;
//  Ru   = 1000.0*RuSI;

//  /*--- Acquire collision integral information ---*/
//  Omega00 = FluidModel->GetCollisionIntegral00();

//  /*--- Calculate species mole fraction ---*/
//  conc = 0.0;
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    Xs[iSpecies] = Primitive[RHOS_INDEX+iSpecies]/Ms[iSpecies];
//    conc        += Xs[iSpecies];
//  }
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    Xs[iSpecies] = Xs[iSpecies]/conc;

//  /*--- Calculate mixture molar mass (kg/mol) ---*/
//  // Note: Species molar masses stored as kg/kmol, need 1E-3 conversion
//  M = 0.0;
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    M += Ms[iSpecies]*Xs[iSpecies];
//  M = M*1E-3;


//  /*---+++                  +++---*/
//  /*--- Diffusion coefficients ---*/
//  /*---+++                  +++---*/

//  /*--- Solve for binary diffusion coefficients ---*/
//  // Note: Dij = Dji, so only loop through req'd indices
//  // Note: Correlation requires kg/mol, hence 1E-3 conversion from kg/kmol
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    Mi = Ms[iSpecies]*1E-3;
//    for (jSpecies = iSpecies; jSpecies < nSpecies; jSpecies++) {
//      Mj = Ms[jSpecies]*1E-3;

//      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
//      Omega_ij = 1E-20/PI_NUMBER * Omega00[iSpecies][jSpecies][3]
//          * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
//          +  Omega00[iSpecies][jSpecies][1]*log(T)
//          +  Omega00[iSpecies][jSpecies][2]);

//      Dij[iSpecies][jSpecies] = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(rho*Omega_ij);
//      Dij[jSpecies][iSpecies] = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(rho*Omega_ij);
//    }
//  }

//  /*--- Calculate species-mixture diffusion coefficient --*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    denom = 0.0;
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      if (jSpecies != iSpecies) {
//        denom += Xs[jSpecies]/Dij[iSpecies][jSpecies];
//      }
//    }
//    DiffusionCoeff[iSpecies] = (1-Xs[iSpecies])/denom;
//    //    DiffusionCoeff[iSpecies] = 0.0;
//  }


//  /*---+++             +++---*/
//  /*--- Laminar viscosity ---*/
//  /*---+++             +++---*/

//  /*--- Get Blottner coefficients ---*/
//  Blottner = FluidModel->GetBlottnerCoeff();

//  /*--- Use Blottner's curve fits for species viscosity ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    mus[iSpecies] = 0.1*exp((Blottner[iSpecies][0]*log(T)  +
//                            Blottner[iSpecies][1])*log(T) +
//        Blottner[iSpecies][2]);

//  /*--- Determine species 'phi' value for Blottner model ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    phis[iSpecies] = 0.0;
//    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//      tmp1 = 1.0 + sqrt(mus[iSpecies]/mus[jSpecies])*pow(Ms[jSpecies]/Ms[iSpecies], 0.25);
//      tmp2 = sqrt(8.0*(1.0+Ms[iSpecies]/Ms[jSpecies]));
//      phis[iSpecies] += Xs[jSpecies]*tmp1*tmp1/tmp2;
//    }
//  }

//  /*--- Calculate mixture laminar viscosity ---*/
//  LaminarViscosity = 0.0;
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
//    LaminarViscosity += Xs[iSpecies]*mus[iSpecies]/phis[iSpecies];


//  /*---+++                +++---*/
//  /*--- Thermal conductivity ---*/
//  /*---+++                +++---*/

//  /*--- Determine species tr & ve conductivities ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    Cves = CalcCvve(Tve, config, iSpecies);
//    ks[iSpecies] = mus[iSpecies]*(15.0/4.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
//    kves[iSpecies] = mus[iSpecies]*Cves;
//  }

//  /*--- Calculate mixture tr & ve conductivities ---*/
//  ThermalCond    = 0.0;
//  ThermalCond_ve = 0.0;
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    ThermalCond    += Xs[iSpecies]*ks[iSpecies]/phis[iSpecies];
//    ThermalCond_ve += Xs[iSpecies]*kves[iSpecies]/phis[iSpecies];
//  }
//}

//CGuptaYos::CGuptaYos(void): CTransportModel() {}


















