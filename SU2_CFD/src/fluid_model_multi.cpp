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
                                   bool val_ionization,
                                   bool val_viscous)  : CFluidModel() {

  nSpecies      = val_nSpecies;
  nDim          = val_nDim;
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

  if (val_viscous){
    LAM_VISC_INDEX   = nSpecies+nDim+8;
    EDDY_VISC_INDEX  = nSpecies+nDim+9;
    DIFF_COEFF_INDEX = nSpecies+nDim+10;
    K_INDEX          = 2*nSpecies+nDim+10;
    KVE_INDEX        = 2*nSpecies+nDim+11;
  }

  ionization    = val_ionization;

  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  DiffusionCoeff  = new su2double[nSpecies];
  Dij = new su2double*[nSpecies];
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Dij[iSpecies] = new su2double[nSpecies];

}

CMultiSpeciesGas::~CMultiSpeciesGas(void) {

  delete [] DiffusionCoeff;
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] Dij[iSpecies];
  delete [] Dij;

}

void CMultiSpeciesGas::InitializeMixture(CConfig *config) {}

CTNE2Gas::CTNE2Gas(unsigned short val_nSpecies,
                   unsigned short val_nDim,
                   bool val_ionization, bool val_viscous) : CMultiSpeciesGas(val_nSpecies,
                                                                             val_nDim,
                                                                             val_ionization,
                                                                             val_viscous) { }

CTNE2Gas::~CTNE2Gas(void) { }

void CTNE2Gas::InitializeMixture(CConfig *config) {

  unsigned short maxEl = 0;
  unsigned short iEl;
  su2double *Gas_Composition;

  Gas_Composition = config->GetGas_Composition();

  switch (config->GetKind_GasModel()) {

  case ONESPECIES:
    /*--- Define parameters of the gas model ---*/
    nSpecies    = 1;
    ionization  = false;

    /*--- Allocate vectors for gas properties ---*/
    Molar_Mass         = new su2double[nSpecies];
    CharVibTemp        = new su2double[nSpecies];
    RotationModes      = new su2double[nSpecies];
    Enthalpy_Formation = new su2double[nSpecies];
    Wall_Catalycity    = new su2double[nSpecies];
    Ref_Temperature    = new su2double[nSpecies];
    nElStates          = new unsigned short[nSpecies];

    MassFrac_FreeStream = new su2double[nSpecies];
    MassFrac_FreeStream[0] = 1.0;

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    // Molar mass [kg/kmol]
    Molar_Mass[0] = 14.0067+15.9994;
    // Characteristic vibrational temperatures for calculating e_vib [K]
    //CharVibTemp[0] = 3395.0;
    CharVibTemp[0] = 1000.0;
    // Formation enthalpy: (JANAF values, [KJ/Kmol])
    Enthalpy_Formation[0] = 0.0;					//N2
    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;

    /*        nElStates[0] = 0;
         CharElTemp   = new double *[nSpecies];
         ElDegeneracy        = new double *[nSpecies];

         OSPthetae    = new double[nElStates[0]];
         OSPthetae[0] = 1.0;
         OSPg         = new double[nElStates[0]];
         OSPg[0]      = 1.0;

         CharElTemp[0] = OSPthetae;
         ElDegeneracy[0] = OSPg;*/

    break;

  case N2:

    /*--- Define parameters of the gas model ---*/
    nReactions  = 2;
    ionization  = false;

    /*--- Allocate vectors for gas properties ---*/
    Wall_Catalycity      = new su2double[nSpecies];
    Molar_Mass           = new su2double[nSpecies];
    CharVibTemp          = new su2double[nSpecies];
    RotationModes        = new su2double[nSpecies];
    Enthalpy_Formation   = new su2double[nSpecies];
    Ref_Temperature      = new su2double[nSpecies];
    Diss                 = new su2double[nSpecies];
    ArrheniusCoefficient = new su2double[nReactions];
    ArrheniusEta         = new su2double[nReactions];
    ArrheniusTheta       = new su2double[nReactions];
    Tcf_a                = new su2double[nReactions];
    Tcf_b                = new su2double[nReactions];
    Tcb_a                = new su2double[nReactions];
    Tcb_b                = new su2double[nReactions];
    nElStates            = new unsigned short[nSpecies];
    Reactions = new int**[nReactions];
    for (unsigned short iRxn = 0; iRxn < nReactions; iRxn++) {
      Reactions[iRxn] = new int*[2];
      for (unsigned short ii = 0; ii < 2; ii++)
        Reactions[iRxn][ii] = new int[6];
    }

    Blottner  = new su2double*[nSpecies];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Blottner[iSpecies] = new su2double[3];

    // Omega[iSpecies][jSpecies][iCoeff]
    Omega00 = new su2double**[nSpecies];
    Omega11 = new su2double**[nSpecies];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Omega00[iSpecies] = new su2double*[nSpecies];
      Omega11[iSpecies] = new su2double*[nSpecies];
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        Omega00[iSpecies][jSpecies] = new su2double[4];
        Omega11[iSpecies][jSpecies] = new su2double[4];
      }
    }

    MassFrac_FreeStream = new su2double[nSpecies];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      MassFrac_FreeStream[iSpecies] = Gas_Composition[iSpecies];

    /*--- Assign gas properties ---*/

    // Wall mass fractions for catalytic boundaries
    Wall_Catalycity[0] = 0.999;
    Wall_Catalycity[1] = 0.001;

    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    RotationModes[1] = 0.0;

    // Molar mass [kg/kmol]
    Molar_Mass[0] = 2.0*14.0067;
    Molar_Mass[1] = 14.0067;

    // Characteristic vibrational temperatures
    CharVibTemp[0] = 3395.0;
    CharVibTemp[1] = 0.0;

    // Formation enthalpy: (JANAF values [KJ/Kmol])
    // J/kg - from Scalabrin
    Enthalpy_Formation[0] = 0.0;					//N2
    Enthalpy_Formation[1] = 3.36E7;		//N

    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;

    // Blottner viscosity coefficients
    // A                        // B                        // C
    Blottner[0][0] = 2.68E-2;   Blottner[0][1] = 3.18E-1;   Blottner[0][2] = -1.13E1;  // N2
    Blottner[1][0] = 1.16E-2;   Blottner[1][1] = 6.03E-1;   Blottner[1][2] = -1.24E1;  // N

    // Number of electron states
    nElStates[0] = 15;                    // N2
    nElStates[1] = 3;                     // N
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate electron data arrays ---*/
    CharElTemp = new su2double*[nSpecies];
    ElDegeneracy      = new su2double*[nSpecies];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      CharElTemp[iSpecies] = new su2double[maxEl];
      ElDegeneracy[iSpecies]      = new su2double[maxEl];
    }

    /*--- Initialize the arrays ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (iEl = 0; iEl < maxEl; iEl++) {
        CharElTemp[iSpecies][iEl] = 0.0;
        ElDegeneracy[iSpecies][iEl] = 0.0;
      }
    }

    /*--- Assign values to data structures ---*/
    // N2: 15 states
    CharElTemp[0][0]  = 0.000000000000000E+00;
    CharElTemp[0][1]  = 7.223156514095200E+04;
    CharElTemp[0][2]  = 8.577862640384000E+04;
    CharElTemp[0][3]  = 8.605026716160000E+04;
    CharElTemp[0][4]  = 9.535118627874400E+04;
    CharElTemp[0][5]  = 9.805635702203200E+04;
    CharElTemp[0][6]  = 9.968267656935200E+04;
    CharElTemp[0][7]  = 1.048976467715200E+05;
    CharElTemp[0][8]  = 1.116489555200000E+05;
    CharElTemp[0][9]  = 1.225836470400000E+05;
    CharElTemp[0][10] = 1.248856873600000E+05;
    CharElTemp[0][11] = 1.282476158188320E+05;
    CharElTemp[0][12] = 1.338060936000000E+05;
    CharElTemp[0][13] = 1.404296391107200E+05;
    CharElTemp[0][14] = 1.504958859200000E+05;
    ElDegeneracy[0][0]  = 1;
    ElDegeneracy[0][1]  = 3;
    ElDegeneracy[0][2]  = 6;
    ElDegeneracy[0][3]  = 6;
    ElDegeneracy[0][4]  = 3;
    ElDegeneracy[0][5]  = 1;
    ElDegeneracy[0][6]  = 2;
    ElDegeneracy[0][7]  = 2;
    ElDegeneracy[0][8]  = 5;
    ElDegeneracy[0][9]  = 1;
    ElDegeneracy[0][10] = 6;
    ElDegeneracy[0][11] = 6;
    ElDegeneracy[0][12] = 10;
    ElDegeneracy[0][13] = 6;
    ElDegeneracy[0][14] = 6;
    // N: 3 states
    CharElTemp[1][0] = 0.000000000000000E+00;
    CharElTemp[1][1] = 2.766469645581980E+04;
    CharElTemp[1][2] = 4.149309313560210E+04;
    ElDegeneracy[1][0] = 4;
    ElDegeneracy[1][1] = 10;
    ElDegeneracy[1][2] = 6;

    /*--- Set Arrhenius coefficients for chemical reactions ---*/
    // Note: Data lists coefficients in (cm^3/mol-s) units, need to convert
    //       to (m^3/kmol-s) to be consistent with the rest of the code
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 3.0E22;
    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    // Characteristic temperature
    ArrheniusTheta[0] = 113200.0;
    ArrheniusTheta[1] = 113200.0;

    /*--- Set reaction maps ---*/
    // N2 + N2 -> 2N + N2
    Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;
    Reactions[0][1][0]=1;		Reactions[0][1][1]=1;		Reactions[0][1][2] =0;
    // N2 + N -> 2N + N
    Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;
    Reactions[1][1][0]=1;		Reactions[1][1][1]=1;		Reactions[1][1][2]=1;

    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:          a = 1,   b = 0
    //   Impact ionization: a = 0,   b = 1
    //
    // Backward Reactions
    //   Recomb ionization:      a = 0, b = 1
    //   Impact ionization:      a = 0, b = 1
    //   N2 impact dissociation: a = 0, b = 1
    //   Others:                 a = 1, b = 0
    Tcf_a[0] = 0.5; Tcf_b[0] = 0.5; Tcb_a[0] = 1;  Tcb_b[0] = 0;
    Tcf_a[1] = 0.5; Tcf_b[1] = 0.5; Tcb_a[1] = 1;  Tcb_b[1] = 0;

    /*--- Dissociation potential [KJ/kg] ---*/
    Diss[0] = 3.36E4;
    Diss[1] = 0.0;

    /*--- Collision integral data ---*/
    Omega00[0][0][0] = -6.0614558E-03;  Omega00[0][0][1] = 1.2689102E-01;   Omega00[0][0][2] = -1.0616948E+00;  Omega00[0][0][3] = 8.0955466E+02;
    Omega00[0][1][0] = -1.0796249E-02;  Omega00[0][1][1] = 2.2656509E-01;   Omega00[0][1][2] = -1.7910602E+00;  Omega00[0][1][3] = 4.0455218E+03;
    Omega00[1][0][0] = -1.0796249E-02;  Omega00[1][0][1] = 2.2656509E-01;   Omega00[1][0][2] = -1.7910602E+00;  Omega00[1][0][3] = 4.0455218E+03;
    Omega00[1][1][0] = -9.6083779E-03;  Omega00[1][1][1] = 2.0938971E-01;   Omega00[1][1][2] = -1.7386904E+00;  Omega00[1][1][3] = 3.3587983E+03;

    Omega11[0][0][0] = -7.6303990E-03;  Omega11[0][0][1] = 1.6878089E-01;   Omega11[0][0][2] = -1.4004234E+00;  Omega11[0][0][3] = 2.1427708E+03;
    Omega11[0][1][0] = -8.3493693E-03;  Omega11[0][1][1] = 1.7808911E-01;   Omega11[0][1][2] = -1.4466155E+00;  Omega11[0][1][3] = 1.9324210E+03;
    Omega11[1][0][0] = -8.3493693E-03;  Omega11[1][0][1] = 1.7808911E-01;   Omega11[1][0][2] = -1.4466155E+00;  Omega11[1][0][3] = 1.9324210E+03;
    Omega11[1][1][0] = -7.7439615E-03;  Omega11[1][1][1] = 1.7129007E-01;   Omega11[1][1][2] = -1.4809088E+00;  Omega11[1][1][3] = 2.1284951E+03;

    break;

  case AIR5:

    /*--- Define parameters of the gas model ---*/
    nReactions  = 17;
    ionization  = false;

    /*--- Allocate vectors for gas properties ---*/
    Wall_Catalycity      = new su2double[nSpecies];
    Molar_Mass           = new su2double[nSpecies];
    CharVibTemp          = new su2double[nSpecies];
    RotationModes        = new su2double[nSpecies];
    Enthalpy_Formation   = new su2double[nSpecies];
    Ref_Temperature      = new su2double[nSpecies];
    ArrheniusCoefficient = new su2double[nReactions];
    ArrheniusEta         = new su2double[nReactions];
    ArrheniusTheta       = new su2double[nReactions];
    Tcf_a                = new su2double[nReactions];
    Tcf_b                = new su2double[nReactions];
    Tcb_a                = new su2double[nReactions];
    Tcb_b                = new su2double[nReactions];
    nElStates            = new unsigned short[nSpecies];
    Reactions            = new int**[nReactions];
    for (unsigned short iRxn = 0; iRxn < nReactions; iRxn++) {
      Reactions[iRxn] = new int*[2];
      for (unsigned short ii = 0; ii < 2; ii++)
        Reactions[iRxn][ii] = new int[6];
    }

    Blottner  = new su2double*[nSpecies];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Blottner[iSpecies] = new su2double[3];

    // Omega[iSpecies][jSpecies][iCoeff]
    Omega00 = new su2double**[nSpecies];
    Omega11 = new su2double**[nSpecies];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Omega00[iSpecies] = new su2double*[nSpecies];
      Omega11[iSpecies] = new su2double*[nSpecies];
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        Omega00[iSpecies][jSpecies] = new su2double[4];
        Omega11[iSpecies][jSpecies] = new su2double[4];
      }
    }

    // Wall mass fractions for catalytic boundaries
    Wall_Catalycity[0] = 0.4;
    Wall_Catalycity[1] = 0.4;
    Wall_Catalycity[2] = 0.1;
    Wall_Catalycity[3] = 0.05;
    Wall_Catalycity[4] = 0.05;

    // Free stream mass fractions
    MassFrac_FreeStream = new su2double[nSpecies];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      MassFrac_FreeStream[iSpecies] = Gas_Composition[iSpecies];

    /*--- Assign gas properties ---*/
    // Rotational modes of energy storage
    RotationModes[0] = 2.0;
    RotationModes[1] = 2.0;
    RotationModes[2] = 2.0;
    RotationModes[3] = 0.0;
    RotationModes[4] = 0.0;

    // Molar mass [kg/kmol]
    Molar_Mass[0] = 2.0*14.0067;
    Molar_Mass[1] = 2.0*15.9994;
    Molar_Mass[2] = 14.0067+15.9994;
    Molar_Mass[3] = 14.0067;
    Molar_Mass[4] = 15.9994;

    //Characteristic vibrational temperatures
    CharVibTemp[0] = 3395.0;
    CharVibTemp[1] = 2239.0;
    CharVibTemp[2] = 2817.0;
    CharVibTemp[3] = 0.0;
    CharVibTemp[4] = 0.0;

    // Formation enthalpy: (Scalabrin values, J/kg)
    Enthalpy_Formation[0] = 0.0;			//N2
    Enthalpy_Formation[1] = 0.0;			//O2
    Enthalpy_Formation[2] = 3.0E6;    //NO
    Enthalpy_Formation[3] = 3.36E7;		//N
    Enthalpy_Formation[4] = 1.54E7;		//O

    // Reference temperature (JANAF values, [K])
    Ref_Temperature[0] = 0.0;
    Ref_Temperature[1] = 0.0;
    Ref_Temperature[2] = 0.0;
    Ref_Temperature[3] = 0.0;
    Ref_Temperature[4] = 0.0;
    //        Ref_Temperature[2] = 298.15;
    //        Ref_Temperature[3] = 298.15;
    //        Ref_Temperature[4] = 298.15;

    // Blottner viscosity coefficients
    // A                        // B                        // C
    Blottner[0][0] = 2.68E-2;   Blottner[0][1] =  3.18E-1;  Blottner[0][2] = -1.13E1;  // N2
    Blottner[1][0] = 4.49E-2;   Blottner[1][1] = -8.26E-2;  Blottner[1][2] = -9.20E0;  // O2
    Blottner[2][0] = 4.36E-2;   Blottner[2][1] = -3.36E-2;  Blottner[2][2] = -9.58E0;  // NO
    Blottner[3][0] = 1.16E-2;   Blottner[3][1] =  6.03E-1;  Blottner[3][2] = -1.24E1;  // N
    Blottner[4][0] = 2.03E-2;   Blottner[4][1] =  4.29E-1;  Blottner[4][2] = -1.16E1;  // O

    // Number of electron states
    nElStates[0] = 15;                    // N2
    nElStates[1] = 7;                     // O2
    nElStates[2] = 16;                    // NO
    nElStates[3] = 3;                     // N
    nElStates[4] = 5;                     // O

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      maxEl = max(maxEl, nElStates[iSpecies]);

    /*--- Allocate electron data arrays ---*/
    CharElTemp = new su2double*[nSpecies];
    ElDegeneracy      = new su2double*[nSpecies];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      CharElTemp[iSpecies] = new su2double[maxEl];
      ElDegeneracy[iSpecies]      = new su2double[maxEl];
    }

    /*--- Initialize the arrays ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (iEl = 0; iEl < maxEl; iEl++) {
        CharElTemp[iSpecies][iEl] = 0.0;
        ElDegeneracy[iSpecies][iEl] = 0.0;
      }
    }

    //N2: 15 states
    CharElTemp[0][0]  = 0.000000000000000E+00;
    CharElTemp[0][1]  = 7.223156514095200E+04;
    CharElTemp[0][2]  = 8.577862640384000E+04;
    CharElTemp[0][3]  = 8.605026716160000E+04;
    CharElTemp[0][4]  = 9.535118627874400E+04;
    CharElTemp[0][5]  = 9.805635702203200E+04;
    CharElTemp[0][6]  = 9.968267656935200E+04;
    CharElTemp[0][7]  = 1.048976467715200E+05;
    CharElTemp[0][8]  = 1.116489555200000E+05;
    CharElTemp[0][9]  = 1.225836470400000E+05;
    CharElTemp[0][10] = 1.248856873600000E+05;
    CharElTemp[0][11] = 1.282476158188320E+05;
    CharElTemp[0][12] = 1.338060936000000E+05;
    CharElTemp[0][13] = 1.404296391107200E+05;
    CharElTemp[0][14] = 1.504958859200000E+05;
    ElDegeneracy[0][0]  = 1;
    ElDegeneracy[0][1]  = 3;
    ElDegeneracy[0][2]  = 6;
    ElDegeneracy[0][3]  = 6;
    ElDegeneracy[0][4]  = 3;
    ElDegeneracy[0][5]  = 1;
    ElDegeneracy[0][6]  = 2;
    ElDegeneracy[0][7]  = 2;
    ElDegeneracy[0][8]  = 5;
    ElDegeneracy[0][9]  = 1;
    ElDegeneracy[0][10] = 6;
    ElDegeneracy[0][11] = 6;
    ElDegeneracy[0][12] = 10;
    ElDegeneracy[0][13] = 6;
    ElDegeneracy[0][14] = 6;

    // O2: 7 states
    CharElTemp[1][0] = 0.000000000000000E+00;
    CharElTemp[1][1] = 1.139156019700800E+04;
    CharElTemp[1][2] = 1.898473947826400E+04;
    CharElTemp[1][3] = 4.755973576639200E+04;
    CharElTemp[1][4] = 4.991242097343200E+04;
    CharElTemp[1][5] = 5.092268575561600E+04;
    CharElTemp[1][6] = 7.189863255967200E+04;
    ElDegeneracy[1][0] = 3;
    ElDegeneracy[1][1] = 2;
    ElDegeneracy[1][2] = 1;
    ElDegeneracy[1][3] = 1;
    ElDegeneracy[1][4] = 6;
    ElDegeneracy[1][5] = 3;
    ElDegeneracy[1][6] = 3;

    // NO: 16 states
    CharElTemp[2][0]  = 0.000000000000000E+00;
    CharElTemp[2][1]  = 5.467345760000000E+04;
    CharElTemp[2][2]  = 6.317139627802400E+04;
    CharElTemp[2][3]  = 6.599450342445600E+04;
    CharElTemp[2][4]  = 6.906120960000000E+04;
    CharElTemp[2][5]  = 7.049998480000000E+04;
    CharElTemp[2][6]  = 7.491055017560000E+04;
    CharElTemp[2][7]  = 7.628875293968000E+04;
    CharElTemp[2][8]  = 8.676188537552000E+04;
    CharElTemp[2][9]  = 8.714431182368000E+04;
    CharElTemp[2][10] = 8.886077063728000E+04;
    CharElTemp[2][11] = 8.981755614528000E+04;
    CharElTemp[2][12] = 8.988445919208000E+04;
    CharElTemp[2][13] = 9.042702132000000E+04;
    CharElTemp[2][14] = 9.064283760000000E+04;
    CharElTemp[2][15] = 9.111763341600000E+04;
    ElDegeneracy[2][0]  = 4;
    ElDegeneracy[2][1]  = 8;
    ElDegeneracy[2][2]  = 2;
    ElDegeneracy[2][3]  = 4;
    ElDegeneracy[2][4]  = 4;
    ElDegeneracy[2][5]  = 4;
    ElDegeneracy[2][6]  = 4;
    ElDegeneracy[2][7]  = 2;
    ElDegeneracy[2][8]  = 4;
    ElDegeneracy[2][9]  = 2;
    ElDegeneracy[2][10] = 4;
    ElDegeneracy[2][11] = 4;
    ElDegeneracy[2][12] = 2;
    ElDegeneracy[2][13] = 2;
    ElDegeneracy[2][14] = 2;
    ElDegeneracy[2][15] = 4;

    // N: 3 states
    CharElTemp[3][0] = 0.000000000000000E+00;
    CharElTemp[3][1] = 2.766469645581980E+04;
    CharElTemp[3][2] = 4.149309313560210E+04;
    ElDegeneracy[3][0] = 4;
    ElDegeneracy[3][1] = 10;
    ElDegeneracy[3][2] = 6;

    // O: 5 states
    CharElTemp[4][0] = 0.000000000000000E+00;
    CharElTemp[4][1] = 2.277077570280000E+02;
    CharElTemp[4][2] = 3.265688785704000E+02;
    CharElTemp[4][3] = 2.283028632262240E+04;
    CharElTemp[4][4] = 4.861993036434160E+04;
    ElDegeneracy[4][0] = 5;
    ElDegeneracy[4][1] = 3;
    ElDegeneracy[4][2] = 1;
    ElDegeneracy[4][3] = 5;
    ElDegeneracy[4][4] = 1;

    /*--- Set reaction maps ---*/
    // N2 dissociation
    Reactions[0][0][0]=0;		Reactions[0][0][1]=0;		Reactions[0][0][2]=nSpecies;		Reactions[0][1][0]=3;		Reactions[0][1][1]=3;		Reactions[0][1][2] =0;
    Reactions[1][0][0]=0;		Reactions[1][0][1]=1;		Reactions[1][0][2]=nSpecies;		Reactions[1][1][0]=3;		Reactions[1][1][1]=3;		Reactions[1][1][2] =1;
    Reactions[2][0][0]=0;		Reactions[2][0][1]=2;		Reactions[2][0][2]=nSpecies;		Reactions[2][1][0]=3;		Reactions[2][1][1]=3;		Reactions[2][1][2] =2;
    Reactions[3][0][0]=0;		Reactions[3][0][1]=3;		Reactions[3][0][2]=nSpecies;		Reactions[3][1][0]=3;		Reactions[3][1][1]=3;		Reactions[3][1][2] =3;
    Reactions[4][0][0]=0;		Reactions[4][0][1]=4;		Reactions[4][0][2]=nSpecies;		Reactions[4][1][0]=3;		Reactions[4][1][1]=3;		Reactions[4][1][2] =4;
    // O2 dissociation
    Reactions[5][0][0]=1;		Reactions[5][0][1]=0;		Reactions[5][0][2]=nSpecies;		Reactions[5][1][0]=4;		Reactions[5][1][1]=4;		Reactions[5][1][2] =0;
    Reactions[6][0][0]=1;		Reactions[6][0][1]=1;		Reactions[6][0][2]=nSpecies;		Reactions[6][1][0]=4;		Reactions[6][1][1]=4;		Reactions[6][1][2] =1;
    Reactions[7][0][0]=1;		Reactions[7][0][1]=2;		Reactions[7][0][2]=nSpecies;		Reactions[7][1][0]=4;		Reactions[7][1][1]=4;		Reactions[7][1][2] =2;
    Reactions[8][0][0]=1;		Reactions[8][0][1]=3;		Reactions[8][0][2]=nSpecies;		Reactions[8][1][0]=4;		Reactions[8][1][1]=4;		Reactions[8][1][2] =3;
    Reactions[9][0][0]=1;		Reactions[9][0][1]=4;		Reactions[9][0][2]=nSpecies;		Reactions[9][1][0]=4;		Reactions[9][1][1]=4;		Reactions[9][1][2] =4;
    // NO dissociation
    Reactions[10][0][0]=2;		Reactions[10][0][1]=0;		Reactions[10][0][2]=nSpecies;		Reactions[10][1][0]=3;		Reactions[10][1][1]=4;		Reactions[10][1][2] =0;
    Reactions[11][0][0]=2;		Reactions[11][0][1]=1;		Reactions[11][0][2]=nSpecies;		Reactions[11][1][0]=3;		Reactions[11][1][1]=4;		Reactions[11][1][2] =1;
    Reactions[12][0][0]=2;		Reactions[12][0][1]=2;		Reactions[12][0][2]=nSpecies;		Reactions[12][1][0]=3;		Reactions[12][1][1]=4;		Reactions[12][1][2] =2;
    Reactions[13][0][0]=2;		Reactions[13][0][1]=3;		Reactions[13][0][2]=nSpecies;		Reactions[13][1][0]=3;		Reactions[13][1][1]=4;		Reactions[13][1][2] =3;
    Reactions[14][0][0]=2;		Reactions[14][0][1]=4;		Reactions[14][0][2]=nSpecies;		Reactions[14][1][0]=3;		Reactions[14][1][1]=4;		Reactions[14][1][2] =4;
    // N2 + O -> NO + N
    Reactions[15][0][0]=0;		Reactions[15][0][1]=4;		Reactions[15][0][2]=nSpecies;		Reactions[15][1][0]=2;		Reactions[15][1][1]=3;		Reactions[15][1][2]= nSpecies;
    // NO + O -> O2 + N
    Reactions[16][0][0]=2;		Reactions[16][0][1]=4;		Reactions[16][0][2]=nSpecies;		Reactions[16][1][0]=1;		Reactions[16][1][1]=3;		Reactions[16][1][2]= nSpecies;

    /*--- Set Arrhenius coefficients for reactions ---*/
    // Pre-exponential factor
    ArrheniusCoefficient[0]  = 7.0E21;
    ArrheniusCoefficient[1]  = 7.0E21;
    ArrheniusCoefficient[2]  = 7.0E21;
    ArrheniusCoefficient[3]  = 3.0E22;
    ArrheniusCoefficient[4]  = 3.0E22;
    ArrheniusCoefficient[5]  = 2.0E21;
    ArrheniusCoefficient[6]  = 2.0E21;
    ArrheniusCoefficient[7]  = 2.0E21;
    ArrheniusCoefficient[8]  = 1.0E22;
    ArrheniusCoefficient[9]  = 1.0E22;
    ArrheniusCoefficient[10] = 5.0E15;
    ArrheniusCoefficient[11] = 5.0E15;
    ArrheniusCoefficient[12] = 5.0E15;
    ArrheniusCoefficient[13] = 1.1E17;
    ArrheniusCoefficient[14] = 1.1E17;
    ArrheniusCoefficient[15] = 6.4E17;
    ArrheniusCoefficient[16] = 8.4E12;

    // Rate-controlling temperature exponent
    ArrheniusEta[0]  = -1.60;
    ArrheniusEta[1]  = -1.60;
    ArrheniusEta[2]  = -1.60;
    ArrheniusEta[3]  = -1.60;
    ArrheniusEta[4]  = -1.60;
    ArrheniusEta[5]  = -1.50;
    ArrheniusEta[6]  = -1.50;
    ArrheniusEta[7]  = -1.50;
    ArrheniusEta[8]  = -1.50;
    ArrheniusEta[9]  = -1.50;
    ArrheniusEta[10] = 0.0;
    ArrheniusEta[11] = 0.0;
    ArrheniusEta[12] = 0.0;
    ArrheniusEta[13] = 0.0;
    ArrheniusEta[14] = 0.0;
    ArrheniusEta[15] = -1.0;
    ArrheniusEta[16] = 0.0;

    // Characteristic temperature
    ArrheniusTheta[0]  = 113200.0;
    ArrheniusTheta[1]  = 113200.0;
    ArrheniusTheta[2]  = 113200.0;
    ArrheniusTheta[3]  = 113200.0;
    ArrheniusTheta[4]  = 113200.0;
    ArrheniusTheta[5]  = 59500.0;
    ArrheniusTheta[6]  = 59500.0;
    ArrheniusTheta[7]  = 59500.0;
    ArrheniusTheta[8]  = 59500.0;
    ArrheniusTheta[9]  = 59500.0;
    ArrheniusTheta[10] = 75500.0;
    ArrheniusTheta[11] = 75500.0;
    ArrheniusTheta[12] = 75500.0;
    ArrheniusTheta[13] = 75500.0;
    ArrheniusTheta[14] = 75500.0;
    ArrheniusTheta[15] = 38400.0;
    ArrheniusTheta[16] = 19450.0;

    /*--- Set rate-controlling temperature exponents ---*/
    //  -----------  Tc = Ttr^a * Tve^b  -----------
    //
    // Forward Reactions
    //   Dissociation:      a = 0.5, b = 0.5  (OR a = 0.7, b =0.3)
    //   Exchange:          a = 1,   b = 0
    //   Impact ionization: a = 0,   b = 1
    //
    // Backward Reactions
    //   Recomb ionization:      a = 0, b = 1
    //   Impact ionization:      a = 0, b = 1
    //   N2 impact dissociation: a = 0, b = 1
    //   Others:                 a = 1, b = 0
    Tcf_a[0]  = 0.5; Tcf_b[0]  = 0.5; Tcb_a[0]  = 1;  Tcb_b[0] = 0;
    Tcf_a[1]  = 0.5; Tcf_b[1]  = 0.5; Tcb_a[1]  = 1;  Tcb_b[1] = 0;
    Tcf_a[2]  = 0.5; Tcf_b[2]  = 0.5; Tcb_a[2]  = 1;  Tcb_b[2] = 0;
    Tcf_a[3]  = 0.5; Tcf_b[3]  = 0.5; Tcb_a[3]  = 1;  Tcb_b[3] = 0;
    Tcf_a[4]  = 0.5; Tcf_b[4]  = 0.5; Tcb_a[4]  = 1;  Tcb_b[4] = 0;

    Tcf_a[5]  = 0.5; Tcf_b[5]  = 0.5; Tcb_a[5]  = 1;  Tcb_b[5] = 0;
    Tcf_a[6]  = 0.5; Tcf_b[6]  = 0.5; Tcb_a[6]  = 1;  Tcb_b[6] = 0;
    Tcf_a[7]  = 0.5; Tcf_b[7]  = 0.5; Tcb_a[7]  = 1;  Tcb_b[7] = 0;
    Tcf_a[8]  = 0.5; Tcf_b[8]  = 0.5; Tcb_a[8]  = 1;  Tcb_b[8] = 0;
    Tcf_a[9]  = 0.5; Tcf_b[9]  = 0.5; Tcb_a[9]  = 1;  Tcb_b[9] = 0;

    Tcf_a[10] = 0.5; Tcf_b[10] = 0.5; Tcb_a[10] = 1;  Tcb_b[10] = 0;
    Tcf_a[11] = 0.5; Tcf_b[11] = 0.5; Tcb_a[11] = 1;  Tcb_b[11] = 0;
    Tcf_a[12] = 0.5; Tcf_b[12] = 0.5; Tcb_a[12] = 1;  Tcb_b[12] = 0;
    Tcf_a[13] = 0.5; Tcf_b[13] = 0.5; Tcb_a[13] = 1;  Tcb_b[13] = 0;
    Tcf_a[14] = 0.5; Tcf_b[14] = 0.5; Tcb_a[14] = 1;  Tcb_b[14] = 0;

    Tcf_a[15] = 1.0; Tcf_b[15] = 0.0; Tcb_a[15] = 1;  Tcb_b[15] = 0;
    Tcf_a[16] = 1.0; Tcf_b[16] = 0.0; Tcb_a[16] = 1;  Tcb_b[16] = 0;

    /*--- Collision integral data ---*/
    // Omega(0,0) ----------------------
    //N2
    Omega00[0][0][0] = -6.0614558E-03;  Omega00[0][0][1] = 1.2689102E-01;   Omega00[0][0][2] = -1.0616948E+00;  Omega00[0][0][3] = 8.0955466E+02;
    Omega00[0][1][0] = -3.7959091E-03;  Omega00[0][1][1] = 9.5708295E-02;   Omega00[0][1][2] = -1.0070611E+00;  Omega00[0][1][3] = 8.9392313E+02;
    Omega00[0][2][0] = -1.9295666E-03;  Omega00[0][2][1] = 2.7995735E-02;   Omega00[0][2][2] = -3.1588514E-01;  Omega00[0][2][3] = 1.2880734E+02;
    Omega00[0][3][0] = -1.0796249E-02;  Omega00[0][3][1] = 2.2656509E-01;   Omega00[0][3][2] = -1.7910602E+00;  Omega00[0][3][3] = 4.0455218E+03;
    Omega00[0][4][0] = -2.7244269E-03;  Omega00[0][4][1] = 6.9587171E-02;   Omega00[0][4][2] = -7.9538667E-01;  Omega00[0][4][3] = 4.0673730E+02;
    //O2
    Omega00[1][0][0] = -3.7959091E-03;  Omega00[1][0][1] = 9.5708295E-02;   Omega00[1][0][2] = -1.0070611E+00;  Omega00[1][0][3] = 8.9392313E+02;
    Omega00[1][1][0] = -8.0682650E-04;  Omega00[1][1][1] = 1.6602480E-02;   Omega00[1][1][2] = -3.1472774E-01;  Omega00[1][1][3] = 1.4116458E+02;
    Omega00[1][2][0] = -6.4433840E-04;  Omega00[1][2][1] = 8.5378580E-03;   Omega00[1][2][2] = -2.3225102E-01;  Omega00[1][2][3] = 1.1371608E+02;
    Omega00[1][3][0] = -1.1453028E-03;  Omega00[1][3][1] = 1.2654140E-02;   Omega00[1][3][2] = -2.2435218E-01;  Omega00[1][3][3] = 7.7201588E+01;
    Omega00[1][4][0] = -4.8405803E-03;  Omega00[1][4][1] = 1.0297688E-01;   Omega00[1][4][2] = -9.6876576E-01;  Omega00[1][4][3] = 6.1629812E+02;
    //NO
    Omega00[2][0][0] = -1.9295666E-03;  Omega00[2][0][1] = 2.7995735E-02;   Omega00[2][0][2] = -3.1588514E-01;  Omega00[2][0][3] = 1.2880734E+02;
    Omega00[2][1][0] = -6.4433840E-04;  Omega00[2][1][1] = 8.5378580E-03;   Omega00[2][1][2] = -2.3225102E-01;  Omega00[2][1][3] = 1.1371608E+02;
    Omega00[2][2][0] = -0.0000000E+00;  Omega00[2][2][1] = -1.1056066E-02;  Omega00[2][2][2] = -5.9216250E-02;  Omega00[2][2][3] = 7.2542367E+01;
    Omega00[2][3][0] = -1.5770918E-03;  Omega00[2][3][1] = 1.9578381E-02;   Omega00[2][3][2] = -2.7873624E-01;  Omega00[2][3][3] = 9.9547944E+01;
    Omega00[2][4][0] = -1.0885815E-03;  Omega00[2][4][1] = 1.1883688E-02;   Omega00[2][4][2] = -2.1844909E-01;  Omega00[2][4][3] = 7.5512560E+01;
    //N
    Omega00[3][0][0] = -1.0796249E-02;  Omega00[3][0][1] = 2.2656509E-01;   Omega00[3][0][2] = -1.7910602E+00;  Omega00[3][0][3] = 4.0455218E+03;
    Omega00[3][1][0] = -1.1453028E-03;  Omega00[3][1][1] = 1.2654140E-02;   Omega00[3][1][2] = -2.2435218E-01;  Omega00[3][1][3] = 7.7201588E+01;
    Omega00[3][2][0] = -1.5770918E-03;  Omega00[3][2][1] = 1.9578381E-02;   Omega00[3][2][2] = -2.7873624E-01;  Omega00[3][2][3] = 9.9547944E+01;
    Omega00[3][3][0] = -9.6083779E-03;  Omega00[3][3][1] = 2.0938971E-01;   Omega00[3][3][2] = -1.7386904E+00;  Omega00[3][3][3] = 3.3587983E+03;
    Omega00[3][4][0] = -7.8147689E-03;  Omega00[3][4][1] = 1.6792705E-01;   Omega00[3][4][2] = -1.4308628E+00;  Omega00[3][4][3] = 1.6628859E+03;
    //O
    Omega00[4][0][0] = -2.7244269E-03;  Omega00[4][0][1] = 6.9587171E-02;   Omega00[4][0][2] = -7.9538667E-01;  Omega00[4][0][3] = 4.0673730E+02;
    Omega00[4][1][0] = -4.8405803E-03;  Omega00[4][1][1] = 1.0297688E-01;   Omega00[4][1][2] = -9.6876576E-01;  Omega00[4][1][3] = 6.1629812E+02;
    Omega00[4][2][0] = -1.0885815E-03;  Omega00[4][2][1] = 1.1883688E-02;   Omega00[4][2][2] = -2.1844909E-01;  Omega00[4][2][3] = 7.5512560E+01;
    Omega00[4][3][0] = -7.8147689E-03;  Omega00[4][3][1] = 1.6792705E-01;   Omega00[4][3][2] = -1.4308628E+00;  Omega00[4][3][3] = 1.6628859E+03;
    Omega00[4][4][0] = -6.4040535E-03;  Omega00[4][4][1] = 1.4629949E-01;   Omega00[4][4][2] = -1.3892121E+00;  Omega00[4][4][3] = 2.0903441E+03;

    // Omega(1,1) ----------------------
    //N2
    Omega11[0][0][0] = -7.6303990E-03;  Omega11[0][0][1] = 1.6878089E-01;   Omega11[0][0][2] = -1.4004234E+00;  Omega11[0][0][3] = 2.1427708E+03;
    Omega11[0][1][0] = -8.0457321E-03;  Omega11[0][1][1] = 1.9228905E-01;   Omega11[0][1][2] = -1.7102854E+00;  Omega11[0][1][3] = 5.2213857E+03;
    Omega11[0][2][0] = -6.8237776E-03;  Omega11[0][2][1] = 1.4360616E-01;   Omega11[0][2][2] = -1.1922240E+00;  Omega11[0][2][3] = 1.2433086E+03;
    Omega11[0][3][0] = -8.3493693E-03;  Omega11[0][3][1] = 1.7808911E-01;   Omega11[0][3][2] = -1.4466155E+00;  Omega11[0][3][3] = 1.9324210E+03;
    Omega11[0][4][0] = -8.3110691E-03;  Omega11[0][4][1] = 1.9617877E-01;   Omega11[0][4][2] = -1.7205427E+00;  Omega11[0][4][3] = 4.0812829E+03;
    //O2
    Omega11[1][0][0] = -8.0457321E-03;  Omega11[1][0][1] = 1.9228905E-01;   Omega11[1][0][2] = -1.7102854E+00;  Omega11[1][0][3] = 5.2213857E+03;
    Omega11[1][1][0] = -6.2931612E-03;  Omega11[1][1][1] = 1.4624645E-01;   Omega11[1][1][2] = -1.3006927E+00;  Omega11[1][1][3] = 1.8066892E+03;
    Omega11[1][2][0] = -6.8508672E-03;  Omega11[1][2][1] = 1.5524564E-01;   Omega11[1][2][2] = -1.3479583E+00;  Omega11[1][2][3] = 2.0037890E+03;
    Omega11[1][3][0] = -1.0608832E-03;  Omega11[1][3][1] = 1.1782595E-02;   Omega11[1][3][2] = -2.1246301E-01;  Omega11[1][3][3] = 8.4561598E+01;
    Omega11[1][4][0] = -3.7969686E-03;  Omega11[1][4][1] = 7.6789981E-02;   Omega11[1][4][2] = -7.3056809E-01;  Omega11[1][4][3] = 3.3958171E+02;
    //NO
    Omega11[2][0][0] = -6.8237776E-03;  Omega11[2][0][1] = 1.4360616E-01;   Omega11[2][0][2] = -1.1922240E+00;  Omega11[2][0][3] = 1.2433086E+03;
    Omega11[2][1][0] = -6.8508672E-03;  Omega11[2][1][1] = 1.5524564E-01;   Omega11[2][1][2] = -1.3479583E+00;  Omega11[2][1][3] = 2.0037890E+03;
    Omega11[2][2][0] = -7.4942466E-03;  Omega11[2][2][1] = 1.6626193E-01;   Omega11[2][2][2] = -1.4107027E+00;  Omega11[2][2][3] = 2.3097604E+03;
    Omega11[2][3][0] = -1.4719259E-03;  Omega11[2][3][1] = 1.8446968E-02;   Omega11[2][3][2] = -2.6460411E-01;  Omega11[2][3][3] = 1.0911124E+02;
    Omega11[2][4][0] = -1.0066279E-03;  Omega11[2][4][1] = 1.1029264E-02;   Omega11[2][4][2] = -2.0671266E-01;  Omega11[2][4][3] = 8.2644384E+01;
    //N
    Omega11[3][0][0] = -8.3493693E-03;  Omega11[3][0][1] = 1.7808911E-01;   Omega11[3][0][2] = -1.4466155E+00;  Omega11[3][0][3] = 1.9324210E+03;
    Omega11[3][1][0] = -1.0608832E-03;  Omega11[3][1][1] = 1.1782595E-02;   Omega11[3][1][2] = -2.1246301E-01;  Omega11[3][1][3] = 8.4561598E+01;
    Omega11[3][2][0] = -1.4719259E-03;  Omega11[3][2][1] = 1.8446968E-02;   Omega11[3][2][2] = -2.6460411E-01;  Omega11[3][2][3] = 1.0911124E+02;
    Omega11[3][3][0] = -7.7439615E-03;  Omega11[3][3][1] = 1.7129007E-01;   Omega11[3][3][2] = -1.4809088E+00;  Omega11[3][3][3] = 2.1284951E+03;
    Omega11[3][4][0] = -5.0478143E-03;  Omega11[3][4][1] = 1.0236186E-01;   Omega11[3][4][2] = -9.0058935E-01;  Omega11[3][4][3] = 4.4472565E+02;
    //O
    Omega11[4][0][0] = -8.3110691E-03;  Omega11[4][0][1] = 1.9617877E-01;   Omega11[4][0][2] = -1.7205427E+00;  Omega11[4][0][3] = 4.0812829E+03;
    Omega11[4][1][0] = -3.7969686E-03;  Omega11[4][1][1] = 7.6789981E-02;   Omega11[4][1][2] = -7.3056809E-01;  Omega11[4][1][3] = 3.3958171E+02;
    Omega11[4][2][0] = -1.0066279E-03;  Omega11[4][2][1] = 1.1029264E-02;   Omega11[4][2][2] = -2.0671266E-01;  Omega11[4][2][3] = 8.2644384E+01;
    Omega11[4][3][0] = -5.0478143E-03;  Omega11[4][3][1] = 1.0236186E-01;   Omega11[4][3][2] = -9.0058935E-01;  Omega11[4][3][3] = 4.4472565E+02;
    Omega11[4][4][0] = -4.2451096E-03;  Omega11[4][4][1] = 9.6820337E-02;   Omega11[4][4][2] = -9.9770795E-01;  Omega11[4][4][3] = 8.3320644E+02;

    break;
  }
}

su2double CTNE2Gas::Calc_CvTraRotSpecies(su2double *Ms, su2double Ru, unsigned short val_Species) {

  su2double Cv_tr_ks;

  Cv_tr_ks = (3.0/2.0 + RotationModes[val_Species]/2.0) * Ru/Ms[val_Species];

  return Cv_tr_ks;

}

su2double CTNE2Gas::Calc_CvVibElSpecies(su2double val_Tve, unsigned short val_Species) {

  unsigned short iEl;
  su2double *Ms, *thetav, **thetae, **g, RuSI, Ru;
  su2double thoTve, exptv, thsqr, Cvvs, Cves;
  su2double Tve;
  su2double num, num2, num3, denom;

  /*--- Read from config ---*/
  thetav    = CharVibTemp;
  thetae    = CharElTemp;
  g         = ElDegeneracy;
  Ms        = Molar_Mass;

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

su2double* CTNE2Gas::Calc_CpTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {}

su2double* CTNE2Gas::Calc_CpVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve) {}

su2double CTNE2Gas::Calc_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve){}

su2double CTNE2Gas::Calc_MixtureEnergy(su2double *cs, su2double sqvel, su2double rho,
                                       su2double T, su2double Tve) {


  unsigned short iEl;
  su2double Ef, rhos, Ev, Ee;

  su2double rhoE = 0.0;
  su2double num  = 0.0, denom = 0.0;

  su2double RuSI = UNIVERSAL_GAS_CONSTANT;
  su2double Ru = 1000*RuSI;

  /*--- Rename for convenience ---*/
  su2double *Ms      = Molar_Mass;
  su2double **g      = ElDegeneracy;
  su2double *hf      = Enthalpy_Formation;
  su2double *Tref    = Ref_Temperature;
  su2double *thetav  = CharVibTemp;
  su2double **thetae = CharElTemp;
  su2double *xi      = RotationModes;


  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    // Species density
    rhos = cs[iSpecies]*rho;

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
  }

  return rhoE;
}

su2double CTNE2Gas::Calc_MixtureEnergies(su2double* cs, su2double rho, su2double T, su2double Tve) {}

su2double CTNE2Gas::Calc_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve) {}

su2double CTNE2Gas::Calc_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve) {}

su2double CTNE2Gas::Calc_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve) {}

su2double CTNE2Gas::Calc_Enthalpies(su2double val_T, su2double val_eves, unsigned short val_Species) {

  su2double RuSI, Ru, *xi, *Ms, *hf, *Tref, T, eve, ef, hs;

  /*--- Read from config ---*/
  xi   = RotationModes;
  Ms   = Molar_Mass;
  hf   = Enthalpy_Formation;
  Tref = Ref_Temperature;

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

su2double CTNE2Gas::Calc_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve){}

su2double CTNE2Gas::Calc_SoundSpeed(su2double *cs, su2double rhoCvtr, su2double rho, su2double P){

  su2double soundspeed, conc = 0.0;
  su2double Ru  = 1000.0*UNIVERSAL_GAS_CONSTANT;
  su2double *Ms = Molar_Mass;

  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    conc += cs[iSpecies]*rho/Ms[iSpecies];
  }
  soundspeed = sqrt((1.0 + Ru/rhoCvtr*conc) * P/rho);

  return soundspeed;
}

su2double CTNE2Gas::Calc_Density(su2double *MassFrac, su2double T, su2double Tve, su2double P){

  su2double denom = 0.0;
  su2double RuSi = UNIVERSAL_GAS_CONSTANT;
  su2double Ru   = 1000*RuSi;

  /*--- Compute Density using mixture quantities ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++)
    denom += MassFrac[iSpecies] * (Ru/Molar_Mass[iSpecies]) * T;
  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    denom += MassFrac[nSpecies-1] * (Ru/Molar_Mass[nSpecies-1]) * Tve;
  Density = P / denom;

  return Density;
}

void CTNE2Gas::CalcdPdU(su2double *V, su2double *val_eves, su2double *val_dPdU) {

  // Note: Requires SetDensity(), SetTemperature(), SetPressure(), & SetGasProperties()
  // Note: Electron energy not included properly.

  unsigned short iDim, iSpecies, iEl;
  su2double *Ms, *Tref, *hf, *xi, *thetav, **thetae, **g;
  su2double RuSI, Ru, RuBAR, CvtrBAR, rhoCvtr, rhoCvve, Cvtrs, rho_el, sqvel, conc;
  su2double rho, rhos, T, Tve, ef;
  //su2double num, denom;

  /*--- Determine the number of heavy species ---*/
  if (ionization) { rho_el = V[RHOS_INDEX+nSpecies-1];
  } else          { rho_el = 0.0; }

  /*--- Read gas mixture properties from config ---*/
  Ms        = Molar_Mass;
  Tref      = Ref_Temperature;
  hf        = Enthalpy_Formation;
  xi        = RotationModes;
  thetav    = CharVibTemp;
  thetae    = CharElTemp;
  g         = ElDegeneracy;

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

su2double CTNE2Gas::CalcEve(su2double val_Tve, unsigned short val_Species) {

  unsigned short iEl;
  su2double *Ms, *thetav, **thetae, **g, *hf, *Tref, RuSI, Ru;
  su2double Tve, Ev, Eel, Ef;
  su2double num, denom;

  /*--- Rename for convenience ---*/
  RuSI   = UNIVERSAL_GAS_CONSTANT;
  Ru     = 1000.0*RuSI;
  Tve    = val_Tve;
  Ms     = Molar_Mass;
  Tref   = Ref_Temperature;
  hf     = Enthalpy_Formation;
  thetav = CharVibTemp;
  thetae = CharElTemp;
  g      = ElDegeneracy;

  /*--- Electron species energy ---*/
  if ((ionization) && (val_Species == nSpecies-1)) {

    /*--- Calculate formation energy ---*/
    Ef = hf[val_Species] - Ru/Ms[val_Species] * Tref[val_Species];

    /*--- Electron t-r mode contributes to mixture vib-el energy ---*/
    Eel = (3.0/2.0) * Ru/Ms[val_Species] * (Tve - Tref[val_Species]) + Ef;
    Ev  = 0.0;

  }

  /*--- Heavy particle energy ---*/
  else {

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

void CTNE2Gas::CalcdTdU(su2double *V, su2double *val_dTdU) {

  unsigned short iDim;
  su2double *Ms, *xi, *Tref, *hf;
  su2double v2, ef, T, Cvtrs, rhoCvtr, RuSI, Ru;

  /*--- Get gas properties from config settings ---*/
  Ms   = Molar_Mass;
  xi   = RotationModes;
  hf   = Enthalpy_Formation;
  Tref = Ref_Temperature;

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

void CTNE2Gas::CalcdTvedU(su2double *V, su2double *val_eves, su2double *val_dTvedU) {

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

void CTNE2Gas::Calc_TransportCoeff(CConfig *config, su2double *V){

  switch (config->GetKind_TransCoeffModel()) {
  case WILKE:{

    su2double Xs[nSpecies], conc;
    su2double *Ms;

    /*--- Rename for convenience ---*/
    Ms   = Molar_Mass;

    /*--- Calculate species mole fraction ---*/
    conc = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Xs[iSpecies] = V[RHOS_INDEX+iSpecies]/Ms[iSpecies];
      conc        += Xs[iSpecies];
    }
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      Xs[iSpecies] = Xs[iSpecies]/conc;

    /*--- Compute transport properties ---*/
    Calc_DiffusionCoeff_WBE(V, Xs);
    Calc_Viscosity_WBE(V, Xs);
    Calc_ThermalConductivity_WBE(V, Xs);

    break;
  }
  case GUPTAYOS:{

    /*--- Compute transport properties ---*/
    Calc_DiffusionCoeff_GY(V);
    Calc_Viscosity_GY(V);
    Calc_ThermalConductivity_GY(V);

    break;
  }
  }
}

void CTNE2Gas::Calc_DiffusionCoeff_WBE(su2double *V, su2double *Xs) {

  su2double *Ms, Mi, Mj, M;
  su2double rho, T;
  su2double Omega_ij;
  su2double denom;

  /*--- Rename for Convenience ---*/
  Ms  = Molar_Mass;
  rho = V[RHO_INDEX];
  T   = V[T_INDEX];

  /*--- Calculate mixture molar mass (kg/mol) ---*/
  // Note: Species molar masses stored as kg/kmol, need 1E-3 conversion
  M = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    M += Ms[iSpecies]*Xs[iSpecies];
  M = M*1E-3;

  /*--- Solve for binary diffusion coefficients ---*/
  // Note: Dij = Dji, so only loop through req'd indices
  // Note: Correlation requires kg/mol, hence 1E-3 conversion from kg/kmol
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Mi = Ms[iSpecies]*1E-3;
    for (jSpecies = iSpecies; jSpecies < nSpecies; jSpecies++) {
      Mj = Ms[jSpecies]*1E-3;

      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
      Omega_ij = 1E-20/PI_NUMBER * Omega00[iSpecies][jSpecies][3]
                                 * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
                                 + Omega00[iSpecies][jSpecies][1]*log(T)
                                 + Omega00[iSpecies][jSpecies][2]);

      Dij[iSpecies][jSpecies] = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(rho*Omega_ij);
      Dij[jSpecies][iSpecies] = 7.1613E-25*M*sqrt(T*(1/Mi+1/Mj))/(rho*Omega_ij);
    }
  }

  /*--- Calculate species-mixture diffusion coefficient --*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    denom = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      if (jSpecies != iSpecies) {
        denom += Xs[jSpecies]/Dij[iSpecies][jSpecies];
      }
    }
    DiffusionCoeff[iSpecies] = (1-Xs[iSpecies])/denom;
  }

  /*--- Fill in primitive variable vector ---*/
  // DELETE ME, MAYBE UNNECESSARY...MORE CLEAR WHAT HAPPENING, THOUGHTs?
  for (iSpecies=0; iSpecies<nSpecies; iSpecies++)
    V[DIFF_COEFF_INDEX+iSpecies] = DiffusionCoeff[iSpecies];

}

void CTNE2Gas::Calc_Viscosity_WBE(su2double *V, su2double *Xs) {

  su2double phis[nSpecies], mus[nSpecies];
  su2double tmp1, tmp2;
  su2double *Ms, T;

  /*--- Rename for convenience ---*/
  Ms = Molar_Mass;
  T  = V[T_INDEX];

  /*--- Use Blottner's curve fits for species viscosity ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    mus[iSpecies] = 0.1*exp((Blottner[iSpecies][0]*log(T)  +
                             Blottner[iSpecies][1])*log(T) +
                             Blottner[iSpecies][2]);

  /*--- Determine species 'phi' value for Blottner model ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    phis[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      tmp1 = 1.0 + sqrt(mus[iSpecies]/mus[jSpecies])*pow(Ms[jSpecies]/Ms[iSpecies], 0.25);
      tmp2 = sqrt(8.0*(1.0+Ms[iSpecies]/Ms[jSpecies]));
      phis[iSpecies] += Xs[jSpecies]*tmp1*tmp1/tmp2;
    }
  }

  /*--- Calculate mixture laminar viscosity ---*/
  LamVisc = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    LamVisc += Xs[iSpecies]*mus[iSpecies]/phis[iSpecies];

  /*---Fill in primitive vector ---*/
  V[LAM_VISC_INDEX] = LamVisc;
}

void CTNE2Gas::Calc_ThermalConductivity_WBE( su2double *V, su2double *Xs) {

  unsigned short iSpecies, jSpecies;
  su2double *Ms, RuSI, Ru, *xi, Tve, T;
  su2double Cves;
  su2double phis[nSpecies], mus[nSpecies], ks[nSpecies], kves[nSpecies];
  su2double tmp1, tmp2;

  /*--- Extract constants ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;

  /*--- Rename for convenence ---*/
  Ms  = Molar_Mass;
  xi  = RotationModes;
  T   = V[T_INDEX];
  Tve = V[TVE_INDEX];

  /*--- Use Blottner's curve fits for species viscosity ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    mus[iSpecies] = 0.1*exp((Blottner[iSpecies][0]*log(T)  +
                             Blottner[iSpecies][1])*log(T) +
                             Blottner[iSpecies][2]);

  /*--- Determine species 'phi' value for Blottner model ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    phis[iSpecies] = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      tmp1 = 1.0 + sqrt(mus[iSpecies]/mus[jSpecies])*pow(Ms[jSpecies]/Ms[iSpecies], 0.25);
      tmp2 = sqrt(8.0*(1.0+Ms[iSpecies]/Ms[jSpecies]));
      phis[iSpecies] += Xs[jSpecies]*tmp1*tmp1/tmp2;
    }
  }

  /*--- Determine species tr & ve conductivities ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Cves = Calc_CvVibElSpecies(Tve,iSpecies);
    ks[iSpecies] = mus[iSpecies]*(15.0/4.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
    kves[iSpecies] = mus[iSpecies]*Cves;
  }

  /*--- Calculate mixture tr & ve conductivities ---*/
  ThermalCond    = 0.0;
  ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    ThermalCond    += Xs[iSpecies]*ks[iSpecies]/phis[iSpecies];
    ThermalCond_ve += Xs[iSpecies]*kves[iSpecies]/phis[iSpecies];
  }

  /*--- Fill in primitive vector ---*/
  V[K_INDEX]   = ThermalCond;
  V[KVE_INDEX] = ThermalCond_ve;
}

void CTNE2Gas::Calc_DiffusionCoeff_GY(su2double *V) {

  unsigned short iSpecies, jSpecies;
  su2double rho, T, Tve, P;
  su2double *Ms, Mi, Mj, pi, RuSI, Ru, kb, gam_i, gam_j, gam_t, Theta_v;
  su2double denom, d1_ij, D_ij;
  su2double Omega_ij;

  /*--- Extract constants ---*/
  pi   = PI_NUMBER;
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  kb   = BOLTZMANN_CONSTANT;

  /*--- Rename for convenience ---*/
  Ms   = Molar_Mass;
  rho  = V[RHO_INDEX];
  T    = V[T_INDEX];
  Tve  = V[TVE_INDEX];
  P    = V[P_INDEX];

  /*--- Calculate mixture gas constant ---*/
  gam_t = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    gam_t += V[RHOS_INDEX+iSpecies] / (rho*Ms[iSpecies]);
  }

  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {

    /*--- Initialize the species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = 0.0;

    /*--- Calculate molar concentration ---*/
    Mi      = Ms[iSpecies];
    gam_i   = V[RHOS_INDEX+iSpecies] / (rho*Mi);
    Theta_v = GetCharVibTemp(iSpecies);

    denom = 0.0;
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      if (jSpecies != iSpecies) {
        Mj    = Ms[jSpecies];
        gam_j = V[RHOS_INDEX+iSpecies] / (rho*Mj);

        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
            * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
            + Omega00[iSpecies][jSpecies][1]*log(T)
            + Omega00[iSpecies][jSpecies][2]);

        /*--- Calculate "delta1_ij" ---*/
        d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;

        /*--- Calculate heavy-particle binary diffusion coefficient ---*/
        D_ij = kb*T/(P*d1_ij);
        denom += gam_j/D_ij;
      }
    }

    if (ionization) {
      jSpecies = nSpecies-1;
      Mj       = GetMolar_Mass(jSpecies);
      gam_j    = V[RHOS_INDEX+iSpecies] / (rho*Mj);

      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
          * pow(Tve, Omega00[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
          + Omega00[iSpecies][jSpecies][1]*log(Tve)
          + Omega00[iSpecies][jSpecies][2]);

      /*--- Calculate "delta1_ij" ---*/
      d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;
    }

    /*--- Assign species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = gam_t*gam_t*Mi*(1-Mi*gam_i) / denom;
  }

  if (ionization) {
    iSpecies = nSpecies-1;
    /*--- Initialize the species diffusion coefficient ---*/
    DiffusionCoeff[iSpecies] = 0.0;

    /*--- Calculate molar concentration ---*/
    Mi      = Ms[iSpecies];
    gam_i   = V[RHOS_INDEX+iSpecies] / (rho*Mi);

    denom = 0.0;
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      if (iSpecies != jSpecies) {
        Mj    = GetMolar_Mass(jSpecies);
        gam_j = V[RHOS_INDEX+iSpecies] / (rho*Mj);

        /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
        Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
            * pow(Tve, Omega00[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
            + Omega00[iSpecies][jSpecies][1]*log(Tve)
            + Omega00[iSpecies][jSpecies][2]);

        /*--- Calculate "delta1_ij" ---*/
        d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;

        /*--- Calculate heavy-particle binary diffusion coefficient ---*/
        D_ij = kb*Tve/(P*d1_ij);
        denom += gam_j/D_ij;
      }
    }
    DiffusionCoeff[iSpecies] = gam_t*gam_t*Ms[iSpecies]*(1-Ms[iSpecies]*gam_i)
        / denom;
  }
}

void CTNE2Gas::Calc_Viscosity_GY(su2double *V) {

  unsigned short iSpecies, jSpecies;
  su2double rho, T, Tve;
  su2double *Ms, Mi, Mj, pi, Ru, RuSI, Na, gam_i, gam_j, denom;
  su2double Omega_ij, d2_ij;

  /*--- Extract constants ---*/
  pi   = PI_NUMBER;
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  Na   = AVOGAD_CONSTANT;

  /*--- Rename for convenience ---*/
  Ms   = Molar_Mass;
  rho  = V[RHO_INDEX];
  T    = V[T_INDEX];
  Tve  = V[TVE_INDEX];

  /*--- Initialize viscosity ---*/
  LamVisc = 0.0;

  /*--- Mixture viscosity via Gupta-Yos approximation ---*/
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    denom = 0.0;

    /*--- Calculate molar concentration ---*/
    Mi    = Ms[iSpecies];
    gam_i = V[RHOS_INDEX+iSpecies] / (rho*Mi);
    for (jSpecies = 0; jSpecies < nHeavy; jSpecies++) {
      Mj    = Ms[jSpecies];
      gam_j = V[RHOS_INDEX+jSpecies] / (rho*Mj);

      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
                       * pow(T, Omega11[iSpecies][jSpecies][0]*log(T)*log(T)
                       + Omega11[iSpecies][jSpecies][1]*log(T)
                       + Omega11[iSpecies][jSpecies][2]);
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;

      /*--- Add to denominator of viscosity ---*/
      denom += gam_j*d2_ij;
    }
    if (ionization) {
      jSpecies = nSpecies-1;
      Mj    = Ms[jSpecies];
      gam_j = V[RHOS_INDEX+jSpecies] / (rho*Mj);

      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
          * pow(Tve, Omega11[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
          + Omega11[iSpecies][jSpecies][1]*log(Tve)
          + Omega11[iSpecies][jSpecies][2]);
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;

      denom += gam_j*d2_ij;
    }

    /*--- Calculate species laminar viscosity ---*/
    LamVisc += (Mi/Na * gam_i) / denom;
  }

  if (ionization) {
    iSpecies = nSpecies-1;
    denom = 0.0;

    /*--- Calculate molar concentration ---*/
    Mi    = Ms[iSpecies];
    gam_i = V[RHOS_INDEX+iSpecies] / (rho*Mi);
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      Mj    = Ms[jSpecies];
      gam_j = V[RHOS_INDEX+jSpecies] / (rho*Mj);

      /*--- Calculate "delta" quantities ---*/
      Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
          * pow(Tve, Omega11[iSpecies][jSpecies][0]*log(Tve)*log(Tve)
          + Omega11[iSpecies][jSpecies][1]*log(Tve)
          + Omega11[iSpecies][jSpecies][2]);
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*Tve*(Mi+Mj))) * Omega_ij;

      /*--- Add to denominator of viscosity ---*/
      denom += gam_j*d2_ij;
    }
    LamVisc += (Mi/Na * gam_i) / denom;
  }

  V[LAM_VISC_INDEX] = LamVisc;
}

void CTNE2Gas::Calc_ThermalConductivity_GY(su2double *V) {

  unsigned short iSpecies, jSpecies;
  su2double rho, T, Tve, Cvve;
  su2double *xi, *Ms, Mi, Mj, mi, mj, pi, R, RuSI, Ru, Na, kb, gam_i, gam_j, Theta_v;
  su2double denom_t, denom_r, d1_ij, d2_ij, a_ij;
  su2double Omega_ij;

  if (ionization) {
    cout << "SetThermalConductivity: NEEDS REVISION w/ IONIZATION" << endl;
    exit(1);
  }

  /*--- Extract constants ---*/
  pi   = PI_NUMBER;
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  Na   = AVOGAD_CONSTANT;
  kb   = BOLTZMANN_CONSTANT;

  /*--- Rename for convenience ---*/
  Ms   = Molar_Mass;
  xi   = GetRotationModes();
  rho  = V[RHO_INDEX];
  T    = V[T_INDEX];
  Tve  = V[TVE_INDEX];
  Cvve = V[RHOCVVE_INDEX]/rho;


  /*--- Calculate mixture gas constant ---*/
  R = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    R += Ru * V[RHOS_INDEX+iSpecies]/rho;
  }

  /*--- Mixture thermal conductivity via Gupta-Yos approximation ---*/
  ThermalCond    = 0.0;
  ThermalCond_ve = 0.0;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {

    /*--- Calculate molar concentration ---*/
    Mi      = Ms[iSpecies];
    mi      = Mi/Na;
    gam_i   = V[RHOS_INDEX+iSpecies] / (rho*Mi);
    Theta_v = GetCharVibTemp(iSpecies);

    denom_t = 0.0;
    denom_r = 0.0;
    for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      Mj    = GetMolar_Mass(jSpecies);
      mj    = Mj/Na;
      gam_j = V[RHOS_INDEX+iSpecies] / (rho*Mj);

      a_ij = 1.0 + (1.0 - mi/mj)*(0.45 - 2.54*mi/mj) / ((1.0 + mi/mj)*(1.0 + mi/mj));

      /*--- Calculate the Omega^(0,0)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega00[iSpecies][jSpecies][3]
          * pow(T, Omega00[iSpecies][jSpecies][0]*log(T)*log(T)
          + Omega00[iSpecies][jSpecies][1]*log(T)
          + Omega00[iSpecies][jSpecies][2]);

      /*--- Calculate "delta1_ij" ---*/
      d1_ij = 8.0/3.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;

      /*--- Calculate the Omega^(1,1)_ij collision cross section ---*/
      Omega_ij = 1E-20 * Omega11[iSpecies][jSpecies][3]
          * pow(T, Omega11[iSpecies][jSpecies][0]*log(T)*log(T)
          + Omega11[iSpecies][jSpecies][1]*log(T)
          + Omega11[iSpecies][jSpecies][2]);

      /*--- Calculate "delta2_ij" ---*/
      d2_ij = 16.0/5.0 * sqrt((2.0*Mi*Mj) / (pi*Ru*T*(Mi+Mj))) * Omega_ij;

      denom_t += a_ij*gam_j*d2_ij;
      denom_r += gam_j*d1_ij;
    }

    /*--- Translational contribution to thermal conductivity ---*/
    ThermalCond    += (15.0/4.0)*kb*gam_i/denom_t;

    /*--- Translational contribution to thermal conductivity ---*/
    if (xi[iSpecies] != 0.0)
      ThermalCond  += kb*gam_i/denom_r;

    /*--- Vibrational-electronic contribution to thermal conductivity ---*/
    ThermalCond_ve += kb*Cvve/R*gam_i / denom_r;
  }
}

CMutationGas::CMutationGas(string Mixfile, string Transport, unsigned short val_nSpecies,
                           unsigned short val_nDim, bool val_ionization, bool val_viscous):
  CMultiSpeciesGas(val_nSpecies, val_nDim, val_ionization, val_viscous)//,
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

