/*!
 * \file fluid_model.hpp
 * \brief Headers of the main thermodynamic subroutines of the SU2 solvers.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna, T. Economon
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

#include "../../Common/include/mpi_structure.hpp"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#define LEN_COMPONENTS 32

#include "stdio.h"
#include "math.h"

using namespace std;

#include "../include/transport_model.hpp"
#include "../../Common/include/config_structure.hpp"
#include "mutation.hpp"

#ifdef HAVE_MUTATIONPP
#include "mutation++.h"
#endif

/*!
 * \class CFluidModel
 * \brief Main class for defining the Thermo-Physical Model
 * a child class for each particular Model (Ideal-Gas, Van der Waals, etc.)
 * \author: S.Vitale, G.Gori, M.Pini, T. Economon
 */
class CFluidModel {
protected:
  su2double
  StaticEnergy,  /*!< \brief Internal Energy. */
  Entropy,       /*!< \brief Entropy. */
  Density,       /*!< \brief Density. */
  Pressure,      /*!< \brief Pressure. */
  SoundSpeed2,   /*!< \brief SpeedSound. */
  Temperature,   /*!< \brief Temperature. */
  dPdrho_e,      /*!< \brief DpDd_e. */
  dPde_rho,      /*!< \brief DpDe_d. */
  dTdrho_e,      /*!< \brief DTDd_e. */
  dTde_rho,      /*!< \brief DTDe_d. */
  dhdrho_P,	     /*!< \brief DhDrho_p. */
  dhdP_rho,	     /*!< \brief DhDp_rho. */
  dsdrho_P,	     /*!< \brief DsDrho_p. */
  dsdP_rho,	     /*!< \brief DsDp_rho. */
  Cp,            /*!< \brief Specific Heat Capacity at constant pressure. */
  Cv,            /*!< \brief Specific Heat Capacity at constant volume. */
  Mu,            /*!< \brief Laminar viscosity. */
  Mu_Turb,       /*!< \brief Eddy viscosity provided by a turbulence model (RANS). */
  dmudrho_T,     /*!< \brief Specific Heat Capacity at constant pressure. */
  dmudT_rho,     /*!< \brief Specific Heat Capacity at constant pressure. */
  Kt,            /*!< \brief Specific Heat Capacity at constant pressure. */
  dktdrho_T,     /*!< \brief Specific Heat Capacity at constant pressure. */
  dktdT_rho;     /*!< \brief Specific Heat Capacity at constant pressure. */

  unsigned short RHO_INDEX, RHOS_INDEX,     /*!< \brief Index for TNE2 Variables */
  T_INDEX, TVE_INDEX, VEL_INDEX,
  P_INDEX, RHOCVTR_INDEX,
  RHOCVVE_INDEX, A_INDEX, H_INDEX,
  LAM_VISC_INDEX, EDDY_VISC_INDEX,
  DIFF_COEFF_INDEX, K_INDEX,
  KVE_INDEX;

  CViscosityModel    *LaminarViscosity;     /*!< \brief Laminar Viscosity Model */
  CConductivityModel *ThermalConductivity;  /*!< \brief Thermal Conductivity Model */

public:

  /*!
   * \brief Constructor of the class.
   */
  CFluidModel(void);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CFluidModel(void);

  /*!
   * \brief Get fluid pressure.
   */
  su2double GetPressure ();

  /*!
   * \brief Get fluid temperature.
   */
  su2double GetTemperature ();

  /*!
   * \brief Get fluid entropy.
   */
  su2double GetEntropy ();

  /*!
   * \brief Get fluid internal energy.
   */
  su2double GetStaticEnergy ();

  /*!
   * \brief Get fluid density.
   */
  su2double GetDensity ();

  /*!
   * \brief Get fluid speed of sound.
   */
  su2double GetSoundSpeed ();

  /*!
   * \brief Get fluid speed of sound squared.
   */
  su2double GetSoundSpeed2 ();

  /*!
   * \brief Get fluid specific heat at constant pressure.
   */
  su2double GetCp ();

  /*!
   * \brief Get fluid specific heat at constant volume.
   */
  su2double GetCv ();

  /*!
   * \brief Get fluid dynamic viscosity
   */
  virtual su2double GetLaminarViscosity ();

  /*!
   * \brief Get fluid thermal conductivity
   */
  virtual su2double GetThermalConductivity ();

  /*!
   * \brief A virtual member.
   * \return Value of the species diffusion coefficient.
   */
  inline su2double* GetDiffusionCoeff(void) {return NULL; }

  /*!
   * \brief Get fluid pressure partial derivative.
   */
  su2double GetdPdrho_e ();

  /*!
   * \brief Get fluid pressure partial derivative.
   */
  su2double GetdPde_rho ();

  /*!
   * \brief Get fluid temperature partial derivative.
   */
  su2double GetdTdrho_e ();

  /*!
   * \brief Get fluid temperature partial derivative.
   */
  su2double GetdTde_rho ();

  /*!
   * \brief Get fluid pressure partial derivative.
   */
  su2double Getdhdrho_P ();

  /*!
   * \brief Get fluid pressure partial derivative.
   */
  su2double GetdhdP_rho ();

  /*!
   * \brief Get fluid temperature partial derivative.
   */
  su2double Getdsdrho_P ();

  /*!
   * \brief Get fluid temperature partial derivative.
   */
  su2double GetdsdP_rho ();

  /*!
   * \brief Get fluid dynamic viscosity partial derivative.
   */
  su2double Getdmudrho_T ();

  /*!
   * \brief Get fluid dynamic viscosity partial derivative.
   */
  su2double GetdmudT_rho ();

  /*!
   * \brief Get fluid thermal conductivity partial derivative.
   */
  su2double Getdktdrho_T ();

  /*!
   * \brief Get fluid thermal conductivity partial derivative.
   */
  su2double GetdktdT_rho ();

  /*!
   * \brief Set specific heat Cp model.
   */
  virtual void SetCpModel (CConfig *config);

  /*!
   * \brief Set viscosity model.
   */
  void SetLaminarViscosityModel (CConfig *config);

  /*!
   * \brief Set thermal conductivity model.
   */
  void SetThermalConductivityModel (CConfig *config);

  /*!
   * \brief Set general transport model for MultiSpecies.
   */
  virtual void Calc_TransportCoeff(CConfig *config, su2double *V);


  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  virtual void SetTDState_rhoe (su2double rho, su2double e );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (T).
   */
  virtual void SetTDState_PT (su2double P, su2double T );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   */
  virtual void SetTDState_Prho (su2double P, su2double rho );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   */
  virtual void SetEnergy_Prho (su2double P, su2double rho );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("hs").
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  virtual void SetTDState_hs (su2double h, su2double s );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
   * \param[in] th1 - first thermodynamic variable (rho).
   * \param[in] th2 - second thermodynamic variable (T).
   */
  virtual void SetTDState_rhoT (su2double rho, su2double T );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  virtual void SetTDState_Ps (su2double P, su2double s );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   */
  virtual void ComputeDerivativeNRBC_Prho (su2double P, su2double rho );

  /*!
   * \brief Virtual member.
   * \param[in] T - Temperature value at the point.
   */
  virtual void SetTDState_T(su2double val_Temperature);

  /*!
   * \brief Set fluid eddy viscosity provided by a turbulence model needed for computing effective thermal conductivity.
   */
  void SetEddyViscosity(su2double val_Mu_Turb);

  /*!
   * \brief Get the flow ionization.
   * \return Flow Ionization (yes or no).
   */
  virtual bool GetIonization(void);

  /*!
   * \brief Get the vector of free stream mass fraction values.
   * \return Ratio of species mass to mixture mass.
   */
  virtual su2double* GetMassFrac_FreeStream(void);

  /*!
   * \brief Provides the degeneracy of electron states for calculating e_el
   * \return: Vector of characteristic vibrational temperatures [K]
   */
  virtual su2double** GetElDegeneracy(void);

  /*!
   * \brief Provides number electron states for calculating e_el
   * \return: Vector of number of electron states for each species
   */
  virtual unsigned short* GetnElStates(void);

  /*!
   * \brief Provides the number of chemical reactions in the chemistry model
   * \return: The number of chemical reactions, read from input file
   */
  virtual unsigned short GetnReactions(void);

  /*!
   * \brief Provides the number of chemical reactions in the chemistry model
   * \return: The number of chemical reactions, read from input file
   */
  virtual su2double GetArrheniusCoeff(unsigned short iReaction);

  /*!
   * \brief Provides the number of chemical reactions in the chemistry model
   * \return: The number of chemical reactions, read from input file
   */
  virtual su2double GetArrheniusEta(unsigned short iReaction);

  /*!
   * \brief Provides the number of chemical reactions in the chemistry model
   * \return: The number of chemical reactions, read from input file
   */
  virtual su2double GetArrheniusTheta(unsigned short iReaction);

  /*!
   * \brief Provides the rate controlling temperature exponents for chemistry.
   * \return: Rate controlling temperature exponents.
   */
  virtual su2double* GetRxnTcf_a(void);

  /*!
   * \brief Provides the rate controlling temperature exponents for chemistry.
   * \return: Rate controlling temperature exponents.
   */
  virtual su2double* GetRxnTcf_b(void);

  /*!
   * \brief Provides the rate controlling temperature exponents for chemistry.
   * \return: Rate controlling temperature exponents.
   */
  virtual su2double* GetRxnTcb_a(void);

  /*!
   * \brief Provides the rate controlling temperature exponents for chemistry.
   * \return: Rate controlling temperature exponents.
   */
  virtual su2double* GetRxnTcb_b(void);

  /*!
   * \brief Dissociation potential of species.
   * \return: Dissociation potential.
   */
  virtual su2double* GetDissociationPot(void);

  /*!
   * \brief Provides the number of rotational modes of energy storage
   * \return: Vector of rotational mode count
   */
  virtual su2double* GetRotationModes(void);

  /*!
   * \brief Provides the characteristic vibrational temperature for calculating e_vib
   * \return: Vector of characteristic vibrational temperatures [K]
   */
  virtual su2double* GetCharVibTemp(void);

  /*!
   * \brief Provides the characteristic electronic temperature for calculating e_el
   * \return: Vector of characteristic vibrational temperatures [K]
   */
  virtual su2double** GetCharElTemp(void);

  /*!
   * \brief Provides the thermodynamic reference temperatures from the JANAF tables
   * \return: Vector of reference temperatures [K]
   */
  virtual su2double* GetRefTemperature(void);

  /*!
   * \brief Provides the characteristic vibrational temperature for calculating e_vib
   * \return: The number of chemical reactions, read from input file
   */
  virtual su2double GetCharVibTemp(unsigned short iSpecies);

  /*!
   * \brief Provides a table of equilibrium constants for a particular chemical reaction for a supplied gas model.
   * \return: Matrix of reaction constants
   */
  virtual void GetChemistryEquilConstants(CConfig *config, su2double **RxnConstantTable, unsigned short iReaction);

  /*!
   * \brief Provides the molar mass of each species present in multi species fluid
   * \return: Vector of molar mass of each species in kg/kmol
   */
  virtual su2double* GetMolar_Mass(void);

  /*!
   * \brief Provides the molar mass of each species present in multi species fluid
   * \return: Mass of each species in Kg
   */
  virtual su2double GetMolar_Mass(unsigned short iSpecies);

  /*!
   * \brief Provides the formation enthalpy of the specified species at standard conditions
   * \return: Enthalpy of formation
   */
  virtual su2double* GetEnthalpy_Formation(void);

  /*!
   * \brief Provides the formation enthalpy of the specified species at standard conditions
   * \return: Enthalpy of formation
   */
  virtual su2double GetEnthalpy_Formation(unsigned short iSpecies);

  /*!
   * \brief Get the array that maps chemical consituents to each chemical reaction.
   * \return Memory location of the triple pointer to the 3-D reaction map array.
   */
  virtual int ***GetReaction_Map(void);

  /*!
   * \brief Get the array containing the curve fit coefficients for the Omega(0,0) collision integrals.
   * \return Memory location of the triple pointer to the 3-D collision integral array.
   */
  virtual su2double ***GetCollisionIntegral00(void);

  /*!
   * \brief Get the array containing the curve fit coefficients for the Omega(1,1) collision integrals.
   * \return Memory location of the triple pointer to the 3-D collision integral array.
   */
  virtual su2double ***GetCollisionIntegral11(void);

  /*!
   * \brief Get the coefficients of the Blottner viscosity model
   * \param[in] val_Species - Index of the species
   * \param[in] val_Coeff - Index of the coefficient (As, Bs, Cs)
   * \return Value of the Blottner coefficient
   */
  virtual su2double **GetBlottnerCoeff(void);

  /*!
   * \brief Get the wall heat flux on a constant heat flux boundary.
   * \return The heat flux.
   */
  virtual su2double *GetWall_Catalycity(void);

  /*!
   * \brief Calculate the Vibrational-Electronic Energy for a species
   */
  virtual su2double CalcEve(su2double val_Tve, unsigned short val_Species);

  /*!
   * \brief Calculate density of the gas
   */
  virtual su2double Calc_Density(su2double *MassFrac,
                                 su2double T, su2double Tve, su2double P);

  virtual su2double Calc_CvVibElSpecies(su2double val_Tve, unsigned short val_Species);
  virtual su2double Calc_CvTraRotSpecies(su2double *Ms, su2double Ru, unsigned short val_Species);

  virtual su2double Calc_SoundSpeed(su2double *cs, su2double rhoCvtr,
                                    su2double rho, su2double P);

  virtual su2double Calc_Enthalpies(su2double val_T, su2double val_eves, unsigned short val_Species);

  virtual su2double  Calc_MixtureEnergy(su2double* cs, su2double sqvel,
                                        su2double rho,
                                        su2double T, su2double Tve);

  virtual void InitializeMixture(CConfig *config);

  /*!
   * \brief A virtual member.
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dPdU
   */
  virtual void CalcdPdU(su2double *V, su2double *val_eves, su2double *dPdU);

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dTdU
   */
  virtual void CalcdTdU(su2double *V, su2double *dTdU);

  /*!
   * \brief Set partial derivative of temperature w.r.t. density \f$\frac{\partial P}{\partial \rho_s}\f$
   * \param[in] V
   * \param[in] config - Configuration settings
   * \param[in] dTdU
   */
  virtual void CalcdTvedU(su2double *V, su2double *val_eves, su2double *dTdU);

};

/*!
 * \class CIdealGas
 * \brief Child class for defining ideal gas model.
 * \author: S.Vitale, M.Pini.
 */
class CIdealGas : public CFluidModel {

protected:
  su2double Gamma,        /*!< \brief Heat Capacity Ratio. */
  Gamma_Minus_One,        /*!< \brief Heat Capacity Ratio Minus One. */
  Gas_Constant;           /*!< \brief Gas Constant. */

  bool  ComputeEntropy;   /*!< \brief Whether or not to compute entropy. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CIdealGas(void);

  /*!
   * \brief Constructor of the class.
   */
  CIdealGas(su2double gamma, su2double R);

  /*!
   * \brief Constructor of the class.
   */
  CIdealGas(su2double gamma, su2double R, bool CompEntropy);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CIdealGas(void);

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  void SetTDState_rhoe (su2double rho, su2double e );

  /*!
   * \brief Set the Dimensionless State using Pressure  and Temperature
   * \param[in] P - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  void SetTDState_PT (su2double P, su2double T );

  /*!
   * \brief Set the Dimensionless State using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetTDState_Prho (su2double P, su2double rho );

  /*!
   * \brief Set the Dimensionless Internal Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetEnergy_Prho (su2double P, su2double rho );

  /*!
   * \brief Set the Dimensionless State using Enthalpy and Entropy
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  void SetTDState_hs (su2double h, su2double s );

  /*!
   * \brief Set the Dimensionless State using Density and Temperature
   * \param[in] th1 - first thermodynamic variable (rho).
   * \param[in] th2 - second thermodynamic variable (T).
   */
  void SetTDState_rhoT (su2double rho, su2double T );

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  void SetTDState_Ps (su2double P, su2double s );

  /*!
   * \brief compute some derivatives of enthalpy and entropy needed for subsonic inflow BC
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   */
  void ComputeDerivativeNRBC_Prho (su2double P, su2double rho );
};

/*!
 * derived class CVanDerWaalsGas
 * \brief Child class for defining the Van der Waals model.
 * \author: S.Vitale, M.Pini
 */
class CVanDerWaalsGas : public CIdealGas {

protected:
  su2double
  a, b, Zed;  /*!< \brief Parameters for the Dimensionless Equation. */

public:

  /*!
     * \brief Constructor of the class.
     */
  CVanDerWaalsGas(void);

  /*!
   * \brief Constructor of the class.
   */
  CVanDerWaalsGas(su2double gamma, su2double R, su2double Pstar, su2double Tstar);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CVanDerWaalsGas(void);

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  void SetTDState_rhoe (su2double rho, su2double e );

  /*!
   * \brief Set the Dimensionless State using Pressure and Temperature
   * \param[in] P - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  void SetTDState_PT (su2double P, su2double T );

  /*!
   * \brief Set the Dimensionless State using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetTDState_Prho (su2double P, su2double rho );

  /*!
   * \brief Set the Dimensionless Internal Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetEnergy_Prho (su2double P, su2double rho );

  /*!
   * \brief Set the Dimensionless state using Enthalpy and Entropy
   * \param[in] h - first thermodynamic variable (h).
   * \param[in] s - second thermodynamic variable (s).
   */
  void SetTDState_hs (su2double h, su2double s );

  /*!
   * \brief Set the Dimensionless state using Density and Temperature
   * \param[in] rho - first thermodynamic variable (rho).
   * \param[in] T - second thermodynamic variable (T).
   */
  void SetTDState_rhoT (su2double rho, su2double T );

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   * \param[in] P - first thermodynamic variable (P).
   * \param[in] s - second thermodynamic variable (s).
   */
  void SetTDState_Ps (su2double P, su2double s );

  /*!
   * \brief compute some derivatives of enthalpy and entropy needed for subsonic inflow BC
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   */
  void ComputeDerivativeNRBC_Prho (su2double P, su2double rho );

};

/*!
 * \derived class CPengRobinson
 * \brief Child class for defining the Peng-Robinson model.
 * \author: S.Vitale, G. Gori
 */
class CPengRobinson : public CIdealGas {

protected:
  su2double a,  /*!< \brief model parameter. */
  b,            /*!< \brief model parameter. */
  k,            /*!< \brief model parameter (computed with acentric factor). */
  Zed,          /*!< \brief compressibility factor. */
  TstarCrit;    /*!< \brief Critical temperature. */

private:

  /*!
   * \brief Internal model parameter.
   */
  su2double  alpha2 (su2double T);

  /*!
   * \brief Internal function for the implicit call hs.
   */
  su2double  T_v_h (su2double v, su2double h);

  /*!
   * \brief Internal function for the implicit call Ps.
   */
  su2double T_P_rho(su2double P, su2double rho);

public:

  /*!
   * \brief Constructor of the class.
   */
  CPengRobinson(void);

  /*!
   * \brief Constructor of the class.
   */
  CPengRobinson(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CPengRobinson(void);

  /*!
   * \brief Set the Dimensionless State using Density and Internal Energy
   * \param[in] rho - first thermodynamic variable.
   * \param[in] e - second thermodynamic variable.
   */
  void SetTDState_rhoe (su2double rho, su2double e );

  /*!
   * \brief Set the Dimensionless State using Pressure and Temperature
   * \param[in] P - first thermodynamic variable.
   * \param[in] T - second thermodynamic variable.
   */
  void SetTDState_PT (su2double P, su2double T );

  /*!
   * \brief Set the Dimensionless State using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetTDState_Prho (su2double P, su2double rho );

  /*!
   * \brief Set the Dimensionless Energy using Pressure and Density
   * \param[in] P - first thermodynamic variable.
   * \param[in] rho - second thermodynamic variable.
   */
  void SetEnergy_Prho (su2double P, su2double rho );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("hs").
   * \param[in] th1 - first thermodynamic variable (h).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  void SetTDState_hs (su2double h, su2double s );

  /*!
   * \brief virtual member that would be different for each gas model implemented
   * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
   * \param[in] th1 - first thermodynamic variable (rho).
   * \param[in] th2 - second thermodynamic variable (T).
   */
  void SetTDState_rhoT (su2double rho, su2double T );

  /*!
   * \brief Set the Dimensionless State using Pressure and Entropy
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (s).
   */
  void SetTDState_Ps (su2double P, su2double s );

  /*!
   * \brief compute some derivatives of enthalpy and entropy needed for subsonic inflow BC
   * \param[in] InputSpec - Input pair for FLP calls ("Pv").
   * \param[in] th1 - first thermodynamic variable (P).
   * \param[in] th2 - second thermodynamic variable (v).
   */
  void ComputeDerivativeNRBC_Prho (su2double P, su2double rho );

};

/*!
 * \derived class CMultiSpeciesGas
 * \brief Child class for defining the two-temperature model.
 * \author: W. Maier
 */
class CMultiSpeciesGas : public CFluidModel {

protected:

  bool ionization;          /*!< \brief Presence of charged species in gas mixture. */
  bool viscous;             /*!< \brief Presence of viscous effects. */
  unsigned short nSpecies;  /*!< \brief Number of species in the gas mixture. */
  unsigned short nDim;      /*!< \brief Number of dimensions. */
  unsigned short nHeavy;    /*!< \brief Number of heavy particles in gas */
  unsigned short nEl;       /*!< \brief Number of electrons in gas */

  unsigned iSpecies;        /*!< \brief Common iteration counter for species */
  unsigned jSpecies;        /*!< \brief Common iteration counter for species */
  unsigned iDim;            /*!< \brief Common iteration counter for dimensions */

  int ***Reactions;           /*!</brief reaction map for chemically reacting flows */
  su2double *Molar_Mass;      /*!< \brief Molar mass of the multi-species fluid [kg/kmol] */
  su2double ***Omega00;       /*!< \brief Collision integrals (Omega(0,0)) */
  su2double ***Omega11;       /*!< \brief Collision integrals (Omega(1,1)) */
  unsigned short *nElStates;  /*!< \brief Number of electron states. */
  su2double **CharElTemp,     /*!< \brief Characteristic temperature of electron states. */
  **ElDegeneracy;             /*!< \brief Degeneracy of electron states. */
  su2double *ArrheniusCoefficient,	/*!< \brief Arrhenius reaction coefficient */
  *ArrheniusEta,					/*!< \brief Arrhenius reaction temperature exponent */
  *ArrheniusTheta,					/*!< \brief Arrhenius reaction characteristic temperature */
  *CharVibTemp,						/*!< \brief Characteristic vibrational temperature for e_vib */
  *RotationModes,			        /*!< \brief Rotational modes of energy storage */
  *Ref_Temperature,   			    /*!< \brief Reference temperature for thermodynamic relations */
  *Tcf_a,                     /*!< \brief Rate controlling temperature exponent (fwd) */
  *Tcf_b,                     /*!< \brief Rate controlling temperature exponent (fwd) */
  *Tcb_a,                     /*!< \brief Rate controlling temperature exponent (bkw) */
  *Tcb_b,                     /*!< \brief Rate controlling temperature exponent (bkw) */
  *Diss;                      /*!< \brief Dissociation potential. */
  unsigned short nReactions;         /*!< \brief Number of reactions in chemical model. */
  su2double   *MassFrac_FreeStream;  /*!< \brief Mixture mass fractions of the fluid. */
  su2double   *Enthalpy_Formation;   /*!< \brief Enthalpy of formation */
  su2double *Wall_Catalycity;        /*!< \brief Specified wall species mass-fractions for catalytic boundaries. */
  su2double *Particle_Mass,          /*!< \brief Mass of all particles present in the plasma */
  Mixture_Molar_mass,                /*!< \brief Molar mass of the multi-species fluid [kg/kmol] */
  **Blottner,                        /*!< \brief Blottner viscosity coefficients */
  *Species_Ref_Temperature,          /*!< \brief Reference Temperature for viscosity of all particles present in the plasma */
  *Species_Ref_Viscosity;            /*!< \brief Reference viscosity  of all particles present in the plasma */
  su2double *DiffusionCoeff;  /*!< \brief Diffusion coefficient of the mixture. */
  su2double **Dij;            /*!< \brief Binary diffusion coefficients. */
  su2double LamVisc;          /*!< \brief Viscosity of the fluid. */
  su2double ThermalCond;      /*!< \brief T-R thermal conductivity of the gas mixture. */
  su2double ThermalCond_ve;   /*!< \brief V-E thermal conductivity of the gas mixture. */

private:

public:

  /*!
   * \brief Constructor of the class.
   */
  CMultiSpeciesGas(unsigned short val_nSpecies, unsigned short val_nDim, bool val_ionization, bool viscous);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CMultiSpeciesGas(void);

  /*!
   * \brief Initialize fluid type/compostion used.
   * \param[in] config - Descritpion of the problem.
   */
  void InitializeMixture(CConfig *config);

  /*!
   * \brief Get the species diffusion coefficient.
   * \return Value of the species diffusion coefficient.
   */
  inline su2double* GetDiffusionCoeff(void) { return DiffusionCoeff; }

  /*!
   * \brief Get the laminar viscosity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetLaminarViscosity(void) { return LamVisc; }

  /*!
   * \brief Get the thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity(void) { return ThermalCond; }

  /*!
   * \brief Get the vib-el. thermal conductivity of the flow.
   * \return Value of the laminar viscosity of the flow.
   */
  inline su2double GetThermalConductivity_ve(void) { return ThermalCond_ve; }


  /*!
   * \brief Initialize fluid type/compostion used.
   * \param[in] config - Descritpion of the problem.
   */
  inline su2double* GetMassFrac_FreeStream(void) { return MassFrac_FreeStream;}

  unsigned short GetnReactions(void) {return nReactions; }

  su2double GetArrheniusCoeff(unsigned short iReaction) { return ArrheniusCoefficient[iReaction]; }

  su2double GetArrheniusEta(unsigned short iReaction) { return ArrheniusEta[iReaction]; }

  su2double GetArrheniusTheta(unsigned short iReaction) { return ArrheniusTheta[iReaction]; }

  su2double* GetRxnTcf_a(void) { return Tcf_a; }

  su2double* GetRxnTcf_b(void) { return Tcf_b; }

  su2double* GetRxnTcb_a(void) { return Tcb_a; }

  su2double* GetRxnTcb_b(void) { return Tcb_b; }

  su2double* GetDissociationPot(void) { return Diss; }

  su2double GetCharVibTemp(unsigned short iSpecies) {return CharVibTemp[iSpecies]; }

  su2double* GetCharVibTemp() {return CharVibTemp; }

  su2double** GetCharElTemp() {return CharElTemp; }

  unsigned short* GetnElStates() {return nElStates; }

  su2double** GetElDegeneracy() {return ElDegeneracy; }

  su2double* GetRotationModes() { return RotationModes; }

  su2double* GetRefTemperature() { return Ref_Temperature; }

  su2double* GetWall_Catalycity() { return Wall_Catalycity; }

  inline su2double* GetMolar_Mass() { return Molar_Mass; }

  inline su2double GetMolar_Mass(unsigned short iSpecies) { return Molar_Mass[iSpecies]; }

  su2double* GetEnthalpy_Formation(void) { return Enthalpy_Formation; }

  su2double GetEnthalpy_Formation(unsigned short iSpecies) { return Enthalpy_Formation[iSpecies]; }

  bool GetIonization(void) { return ionization; }

  int ***GetReaction_Map(void) { return Reactions; }

  su2double ***GetCollisionIntegral00(void) { return Omega00; }

  su2double ***GetCollisionIntegral11(void) { return Omega11; }

  su2double **GetBlottnerCoeff(void) { return Blottner; }

  void GetChemistryEquilConstants(CConfig *config, su2double **RxnConstantTable, unsigned short iReaction) {

    switch (config->GetKind_GasModel()) {

    case O2:

      //O2 + M -> 2O + M
      RxnConstantTable[0][0] = 1.8103;  RxnConstantTable[0][1] = 1.9607;  RxnConstantTable[0][2] = 3.5716;  RxnConstantTable[0][3] = -7.3623;   RxnConstantTable[0][4] = 0.083861;
      RxnConstantTable[1][0] = 0.91354; RxnConstantTable[1][1] = 2.3160;  RxnConstantTable[1][2] = 2.2885;  RxnConstantTable[1][3] = -6.7969;   RxnConstantTable[1][4] = 0.046338;
      RxnConstantTable[2][0] = 0.64183; RxnConstantTable[2][1] = 2.4253;  RxnConstantTable[2][2] = 1.9026;  RxnConstantTable[2][3] = -6.6277;   RxnConstantTable[2][4] = 0.035151;
      RxnConstantTable[3][0] = 0.55388; RxnConstantTable[3][1] = 2.4600;  RxnConstantTable[3][2] = 1.7763;  RxnConstantTable[3][3] = -6.5720;   RxnConstantTable[3][4] = 0.031445;
      RxnConstantTable[4][0] = 0.52455; RxnConstantTable[4][1] = 2.4715;  RxnConstantTable[4][2] = 1.7342;  RxnConstantTable[4][3] = -6.55534;  RxnConstantTable[4][4] = 0.030209;
      RxnConstantTable[5][0] = 0.50989; RxnConstantTable[5][1] = 2.4773;  RxnConstantTable[5][2] = 1.7132;  RxnConstantTable[5][3] = -6.5441;   RxnConstantTable[5][4] = 0.029591;

      break;

    case N2:

      //N2 + M -> 2N + M
      RxnConstantTable[0][0] = 3.4907;  RxnConstantTable[0][1] = 0.83133; RxnConstantTable[0][2] = 4.0978;  RxnConstantTable[0][3] = -12.728; RxnConstantTable[0][4] = 0.07487;   //n = 1E14
      RxnConstantTable[1][0] = 2.0723;  RxnConstantTable[1][1] = 1.38970; RxnConstantTable[1][2] = 2.0617;  RxnConstantTable[1][3] = -11.828; RxnConstantTable[1][4] = 0.015105;  //n = 1E15
      RxnConstantTable[2][0] = 1.6060;  RxnConstantTable[2][1] = 1.57320; RxnConstantTable[2][2] = 1.3923;  RxnConstantTable[2][3] = -11.533; RxnConstantTable[2][4] = -0.004543; //n = 1E16
      RxnConstantTable[3][0] = 1.5351;  RxnConstantTable[3][1] = 1.60610; RxnConstantTable[3][2] = 1.2993;  RxnConstantTable[3][3] = -11.494; RxnConstantTable[3][4] = -0.00698;  //n = 1E17
      RxnConstantTable[4][0] = 1.4766;  RxnConstantTable[4][1] = 1.62910; RxnConstantTable[4][2] = 1.2153;  RxnConstantTable[4][3] = -11.457; RxnConstantTable[4][4] = -0.00944;  //n = 1E18
      RxnConstantTable[5][0] = 1.4766;  RxnConstantTable[5][1] = 1.62910; RxnConstantTable[5][2] = 1.2153;  RxnConstantTable[5][3] = -11.457; RxnConstantTable[5][4] = -0.00944;  //n = 1E19

      break;

    case ARGON_SID:

      //N2 + M -> 2N + M
      RxnConstantTable[0][0] = 3.4907;  RxnConstantTable[0][1] = 0.83133; RxnConstantTable[0][2] = 4.0978;  RxnConstantTable[0][3] = -12.728; RxnConstantTable[0][4] = 0.07487;   //n = 1E14
      RxnConstantTable[1][0] = 2.0723;  RxnConstantTable[1][1] = 1.38970; RxnConstantTable[1][2] = 2.0617;  RxnConstantTable[1][3] = -11.828; RxnConstantTable[1][4] = 0.015105;  //n = 1E15
      RxnConstantTable[2][0] = 1.6060;  RxnConstantTable[2][1] = 1.57320; RxnConstantTable[2][2] = 1.3923;  RxnConstantTable[2][3] = -11.533; RxnConstantTable[2][4] = -0.004543; //n = 1E16
      RxnConstantTable[3][0] = 1.5351;  RxnConstantTable[3][1] = 1.60610; RxnConstantTable[3][2] = 1.2993;  RxnConstantTable[3][3] = -11.494; RxnConstantTable[3][4] = -0.00698;  //n = 1E17
      RxnConstantTable[4][0] = 1.4766;  RxnConstantTable[4][1] = 1.62910; RxnConstantTable[4][2] = 1.2153;  RxnConstantTable[4][3] = -11.457; RxnConstantTable[4][4] = -0.00944;  //n = 1E18
      RxnConstantTable[5][0] = 1.4766;  RxnConstantTable[5][1] = 1.62910; RxnConstantTable[5][2] = 1.2153;  RxnConstantTable[5][3] = -11.457; RxnConstantTable[5][4] = -0.00944;  //n = 1E19

      break;

    case AIR5:

      if (iReaction <= 4) {

        //N2 + M -> 2N + M
        RxnConstantTable[0][0] = 3.4907;  RxnConstantTable[0][1] = 0.83133; RxnConstantTable[0][2] = 4.0978;  RxnConstantTable[0][3] = -12.728; RxnConstantTable[0][4] = 0.07487;   //n = 1E14
        RxnConstantTable[1][0] = 2.0723;  RxnConstantTable[1][1] = 1.38970; RxnConstantTable[1][2] = 2.0617;  RxnConstantTable[1][3] = -11.828; RxnConstantTable[1][4] = 0.015105;  //n = 1E15
        RxnConstantTable[2][0] = 1.6060;  RxnConstantTable[2][1] = 1.57320; RxnConstantTable[2][2] = 1.3923;  RxnConstantTable[2][3] = -11.533; RxnConstantTable[2][4] = -0.004543; //n = 1E16
        RxnConstantTable[3][0] = 1.5351;  RxnConstantTable[3][1] = 1.60610; RxnConstantTable[3][2] = 1.2993;  RxnConstantTable[3][3] = -11.494; RxnConstantTable[3][4] = -0.00698;  //n = 1E17
        RxnConstantTable[4][0] = 1.4766;  RxnConstantTable[4][1] = 1.62910; RxnConstantTable[4][2] = 1.2153;  RxnConstantTable[4][3] = -11.457; RxnConstantTable[4][4] = -0.00944;  //n = 1E18
        RxnConstantTable[5][0] = 1.4766;  RxnConstantTable[5][1] = 1.62910; RxnConstantTable[5][2] = 1.2153;  RxnConstantTable[5][3] = -11.457; RxnConstantTable[5][4] = -0.00944;  //n = 1E19

      } else if (iReaction > 4 && iReaction <= 9) {

        //O2 + M -> 2O + M
        RxnConstantTable[0][0] = 1.8103;  RxnConstantTable[0][1] = 1.9607;  RxnConstantTable[0][2] = 3.5716;  RxnConstantTable[0][3] = -7.3623;   RxnConstantTable[0][4] = 0.083861;
        RxnConstantTable[1][0] = 0.91354; RxnConstantTable[1][1] = 2.3160;  RxnConstantTable[1][2] = 2.2885;  RxnConstantTable[1][3] = -6.7969;   RxnConstantTable[1][4] = 0.046338;
        RxnConstantTable[2][0] = 0.64183; RxnConstantTable[2][1] = 2.4253;  RxnConstantTable[2][2] = 1.9026;  RxnConstantTable[2][3] = -6.6277;   RxnConstantTable[2][4] = 0.035151;
        RxnConstantTable[3][0] = 0.55388; RxnConstantTable[3][1] = 2.4600;  RxnConstantTable[3][2] = 1.7763;  RxnConstantTable[3][3] = -6.5720;   RxnConstantTable[3][4] = 0.031445;
        RxnConstantTable[4][0] = 0.52455; RxnConstantTable[4][1] = 2.4715;  RxnConstantTable[4][2] = 1.7342;  RxnConstantTable[4][3] = -6.55534;  RxnConstantTable[4][4] = 0.030209;
        RxnConstantTable[5][0] = 0.50989; RxnConstantTable[5][1] = 2.4773;  RxnConstantTable[5][2] = 1.7132;  RxnConstantTable[5][3] = -6.5441;   RxnConstantTable[5][4] = 0.029591;

      } else if (iReaction > 9 && iReaction <= 14) {

        //NO + M -> N + O + M
        RxnConstantTable[0][0] = 2.1649;  RxnConstantTable[0][1] = 0.078577;  RxnConstantTable[0][2] = 2.8508;  RxnConstantTable[0][3] = -8.5422; RxnConstantTable[0][4] = 0.053043;
        RxnConstantTable[1][0] = 1.0072;  RxnConstantTable[1][1] = 0.53545;   RxnConstantTable[1][2] = 1.1911;  RxnConstantTable[1][3] = -7.8098; RxnConstantTable[1][4] = 0.004394;
        RxnConstantTable[2][0] = 0.63817; RxnConstantTable[2][1] = 0.68189;   RxnConstantTable[2][2] = 0.66336; RxnConstantTable[2][3] = -7.5773; RxnConstantTable[2][4] = -0.011025;
        RxnConstantTable[3][0] = 0.55889; RxnConstantTable[3][1] = 0.71558;   RxnConstantTable[3][2] = 0.55396; RxnConstantTable[3][3] = -7.5304; RxnConstantTable[3][4] = -0.014089;
        RxnConstantTable[4][0] = 0.5150;  RxnConstantTable[4][1] = 0.73286;   RxnConstantTable[4][2] = 0.49096; RxnConstantTable[4][3] = -7.5025; RxnConstantTable[4][4] = -0.015938;
        RxnConstantTable[5][0] = 0.50765; RxnConstantTable[5][1] = 0.73575;   RxnConstantTable[5][2] = 0.48042; RxnConstantTable[5][3] = -7.4979; RxnConstantTable[5][4] = -0.016247;

      } else if (iReaction == 15) {

        //N2 + O -> NO + N
        RxnConstantTable[0][0] = 1.3261;  RxnConstantTable[0][1] = 0.75268; RxnConstantTable[0][2] = 1.2474;  RxnConstantTable[0][3] = -4.1857; RxnConstantTable[0][4] = 0.02184;
        RxnConstantTable[1][0] = 1.0653;  RxnConstantTable[1][1] = 0.85417; RxnConstantTable[1][2] = 0.87093; RxnConstantTable[1][3] = -4.0188; RxnConstantTable[1][4] = 0.010721;
        RxnConstantTable[2][0] = 0.96794; RxnConstantTable[2][1] = 0.89131; RxnConstantTable[2][2] = 0.7291;  RxnConstantTable[2][3] = -3.9555; RxnConstantTable[2][4] = 0.006488;
        RxnConstantTable[3][0] = 0.97646; RxnConstantTable[3][1] = 0.89043; RxnConstantTable[3][2] = 0.74572; RxnConstantTable[3][3] = -3.9642; RxnConstantTable[3][4] = 0.007123;
        RxnConstantTable[4][0] = 0.96188; RxnConstantTable[4][1] = 0.89617; RxnConstantTable[4][2] = 0.72479; RxnConstantTable[4][3] = -3.955;  RxnConstantTable[4][4] = 0.006509;
        RxnConstantTable[5][0] = 0.96921; RxnConstantTable[5][1] = 0.89329; RxnConstantTable[5][2] = 0.73531; RxnConstantTable[5][3] = -3.9596; RxnConstantTable[5][4] = 0.006818;

      } else if (iReaction == 16) {

        //NO + O -> O2 + N
        RxnConstantTable[0][0] = 0.35438;   RxnConstantTable[0][1] = -1.8821; RxnConstantTable[0][2] = -0.72111;  RxnConstantTable[0][3] = -1.1797;   RxnConstantTable[0][4] = -0.030831;
        RxnConstantTable[1][0] = 0.093613;  RxnConstantTable[1][1] = -1.7806; RxnConstantTable[1][2] = -1.0975;   RxnConstantTable[1][3] = -1.0128;   RxnConstantTable[1][4] = -0.041949;
        RxnConstantTable[2][0] = -0.003732; RxnConstantTable[2][1] = -1.7434; RxnConstantTable[2][2] = -1.2394;   RxnConstantTable[2][3] = -0.94952;  RxnConstantTable[2][4] = -0.046182;
        RxnConstantTable[3][0] = 0.004815;  RxnConstantTable[3][1] = -1.7443; RxnConstantTable[3][2] = -1.2227;   RxnConstantTable[3][3] = -0.95824;  RxnConstantTable[3][4] = -0.045545;
        RxnConstantTable[4][0] = -0.009758; RxnConstantTable[4][1] = -1.7386; RxnConstantTable[4][2] = -1.2436;   RxnConstantTable[4][3] = -0.949;    RxnConstantTable[4][4] = -0.046159;
        RxnConstantTable[5][0] = -0.002428; RxnConstantTable[5][1] = -1.7415; RxnConstantTable[5][2] = -1.2331;   RxnConstantTable[5][3] = -0.95365;  RxnConstantTable[5][4] = -0.04585;
      }

      break;

    case AIR7:

      if (iReaction <= 6) {

        //N2 + M -> 2N + M
        RxnConstantTable[0][0] = 3.4907;  RxnConstantTable[0][1] = 0.83133; RxnConstantTable[0][2] = 4.0978;  RxnConstantTable[0][3] = -12.728; RxnConstantTable[0][4] = 0.07487;   //n = 1E14
        RxnConstantTable[1][0] = 2.0723;  RxnConstantTable[1][1] = 1.38970; RxnConstantTable[1][2] = 2.0617;  RxnConstantTable[1][3] = -11.828; RxnConstantTable[1][4] = 0.015105;  //n = 1E15
        RxnConstantTable[2][0] = 1.6060;  RxnConstantTable[2][1] = 1.57320; RxnConstantTable[2][2] = 1.3923;  RxnConstantTable[2][3] = -11.533; RxnConstantTable[2][4] = -0.004543; //n = 1E16
        RxnConstantTable[3][0] = 1.5351;  RxnConstantTable[3][1] = 1.60610; RxnConstantTable[3][2] = 1.2993;  RxnConstantTable[3][3] = -11.494; RxnConstantTable[3][4] = -0.00698;  //n = 1E17
        RxnConstantTable[4][0] = 1.4766;  RxnConstantTable[4][1] = 1.62910; RxnConstantTable[4][2] = 1.2153;  RxnConstantTable[4][3] = -11.457; RxnConstantTable[4][4] = -0.00944;  //n = 1E18
        RxnConstantTable[5][0] = 1.4766;  RxnConstantTable[5][1] = 1.62910; RxnConstantTable[5][2] = 1.2153;  RxnConstantTable[5][3] = -11.457; RxnConstantTable[5][4] = -0.00944;  //n = 1E19

      } else if (iReaction > 6 && iReaction <= 13) {

        //O2 + M -> 2O + M
        RxnConstantTable[0][0] = 1.8103;  RxnConstantTable[0][1] = 1.9607;  RxnConstantTable[0][2] = 3.5716;  RxnConstantTable[0][3] = -7.3623;   RxnConstantTable[0][4] = 0.083861;
        RxnConstantTable[1][0] = 0.91354; RxnConstantTable[1][1] = 2.3160;  RxnConstantTable[1][2] = 2.2885;  RxnConstantTable[1][3] = -6.7969;   RxnConstantTable[1][4] = 0.046338;
        RxnConstantTable[2][0] = 0.64183; RxnConstantTable[2][1] = 2.4253;  RxnConstantTable[2][2] = 1.9026;  RxnConstantTable[2][3] = -6.6277;   RxnConstantTable[2][4] = 0.035151;
        RxnConstantTable[3][0] = 0.55388; RxnConstantTable[3][1] = 2.4600;  RxnConstantTable[3][2] = 1.7763;  RxnConstantTable[3][3] = -6.5720;   RxnConstantTable[3][4] = 0.031445;
        RxnConstantTable[4][0] = 0.52455; RxnConstantTable[4][1] = 2.4715;  RxnConstantTable[4][2] = 1.7342;  RxnConstantTable[4][3] = -6.55534;  RxnConstantTable[4][4] = 0.030209;
        RxnConstantTable[5][0] = 0.50989; RxnConstantTable[5][1] = 2.4773;  RxnConstantTable[5][2] = 1.7132;  RxnConstantTable[5][3] = -6.5441;   RxnConstantTable[5][4] = 0.029591;

      } else if (iReaction > 13 && iReaction <= 20) {

        //NO + M -> N + O + M
        RxnConstantTable[0][0] = 2.1649;  RxnConstantTable[0][1] = 0.078577;  RxnConstantTable[0][2] = 2.8508;  RxnConstantTable[0][3] = -8.5422; RxnConstantTable[0][4] = 0.053043;
        RxnConstantTable[1][0] = 1.0072;  RxnConstantTable[1][1] = 0.53545;   RxnConstantTable[1][2] = 1.1911;  RxnConstantTable[1][3] = -7.8098; RxnConstantTable[1][4] = 0.004394;
        RxnConstantTable[2][0] = 0.63817; RxnConstantTable[2][1] = 0.68189;   RxnConstantTable[2][2] = 0.66336; RxnConstantTable[2][3] = -7.5773; RxnConstantTable[2][4] = -0.011025;
        RxnConstantTable[3][0] = 0.55889; RxnConstantTable[3][1] = 0.71558;   RxnConstantTable[3][2] = 0.55396; RxnConstantTable[3][3] = -7.5304; RxnConstantTable[3][4] = -0.014089;
        RxnConstantTable[4][0] = 0.5150;  RxnConstantTable[4][1] = 0.73286;   RxnConstantTable[4][2] = 0.49096; RxnConstantTable[4][3] = -7.5025; RxnConstantTable[4][4] = -0.015938;
        RxnConstantTable[5][0] = 0.50765; RxnConstantTable[5][1] = 0.73575;   RxnConstantTable[5][2] = 0.48042; RxnConstantTable[5][3] = -7.4979; RxnConstantTable[5][4] = -0.016247;

      } else if (iReaction == 21) {

        //N2 + O -> NO + N
        RxnConstantTable[0][0] = 1.3261;  RxnConstantTable[0][1] = 0.75268; RxnConstantTable[0][2] = 1.2474;  RxnConstantTable[0][3] = -4.1857; RxnConstantTable[0][4] = 0.02184;
        RxnConstantTable[1][0] = 1.0653;  RxnConstantTable[1][1] = 0.85417; RxnConstantTable[1][2] = 0.87093; RxnConstantTable[1][3] = -4.0188; RxnConstantTable[1][4] = 0.010721;
        RxnConstantTable[2][0] = 0.96794; RxnConstantTable[2][1] = 0.89131; RxnConstantTable[2][2] = 0.7291;  RxnConstantTable[2][3] = -3.9555; RxnConstantTable[2][4] = 0.006488;
        RxnConstantTable[3][0] = 0.97646; RxnConstantTable[3][1] = 0.89043; RxnConstantTable[3][2] = 0.74572; RxnConstantTable[3][3] = -3.9642; RxnConstantTable[3][4] = 0.007123;
        RxnConstantTable[4][0] = 0.96188; RxnConstantTable[4][1] = 0.89617; RxnConstantTable[4][2] = 0.72479; RxnConstantTable[4][3] = -3.955;  RxnConstantTable[4][4] = 0.006509;
        RxnConstantTable[5][0] = 0.96921; RxnConstantTable[5][1] = 0.89329; RxnConstantTable[5][2] = 0.73531; RxnConstantTable[5][3] = -3.9596; RxnConstantTable[5][4] = 0.006818;

      } else if (iReaction == 22) {

        //NO + O -> O2 + N
        RxnConstantTable[0][0] = 0.35438;   RxnConstantTable[0][1] = -1.8821; RxnConstantTable[0][2] = -0.72111;  RxnConstantTable[0][3] = -1.1797;   RxnConstantTable[0][4] = -0.030831;
        RxnConstantTable[1][0] = 0.093613;  RxnConstantTable[1][1] = -1.7806; RxnConstantTable[1][2] = -1.0975;   RxnConstantTable[1][3] = -1.0128;   RxnConstantTable[1][4] = -0.041949;
        RxnConstantTable[2][0] = -0.003732; RxnConstantTable[2][1] = -1.7434; RxnConstantTable[2][2] = -1.2394;   RxnConstantTable[2][3] = -0.94952;  RxnConstantTable[2][4] = -0.046182;
        RxnConstantTable[3][0] = 0.004815;  RxnConstantTable[3][1] = -1.7443; RxnConstantTable[3][2] = -1.2227;   RxnConstantTable[3][3] = -0.95824;  RxnConstantTable[3][4] = -0.045545;
        RxnConstantTable[4][0] = -0.009758; RxnConstantTable[4][1] = -1.7386; RxnConstantTable[4][2] = -1.2436;   RxnConstantTable[4][3] = -0.949;    RxnConstantTable[4][4] = -0.046159;
        RxnConstantTable[5][0] = -0.002428; RxnConstantTable[5][1] = -1.7415; RxnConstantTable[5][2] = -1.2331;   RxnConstantTable[5][3] = -0.95365;  RxnConstantTable[5][4] = -0.04585;

      } else if (iReaction == 23) {

        //N + O -> NO+ + e-
        RxnConstantTable[0][0] = -2.1852;   RxnConstantTable[0][1] = -6.6709; RxnConstantTable[0][2] = -4.2968; RxnConstantTable[0][3] = -2.2175; RxnConstantTable[0][4] = -0.050748;
        RxnConstantTable[1][0] = -1.0276;   RxnConstantTable[1][1] = -7.1278; RxnConstantTable[1][2] = -2.637;  RxnConstantTable[1][3] = -2.95;   RxnConstantTable[1][4] = -0.0021;
        RxnConstantTable[2][0] = -0.65871;  RxnConstantTable[2][1] = -7.2742; RxnConstantTable[2][2] = -2.1096; RxnConstantTable[2][3] = -3.1823; RxnConstantTable[2][4] = 0.01331;
        RxnConstantTable[3][0] = -0.57924;  RxnConstantTable[3][1] = -7.3079; RxnConstantTable[3][2] = -1.9999; RxnConstantTable[3][3] = -3.2294; RxnConstantTable[3][4] = 0.016382;
        RxnConstantTable[4][0] = -0.53538;  RxnConstantTable[4][1] = -7.3252; RxnConstantTable[4][2] = -1.937;  RxnConstantTable[4][3] = -3.2572; RxnConstantTable[4][4] = 0.01823;
        RxnConstantTable[5][0] = -0.52801;  RxnConstantTable[5][1] = -7.3281; RxnConstantTable[5][2] = -1.9264; RxnConstantTable[5][3] = -3.2618; RxnConstantTable[5][4] = 0.01854;
      }
      break;
    }
  }

  void Calc_DiffusionCoeff_WBE(su2double *V, su2double *Xs);

  void Calc_Viscosity_WBE(su2double *V, su2double *Xs);

  void Calc_ThermalConductivity_WBE( su2double *V, su2double *Xs);

  void Calc_DiffusionCoeff_GY(su2double *V);

  void Calc_Viscosity_GY(su2double *V);

  void Calc_ThermalConductivity_GY(su2double *V);

};

/*!
 * \derived class CTNE2Gas
 * \brief Child class for defining the two-temperature gas model.
 * \author: W. Maier
 */
class CTNE2Gas : public CMultiSpeciesGas {

protected:

private:

public:

  /*!
   * \brief Constructor of the class.
   */
  CTNE2Gas(unsigned short val_nSpecies, unsigned short val_nDim, bool val_ionization, bool val_viscous);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CTNE2Gas(void);

  void InitializeMixture(CConfig *Config);

  su2double Calc_CvTraRotSpecies(su2double *Ms, su2double Ru, unsigned short val_Species);

  su2double Calc_CvVibElSpecies(su2double val_Tve, unsigned short val_Species);

  su2double* Calc_CpTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  su2double* Calc_CpVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  su2double  Calc_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Calc_MixtureEnergy(su2double* cs,
                                su2double sqvel,
                                su2double rho,
                                su2double T, su2double Tve);

  su2double Calc_MixtureEnergies(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double Calc_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve);

  su2double Calc_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double Calc_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Calc_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double Calc_Enthalpies(su2double val_T, su2double val_eves, unsigned short val_Species);

  su2double* Calc_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double  Calc_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double Calc_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve);

  su2double Calc_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve);

  su2double Calc_SoundSpeed(su2double *cs, su2double rhoCvtr, su2double rho, su2double P);

  su2double Calc_Density(su2double *MassFrac,
                         su2double T, su2double Tve, su2double P);

  void CalcdPdU(su2double *V, su2double *val_eves, su2double *val_dPdU);

  su2double CalcEve(su2double val_Tve,
                    unsigned short val_Species);

  void CalcdTdU(su2double *V, su2double *val_dTdU);

  void CalcdTvedU(su2double *V, su2double *val_eves, su2double *val_dTvedU);

  void Calc_TransportCoeff(CConfig *config, su2double *V);

  void Calc_DiffusionCoeff_WBE(su2double *V, su2double *Xs);

  void Calc_Viscosity_WBE(su2double *V, su2double *Xs );

  void Calc_ThermalConductivity_WBE(su2double *V, su2double *Xs);

  void Calc_DiffusionCoeff_GY(su2double *V);

  void Calc_Viscosity_GY(su2double *V);

  void Calc_ThermalConductivity_GY(su2double *V);

};

/*!
 * \derived class CMutationGas
 * \brief Child class for defining the Mutation++ Link.
 * \author: C. Garbacz
 */
class CMutationGas : public CMultiSpeciesGas {

protected:

  string MixtureFile; /*!< \brief Presence of charged species in gas mixture. */
  CMutation *mutation; /*!< \brief Presence of charged species in gas mixture. */
  su2double *comp, *Cp_ks,  *Cp_trs, *Cp_ves, E, gamma, gammaFrozen, gammaEquilibrium, *Ws, Tref, *hs; /*!< \brief Presence of charged species in gas mixture. */
  su2double  mu, *Ds; /*!< \brief Presence of charged species in gas mixture. */
  su2double  a, *Xs; /*!< \brief Presence of charged species in gas mixture. */
  vector<su2double> Ms, Cv_trs, Cv_ks, Energies, OmegaVT, hf, Energies_Species,  Cv_ves, Temp, lambda; /*!< \brief Presence of charged species in gas mixture. */
  //Mutation::MixtureOptions opt;, THIS DONT WORK, ASK CATARINA, DELETE ME

private:

public:

  /*!
   * \brief Constructor of the class.
   */
  CMutationGas(string MixFile, string Transport, unsigned short val_nSpecies,
               unsigned short val_nDim, bool val_ionization, bool val_viscous);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CMutationGas(void);

  /*!
   * \brief Initial gas mizture for mutation++
   * \param[in] config - problem setup information.
   */
  void InitializeMixture(CConfig *config);

  /*!
   * \brief Get flow ionization bool
   */
  bool Calc_Ionization();

  /*!
   * \brief Calculate molar mass of gas/species
   */
  vector<su2double> Calc_MolarMass();

  /*!
   * \brief Calculate the Cv of the T-R mode
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  vector<su2double> Calc_CvTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the Cv of the V-E mode
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  vector<su2double> Calc_CvVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the Cp of the T-R mode
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double* Calc_CpTraRotSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the Cp of the T-R mode
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double* Calc_CpVibElSpecies(su2double* cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate gamma of gas
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double  Calc_Gamma(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the frozen gamma
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double  Calc_GammaFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the equilibrium gamma.
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double  Calc_GammaEquilibrium(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the total energy of gas
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double  Calc_MixtureEnergy(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the t-r and v-e energies.
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  vector<su2double> Calc_MixtureEnergies(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the species energies
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  vector<su2double> Calc_SpeciesEnergies(su2double* cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the net production rates of species
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double* Calc_NetProductionRates(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the vibrational energy source term
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  vector<su2double> Calc_VTEnergysourceTerm(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the reference temperature of the gas.
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double  Calc_ReferenceTemperature(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the formation enthalpies of the species.
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  vector<su2double> Calc_EnthalpiesFormation(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the enthalpy
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double* Calc_Enthalpies(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the diffusion coefficients.
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double* Calc_DiffusionCoeff(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the viscosity
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double  Calc_Viscosity(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the thermal conductivity.
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  vector<su2double> Calc_ThermalConductivity(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the temperatures
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] rhoE - Total energy.
   * \param[in] rhoEve - Vibe-Elec energy.
   */
  vector<su2double> Calc_Temperatures(su2double *cs, su2double rho, su2double rhoE, su2double rhoEve);

  /*!
   * \brief Calculate the frozen sound speed.
   * \param[in] cs - Mass fractions of species.
   * \param[in] rho - Gas total density.
   * \param[in] T - Temperature.
   * \param[in] Tve - Vibe-Elec Temperature.
   */
  su2double  Calc_SoundSpeedFrozen(su2double *cs, su2double rho, su2double T, su2double Tve);

  /*!
   * \brief Calculate the total density.
   * \param[in] T - Temperature.
   * \param[in] Xs - Rotational modes.
   * \param[in] P - Pressure.
   */
  su2double  Calc_Density(su2double T, su2double *Xs, su2double P);

};

/*!
 * \class CConstantDensity
 * \brief Child class for defining a constant density gas model (incompressible only).
 * \author: T. Economon
 */
class CConstantDensity : public CFluidModel {

protected:

public:

  /*!
   * \brief Constructor of the class.
   */
  CConstantDensity(void);

  /*!
   * \brief Constructor of the class.
   */
  CConstantDensity(su2double val_Density, su2double val_Cp);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CConstantDensity(void);

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] T - Temperature value at the point.
   */
  void SetTDState_T(su2double val_Temperature);

};

/*!
 * \class CIncIdealGas
 * \brief Child class for defining an incompressible ideal gas model.
 * \author: T. Economon
 */
class CIncIdealGas : public CFluidModel {

protected:
  su2double Gas_Constant,  /*!< \brief Gas Constant. */
  Gamma;                   /*!< \brief Heat Capacity Ratio. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CIncIdealGas(void);

  /*!
   * \brief Constructor of the class.
   */
  CIncIdealGas(su2double val_Cp, su2double val_gas_constant, su2double val_operating_pressure);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CIncIdealGas(void);

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] T - Temperature value at the point.
   */

  void SetTDState_T(su2double val_Temperature);

};

/*!
 * \class CIncIdealGasPolynomial
 * \brief Child class for defining a custom incompressible ideal gas model.
 * \author: T. Economon
 */
class CIncIdealGasPolynomial : public CFluidModel {

protected:
  unsigned short nPolyCoeffs; /*!< \brief Number of coefficients in the temperature polynomial. */
  su2double Gas_Constant,     /*!< \brief Specific Gas Constant. */
  *b,                         /*!< \brief Polynomial coefficients for Cp as a function of temperature. */
  Gamma;                      /*!< \brief Ratio of specific heats. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CIncIdealGasPolynomial(void);

  /*!
   * \brief Constructor of the class.
   */
  CIncIdealGasPolynomial(su2double val_gas_constant, su2double val_operating_pressure);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CIncIdealGasPolynomial(void);

  /*!
   * \brief Set the temperature polynomial coefficients for variable Cp.
   * \param[in] config - configuration container for the problem.
   */
  void SetCpModel(CConfig *config);

  /*!
   * \brief Set the Dimensionless State using Temperature.
   * \param[in] T - Temperature value at the point.
   */
  void SetTDState_T(su2double val_temperature);

};

#include "fluid_model.inl"
