/*!
 * \file CFlameletSolver.cpp
 * \brief Main subroutines for the flamelet model solver.
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

#include "../../include/solvers/CFlameletSolver.hpp"
#include "../../include/variables/CFlameletVariable.hpp"
#include "../../include/CFluidFlamelet.hpp"

CFlameletSolver::CFlameletSolver(void) : CScalarSolver() {
  
  Inlet_ScalarVars = NULL;
  
}

CFlameletSolver::CFlameletSolver(CGeometry *geometry,
                                                 CConfig *config,
                                                 unsigned short iMesh)
: CScalarSolver(geometry, config) {
  
  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  su2double Density_Inf, Viscosity_Inf;
  
  bool turbulent = ((config->GetKind_Solver() == RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));
  bool turb_SST  = ((turbulent) && (config->GetKind_Turb_Model() == SST));
  bool turb_SA   = ((turbulent) && (config->GetKind_Turb_Model() == SA));
  vector <string> variable_names;
  
  /*--- Dimension of the problem --> passive scalar will only ever
   have a single equation. Other child classes of CScalarSolver
   will have variable numbers of equations. ---*/
  
  nVar     = 2;
  nPrimVar = 2;
  variable_names.push_back("progress_variable");
  variable_names.push_back("enthalpy");
  
    /*--- Specify the names of variables (scalars) in the same order as they are being solved --- */
  /*--- Please set SCALAR_VARIABLES in option_structure.hpp accordingly */
  
  SetScalarNames(variable_names);
  
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim = geometry->GetnDim();
  
  /*--- Fluid model pointer initialization ---*/
  
  FluidModel = NULL;
  
  /*--- Define some auxiliar vector related with the residual ---*/
  
  Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual    [iVar] = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar] = 0.0;
  Residual_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i  [iVar] = 0.0;
  Residual_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j  [iVar] = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar] = 0.0;
  Res_Conv     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv    [iVar] = 0.0;
  Res_Visc     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc    [iVar] = 0.0;
  
  
  /*--- Define some structures for locating max residuals ---*/
  
  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  /*--- Define some auxiliar vector related with the solution ---*/
  
  Solution = new su2double[nVar];
  Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];
  
  /*--- Define some auxiliar vector related with the geometry ---*/
  
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
  
  /*--- Define some auxiliar vector related with the flow solution ---*/
  
  FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Initialization of the structure of the whole Jacobian ---*/
  
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Combustion Scalar). MG level: " << iMesh <<"." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  
  if (config->GetKind_Linear_Solver_Prec() == LINELET) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Computation of gradients by least squares ---*/
  
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*transpose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    
    /*--- c vector := transpose(WA)*(Wb) ---*/
    
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Initialize lower and upper limits---*/
  
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];
  
  lowerlimit[0] = -1.0e15;
  upperlimit[0] =  1.0e15;
  
  /*--- Read farfield conditions from config ---*/
  
  Density_Inf   = config->GetDensity_FreeStreamND();
  Viscosity_Inf = config->GetViscosity_FreeStreamND();
  
  /*--- Set up fluid model for the diffusivity ---*/
  
  su2double Diffusivity_Ref = 1.0;
  
  su2double DiffusivityND = config->GetDiffusivity_Constant()/Diffusivity_Ref;
  
  config->SetDiffusivity_ConstantND(DiffusivityND);
  
  /*--- Scalar variable state at the far-field. ---*/
  
  Scalar_Inf = new su2double[nVar];
  
  //FIXME daniel: scalar_inf should be set depending on inlet temperature
  //              can be done in the initial condition, also set Inlet_Scalar_vars!
  for (iVar = 0; iVar < nVar; iVar++)
    Scalar_Inf[iVar] = config->GetScalar_Init()[iVar];
  
  /*--- Initialize the solution to the far-field state everywhere. ---*/
  
  nodes = new CFlameletVariable(Scalar_Inf, nPoint, nDim, nVar, config);
  SetBaseClassPointerToNodes();

  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Initialize quantities for SlidingMesh Interface ---*/
  
  unsigned long iMarker;
  
  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++){
    
    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;
    
    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){
      
      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];
      
      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];
        
        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }
      
    }
  }
  
  /*-- Allocation of inlets has to happen in derived classes
   (not CScalarSolver), due to arbitrary number of scalar variables.
   First, we also set the column index for any inlet profiles. ---*/
  
  Inlet_Position = nDim*2+2;
  if (turbulent) {
    if (turb_SA) Inlet_Position += 1;
    else if (turb_SST) Inlet_Position += 2;
  }
  
  Inlet_ScalarVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    Inlet_ScalarVars[iMarker] = new su2double*[nVertex[iMarker]];
    for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
      Inlet_ScalarVars[iMarker][iVertex] = new su2double[nVar];
      for (unsigned short iVar = 0; iVar < nVar; iVar++)
        Inlet_ScalarVars[iMarker][iVertex][0] = Scalar_Inf[iVar];
    }
  }
  
}

CFlameletSolver::~CFlameletSolver(void) {
  
  unsigned long iMarker, iVertex;
  unsigned short iVar;
  
  if ( SlidingState != NULL ) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != NULL ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }
  
  if ( SlidingStateNodes != NULL ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      if (SlidingStateNodes[iMarker] != NULL)
        delete [] SlidingStateNodes[iMarker];
    }
    delete [] SlidingStateNodes;
  }
  
}

void CFlameletSolver::Preprocessing(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CConfig *config,
                                            unsigned short iMesh,
                                            unsigned short iRKStep,
                                            unsigned short RunTime_EqSystem,
                                            bool Output) {
  
  unsigned long ErrorCounter = 0;
  unsigned long InnerIter = config->GetInnerIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) &&
                           (InnerIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_scalar   = ((config->GetKind_SlopeLimit_Scalar() != NO_LIMITER) &&
                           (InnerIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  
  /*--- Set the primitive variables ---*/
  
  ErrorCounter = SetPrimitive_Variables(solver_container, config, Output);
  
  /*--- Initialize the Jacobian matrices ---*/
  
  if (!disc_adjoint) Jacobian.SetValZero();
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Upwind second order reconstruction ---*/
  
  if (limiter_scalar) SetSolution_Limiter(geometry, config);
  
  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);
  
}

void CFlameletSolver::Postprocessing(CGeometry *geometry,
                                             CSolver **solver_container,
                                             CConfig *config,
                                             unsigned short iMesh) { }

unsigned long CFlameletSolver::SetPrimitive_Variables(CSolver **solver_container,
                                                              CConfig *config,
                                                              bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double Density, Temperature, lam_visc = 0.0, eddy_visc = 0.0, Cp;
  unsigned short turb_model = config->GetKind_Turb_Model();
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Retrieve the density, temperature, Cp, and laminar viscosity. ---*/
    
    Density     = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    Cp          = solver_container[FLOW_SOL]->GetNodes()->GetSpecificHeatCp(iPoint);
    Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    lam_visc    = solver_container[FLOW_SOL]->GetNodes()->GetLaminarViscosity(iPoint);
    
    /*--- Retrieve the value of the kinetic energy (if needed) ---*/
    
    if (turb_model != NONE) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
    }
    
    /*--- Compute and store the mass diffusivity. ---*/
    
    CFluidModel * fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
    su2double * scalars = nodes->GetSolution(iPoint);
    fluid_model_local->SetTDState_T(Temperature,scalars);
    
    nodes->SetDiffusivity(iPoint,
                                 fluid_model_local->GetMassDiffusivity(), I_PROG_VAR);
    
    nodes->SetDiffusivity(iPoint,
                                 fluid_model_local->GetMassDiffusivity(), I_ENTHALPY);
    
    /*--- Initialize the convective, source and viscous residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  return ErrorCounter;
  
}

void CFlameletSolver::SetInitialCondition(CGeometry **geometry,
                                                  CSolver ***solver_container,
                                                  CConfig *config,
                                                  unsigned long ExtIter) {
  su2double *coords;
  bool Restart   = (config->GetRestart() || config->GetRestart_Flow());
  
  
  if ((!Restart) && ExtIter == 0) {
    if (rank == MASTER_NODE)
      cout << "Initializing progress variable and temperature (initial condition)." << endl;
    
    su2double *scalar_init    = new su2double[nVar];
    su2double temp_init;
    su2double prog_unburnt    = 0.0;
    su2double prog_burnt      = 2.09379237374982e-01;
    su2double flame_offset    = config->GetFlameOffset();
    su2double flame_thickness = config->GetFlameThickness();
    su2double burnt_thickness = config->GetBurntThickness();

    su2double temp_inlet = 300.;
    su2double prog_inlet = 0.0;
    
    CFluidModel *fluid_model_local;
    
    for (unsigned long i_mesh = 0; i_mesh <= config->GetnMGLevels(); i_mesh++) {
      for (unsigned long i_point = 0; i_point < geometry[i_mesh]->GetnPoint(); i_point++) {
        
        for (unsigned long i_var = 0; i_var < nVar; i_var++)
          Solution[i_var] = 0.0;
        
        // FIXME dan: This works only in x direction and it should be possible to specify a normal for the flame front
        coords = geometry[i_mesh]->node[i_point]->GetCoord();

        if (coords[0] < flame_offset){ /* unburnt region */

          scalar_init[I_PROG_VAR] = prog_unburnt;

        } else if (coords[0] < flame_offset + flame_thickness){ /* flame zone */

          scalar_init[I_PROG_VAR] = prog_unburnt + (prog_burnt - prog_unburnt) / flame_thickness * (coords[0] - flame_offset);

        } else if (coords[0] < flame_offset + flame_thickness + burnt_thickness){ /* burnt region */

          scalar_init[I_PROG_VAR] = prog_burnt;

        } else {

          scalar_init[I_PROG_VAR] = prog_unburnt;

        }
        

        fluid_model_local = solver_container[i_mesh][FLOW_SOL]->GetFluidModel();
        su2double enth_inlet = fluid_model_local->GetEnthFromTemp(prog_inlet,temp_inlet);
        
        scalar_init[I_ENTHALPY] = enth_inlet;
        
        solver_container[i_mesh][SCALAR_SOL]->GetNodes()->SetSolution(i_point, scalar_init);

      }
      
      solver_container[i_mesh][SCALAR_SOL]->InitiateComms(geometry[i_mesh], config, SOLUTION);
      solver_container[i_mesh][SCALAR_SOL]->CompleteComms(geometry[i_mesh], config, SOLUTION);
      
      solver_container[i_mesh][FLOW_SOL]->InitiateComms(geometry[i_mesh], config, SOLUTION);
      solver_container[i_mesh][FLOW_SOL]->CompleteComms(geometry[i_mesh], config, SOLUTION);
      
      solver_container[i_mesh][FLOW_SOL]->Preprocessing( geometry[i_mesh], solver_container[i_mesh], config, i_mesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
      
    }
    delete[] scalar_init;
  }
}

void CFlameletSolver::SetPreconditioner(CGeometry *geometry, CSolver **solver_container,  CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  
  su2double  BetaInc2, Density, dRhodT, dRhodC, Temperature, Cp, Delta;
  
  bool variable_density = (config->GetKind_DensityModel() == VARIABLE);
  bool implicit         = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Access the primitive variables at this node. ---*/
    
    Density     = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
    BetaInc2    = solver_container[FLOW_SOL]->GetNodes()->GetBetaInc2(iPoint);
    Cp          = solver_container[FLOW_SOL]->GetNodes()->GetSpecificHeatCp(iPoint);
    Temperature = solver_container[FLOW_SOL]->GetNodes()->GetTemperature(iPoint);
    
    /*--- We need the derivative of the equation of state to build the
     preconditioning matrix. For now, the only option is the ideal gas
     law, but in the future, dRhodT should be in the fluid model. ---*/
    
    if (variable_density) {
      dRhodT = -Density/Temperature;
    } else {
      dRhodT = 0.0;
    }
    
    /*--- Passive scalars have no impact on the density. ---*/
    
    dRhodC = 0.0;
    
    /*--- Modify matrix diagonal with term including volume and time step. ---*/
    
    su2double Vol = geometry->node[iPoint]->GetVolume();
    Delta = Vol / (config->GetCFLRedCoeff_Scalar()*
                   solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));
    
    /*--- Calculating the inverse of the preconditioning matrix
     that multiplies the time derivative during time integration. ---*/
    
    if (implicit) {
      
      for (int i_var = 0; i_var < nVar; i_var++) {
        for (int j_var = 0; j_var < nVar; j_var++) {
          Jacobian_i[i_var][j_var] = 0.0;
        }
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        total_index = iPoint*nVar+iVar;
        
        su2double scalar = nodes->GetSolution(iPoint, iVar);
        
        /*--- Add the extra Jacobian term to the scalar system. ---*/
        
        su2double Jaccomp = scalar * dRhodC + Density;
        
        Jacobian_i[iVar][iVar] = Jaccomp*Delta;
        
      }
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
    
  }
  
}

void CFlameletSolver::Source_Residual(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CNumerics *numerics,
                                              CNumerics *second_numerics,
                                              CConfig *config,
                                              unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  
  bool implicit       = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool axisymmetric   = config->GetAxisymmetric();
  bool viscous        = config->GetViscous();
  
  /*--- Combustion source term implementation. (active by default)  ---*/
  
  su2double  temperature_dummy = 73                 ;
  su2double* scalars           = new su2double[nVar];
  su2double* source            = new su2double[nVar];
  pair<su2double, su2double> limits_lut_prog;
  
  CFluidModel *fluid_model_local;
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Get volume of the dual cell. ---*/
    
    su2double Volume = geometry->node[iPoint]->GetVolume();
    
    /*--- Compute the production term. ---*/
    for (iVar = 0; iVar < nVar; iVar ++)
      scalars[iVar] = nodes->GetSolution(iPoint, iVar);
    
    fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
    
    /*--- clip progress variable to table limits ---*/
    
    limits_lut_prog  = fluid_model_local->GetTableLimitsProg();
    
    // FIXME TDE: we can disable this lookup to save computational cost
    // since the source should already be saved in the fluid model.
    fluid_model_local-> SetTDState_T(temperature_dummy, scalars) ;
    
    source [ I_PROG_VAR ] = fluid_model_local->GetSourceProg()      ;
    source [ I_ENTHALPY ] = 0.0                                     ;

    
    for (iVar = 0; iVar < nVar; iVar ++)
      Residual[iVar] = source[iVar] * Volume;
    
    /*---  set progress variable source for paraview viewing ---*/
    // nijso FIXME: possibility to write all source terms to file
    // FIXME daniel: this should be changed somehow. At themoment we store all sources in source prog which is confusing
    // The following 5 lines  should only be called in an iteration in which the VTK files are dumped.
    nodes->SetSourceProg(iPoint, source[I_PROG_VAR], I_PROG_VAR);
    nodes->SetSourceProg(iPoint, source[I_ENTHALPY], I_ENTHALPY);
    
    if (config->GetnFlameletTableOutput() > 0){
    fluid_model_local-> SetFluidFlameletTableOutput(config, temperature_dummy, scalars) ;
    nodes->SetScalarTableParaview(iPoint, fluid_model_local->GetScalarTable(),config->GetnFlameletTableOutput()); }   
    
    /*--- Implicit part for production term (to do). ---*/
    
    for (int i_var = 0; i_var < nVar; i_var++) {
      for (int j_var = 0; j_var < nVar; j_var++) {
        Jacobian_i[i_var][j_var] = 0.0;
      }
    }
    Jacobian_i[0][0] = Volume*fluid_model_local->GetdSourcePVdPV();
    
    
    /*--- Add Residual ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    
    /*--- Implicit part ---*/
    
    if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    
  }
  
  /*--- Axisymmetry source term for the scalar equation. ---*/
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Primitive variables w/o reconstruction ---*/
      
      second_numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint), NULL);
      
      /*--- Scalar variables w/o reconstruction ---*/
      
      second_numerics->SetScalarVar(nodes->GetSolution(iPoint), NULL);
      
      /*--- Mass diffusivity coefficients. ---*/
      
      second_numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint),
                                         NULL);
      
      /*--- Set control volume ---*/
      
      second_numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set y coordinate ---*/
      
      second_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                                NULL);
      
      /*--- If viscous, we need gradients for extra terms. ---*/
      
      if (viscous) {
        
        /*--- Gradient of the scalar variables ---*/
        
        second_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint), NULL);
        
      }
      
      /*--- Compute Source term Residual ---*/
      
      second_numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
}

void CFlameletSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                               CNumerics *conv_numerics,
                               CNumerics *visc_numerics, CConfig *config,
                               unsigned short val_marker) {

  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint, total_index;
  su2double *V_inlet, *V_domain, *Normal, *inlet_scalar;
  su2double *Coords;
  
  Normal = new su2double[nDim];
  inlet_scalar = new su2double[nVar];
  
  bool        grid_movement = config->GetGrid_Movement      (          );
  string      Marker_Tag    = config->GetMarker_All_TagBound(val_marker);
  su2double   temp_inlet    = config->GetInlet_Ttotal       (Marker_Tag);
              inlet_scalar  = config->GetInlet_ScalarVal    (Marker_Tag);
 
  CFluidModel  *fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();
 
  inlet_scalar[I_ENTHALPY] = fluid_model_local->GetEnthFromTemp(inlet_scalar[I_PROG_VAR], temp_inlet);
 
  /*--- Loop over all the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

   /* Dirichlet boundary condition at the inlet for scalars */

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/

       if (geometry->node[iPoint]->GetDomain()) {
    
        nodes->SetSolution_Old(iPoint, inlet_scalar);

        for (iDim = 0; iDim < nVar; iDim++) {
          LinSysRes.SetBlock_Zero(iPoint, iDim);
          nodes->SetVal_ResTruncError_Zero(iPoint, iDim);
        }

        /*--- Includes 1 in the diagonal ---*/
         for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }      
    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;
}

void CFlameletSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config,
                                unsigned short val_marker) {

  unsigned long iPoint, iVertex, Point_Normal, total_index;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  CFluidModel *fluid_model_local = solver_container[FLOW_SOL]->GetFluidModel();  
  Normal = new su2double[nDim];
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    /* strong zero flux Neumann boundary condition at the outlet */

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the outlet ---*/
        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor(); 
          
        for (iVar = 0; iVar < nVar; iVar++)
          Solution[iVar] = nodes->GetSolution(Point_Normal, iVar);
        nodes->SetSolution_Old(iPoint, Solution);
    
        for (iDim = 0; iDim < nVar; iDim++){
          LinSysRes.SetBlock_Zero(iPoint, iDim);
          nodes->SetVal_ResTruncError_Zero(iPoint, iDim);
        }

        /*--- Includes 1 in the diagonal ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
    }
  }
  delete [] Normal;
}

void CFlameletSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CNumerics *conv_numerics,
                                              CNumerics *visc_numerics,
                                              CConfig *config,
                                              unsigned short val_marker) {
  unsigned short iVar, jVar, iDim;
  unsigned long iVertex, iPoint, total_index;

  bool implicit                             = config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT;
  string Marker_Tag                         = config->GetMarker_All_TagBound(val_marker);
  su2double temp_wall                       = config->GetIsothermal_Temperature(Marker_Tag);
  CFluidModel *fluid_model_local            = solver_container[FLOW_SOL]->GetFluidModel();    
  su2double enth_wall, prog_wall;
  

  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv[iVar] = 0.0;
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Set enthalpy on the wall ---*/
            
      prog_wall = solver_container[SCALAR_SOL]->GetNodes()->GetSolution(iPoint)[I_PROG_VAR];
      enth_wall = fluid_model_local->GetEnthFromTemp(prog_wall,temp_wall);      

      /*--- Impose the value of the enthalpy as a strong boundary
       condition (Dirichlet) and remove any 
       contribution to the residual at this node. ---*/

      nodes->SetSolution_Old(iPoint, I_ENTHALPY, enth_wall);  
          

      LinSysRes.SetBlock_Zero(iPoint, I_ENTHALPY);
      nodes->SetVal_ResTruncError_Zero(iPoint, I_ENTHALPY);

      if (implicit) {

        total_index = iPoint*nVar+I_ENTHALPY;

        Jacobian.DeleteValsRowi(total_index);
      }    
    }

  }
  
}