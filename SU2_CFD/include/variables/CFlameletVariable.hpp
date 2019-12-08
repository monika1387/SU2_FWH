/*!
 * \file CFlameletVariable.hpp
 * \brief Declaration of the variables of the flamelet model.
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

#include "CScalarVariable.hpp"

/*!
 * \class CFlameletVariable
 * \brief Main class for defining the variables of the flamelet model.
 * \author D. Mayer, T. Economon, N.Beishuizen, G. Tamanampudi
 */
class CFlameletVariable final : public CScalarVariable {
protected:
  MatrixType source_prog;            /*!< \brief Vector of the source terms from the lookup table for each sclar equation */
  MatrixType scalar_table_paraview;  /*!< \brief Vector of the values of the scalar from the lookup table for visualization */
  
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_scalar_inf - Scalar variable value (initialization value).
   * \param[in] npoint         - Number of points/nodes/vertices in the domain.
   * \param[in] ndim           - Number of dimensions of the problem.
   * \param[in] nvar           - Number of variables of the problem.
   * \param[in] config         - Definition of the particular problem.
   */
  CFlameletVariable(su2double     *val_scalar_inf,
                    unsigned long npoint,
                    unsigned long ndim,
                    unsigned long nvar,
                    CConfig       *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CFlameletVariable() = default;
  
  /*!
   * \brief Set the value of the progress variable source term
   * \param[in] val_- the .
   * \param[in] val_ivar        - eqn. index to the .
   */
  inline void SetSourceProg(unsigned long  iPoint,
                            su2double      val_source_prog,
                            unsigned short val_ivar) override {
    source_prog(iPoint,val_ivar) = val_source_prog;
  }
  
  /*!
   * \brief Get the value of the progress variable source term
   * \param[in] val_ivar - eqn. index to the progress variable source term
   * \return Value of the progress variable source term
   */
  inline su2double GetSourceProg(unsigned long  iPoint,
                                 unsigned short val_ivar) override {
    return source_prog(iPoint,val_ivar);
  }
  
  /*!
   * \brief Get the value of the progress variable source term
   * \return Pointer to the progress variable source term
   */
  inline su2double *GetSourceProg(unsigned long iPoint) override {
    return source_prog[iPoint];
  }
  
  /*!
   * \brief Store the scalar variables interpolated from the flamelet table to output in paraview
   * \param[out] stores the interpolated scalar variables in the array scalar_table_paraview
   */
  inline void SetScalarTableParaview(unsigned long  iPoint,
                                     su2double      *val_scalar_table,
                                     unsigned short nvar) override {
    for (int ivar=0;ivar<nvar;ivar++)
    scalar_table_paraview(iPoint,ivar) = val_scalar_table[ivar];
  }
  
  /*!
   * \brief Retrieve the stored scalar variables interpolated from the table to output in paraview
   * \param[out] returns the scalar variables  stored in the array sclar_table_paraview
   */
  inline su2double GetScalarTableParaview(unsigned long  iPoint,unsigned short val_ivar) override {
    return scalar_table_paraview(iPoint,val_ivar);
  }
  
};
