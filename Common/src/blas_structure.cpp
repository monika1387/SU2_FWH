/*!
 * \file blas_structure.cpp
 * \brief Implementation of the functions that either simulate BLAS functionality
          or interface to an actual BLAS implementation.
 * \author E. van der Weide
 * \version 6.0.1 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/blas_structure.hpp"
#include <Eigen/Core>

using namespace Eigen;

/* Constructor. */
CBlasStructure::CBlasStructure(void) {}

/* Destructor. Nothing to be done. */
CBlasStructure::~CBlasStructure(void) {}

/* Dense matrix multiplication, gemm functionality. */
void CBlasStructure::gemm(const int M,        const int N,        const int K,
                          const su2double *A, const su2double *B, su2double *C,
                          CConfig *config) {

  /* Initialize the variable for the timing, if profiling is active. */
#ifdef PROFILE
  double timeGemm;
  if( config ) config->GEMM_Tick(&timeGemm);
#endif

  /* The matrices are stored in row major order, but the actual multiplication
     is carried out in column major order. Therefore m and n are swapped as well
     as A and B. I.e. the end result is C = B*A, where the dimensions of
     the matrices are: A(k,n), B(m,k), C(m,n).
     First map the the given arrays to the matrix structure. */
  const Map<const Matrix<su2double,Dynamic,Dynamic> > matA(A, K, M);
  const Map<const Matrix<su2double,Dynamic,Dynamic> > matB(B, N, K);

  Map<Matrix<su2double,Dynamic,Dynamic> > matC(C, N, M);

  /* Carry out the matrix multiplication. */
  matC.noalias() = matB * matA;

  /* Store the profiling information, if needed. */
#ifdef PROFILE
  if( config ) config->GEMM_Tock(timeGemm, M, N, K);
#endif
}

/* Dense matrix vector multiplication, gemv functionality. */
void CBlasStructure::gemv(const int M,        const int N,   const su2double *A,
                          const su2double *x, su2double *y) {

  /* Native implementation of the matix vector product.
     Initialize the elements of y to zero. */
  memset(y, 0, M*sizeof(su2double));  

  /* Carry out the matrix vector product. */
  for(int k=0; k<M; ++k) {
    const su2double *AA = A + k*N;
    for(int l=0; l<N; ++l)
      y[k] += AA[l]*x[l];
  }
}
