/**************************************************
 * FEvoR: Fabric Evolution with Recrystallization *
 **************************************************
 * Copyright (C) 2009-2014  Joseph H Kennedy
 *
 * This file is part of FEvoR.
 *
 * FEvoR is free software: you can redistribute it and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * FEvoR is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with FEvoR.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Additional permission under GNU GPL version 3 section 7
 *
 * If you modify FEvoR, or any covered work, to interface with
 * other modules (such as MATLAB code and MEX-files) available in a
 * MATLAB(R) or comparable environment containing parts covered
 * under other licensing terms, the licensors of FEvoR grant
 * you additional permission to convey the resulting work.
 */
/*
 */

#ifndef VECTOR_TENSOR_OPPERATIONS
#define VECTOR_TENSOR_OPPERATIONS
 
#include <vector>

namespace FEvoR {
    
// Vector containing the independent elements of the outer product between
// two vectors (upper triangle of resultant tensor) in ROW-MAJOR order.
std::vector<double> vectorOuter(const std::vector<double> &v1, const std::vector<double> &v2);

// Vector containing the elements of the outer product between two tensors 
// in ROW-MAJOR order (M_ijkl = M_1111, M_1112, M_1113,...)
std::vector<double> tensorOuter(const std::vector<double> &t1, const std::vector<double> &t2);

// Vector containing the elements of the inner product between two tensors 
// in ROW-MAJOR order.
std::vector<double> tensorMixedInner(const std::vector<double> &t1, const std::vector<double> &t2);

// double containing the magnitude of a tensor in ROW-MAJOR order.
double tensorMagnitude(const std::vector<double> &t1);

// transpose a matrix
std::vector<double> matrixTranspose(const std::vector<double> &matrix, const int rows, const int columns);

// get the taylor series of expm()
std::vector<double> vectorTimesExpm(const std::vector<double> &vec, const std::vector<double> &m, int terms);

std::vector<double> matrixPowTimesVector(std::vector<double> m, int power, std::vector<double> vec);

int factorial(int n);

// display a tensor
void tensorDisplay(const std::vector<double> &tensor, const int rows, const int columns);

} // end of namespace FEvoR

#endif
