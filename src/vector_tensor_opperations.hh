/*
 */

#ifndef VECTOR_TENSOR_OPPERATIONS
#define VECTOR_TENSOR_OPPERATIONS
 
#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>

// Vector containing the independent elements of the outer product between
// two vectors (upper triangle of resultant tensor) in ROW-MAJOR order.
std::vector<double> vectorOuter(const std::vector<double> &v1, const std::vector<double> &v2);

// Vector containing the elements of the outer product between two tensors 
// in ROW-MAJOR order (M_ijkl = M_1111, M_1112, M_1113,...)
std::vector<double> tensorOuter(std::vector<double> &t1, std::vector<double> &t2);

// Vector containing the elements of the inner product between two tensors 
// in ROW-MAJOR order.
std::vector<double> tensorMixedInner(std::vector<double> &t1, std::vector<double> &t2);

// double containing the magnitude of a tensor in ROW-MAJOR order.
double tensorMagnitude(std::vector<double> &t1);


#endif
