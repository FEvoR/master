/*
 */

#ifndef VECTOR_TENSOR_OPPERATIONS
#define VECTOR_TENSOR_OPPERATIONS
 
#include <vector>

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



#endif
