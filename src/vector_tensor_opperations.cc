/*
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "vector_tensor_opperations.hh"

// Vector containing the independent elements of the outer product between
// two vectors (upper triangle of resultant tensor) in ROW-MAJOR order.
std::vector<double> vectorOuter(const std::vector<double> &v1, const std::vector<double> &v2) {
    std::vector<double> tensor;



    for (decltype(v1.size()) ii = 0; ii != v1.size(); ++ii) {
        for (decltype(v1.size()) jj = ii; jj != v2.size(); ++jj) {
            tensor.push_back(v1[ii]*v2[jj]);
        }
    }
    
    return tensor;
}

// Vector containing the elements of the outer product between two tensors 
// in ROW-MAJOR order (M_ijkl = M_1111, M_1112, M_1113,...)
std::vector<double> tensorOuter(const std::vector<double> &t1, const std::vector<double> &t2) {
    std::vector<double> tensor;

    for (auto &ii : t1) {
        for (auto &jj : t2) {
            tensor.push_back(ii*jj);
        }
    }
    
    return tensor;
}


// Vector containing the elements of the inner product between two tensors  
// in ROW-MAJOR order of mixed rank.
std::vector<double> tensorMixedInner(const std::vector<double> &t1, const std::vector<double> &t2) {
    std::vector<double> tensor;
    decltype(t2.size()) sz2 = t2.size();
    
    for (decltype(t1.size()) ii = 0; ii != t1.size(); ii+=sz2) {
        
        tensor.push_back(0);
        
        for (decltype(t2.size()) jj = 0; jj != t2.size(); ++jj) {
        
            tensor.back() += t1[ii+jj]*t2[jj];
            
        }
    }
    
    return tensor;
}

// magnitude of a full tensor. 
double tensorMagnitude(const std::vector<double> &t1) {
    double mag = 0;
    for (auto &ii : t1)
        mag += ii*ii;
        
    return sqrt(mag);
}

std::vector<double> matrixTranspose(const std::vector<double> &m1, const int rows, const int columns) {
    std::vector<double> matrix;
    for (int ii = 0; ii != columns; ++ii) {
        for (int jj = 0; jj != rows; ++jj) {
            matrix.push_back(m1[ii+jj*3]);
        }
    }
    
    return matrix;
    
}

std::vector<double> vectorTimesExpm(const std::vector<double> &vec, const std::vector<double> &m, int terms){
    
    std::vector<double> sum = vec;
    std::vector<double> tempMat;
    
    for (int ii = 1; ii != terms; ++ii ) {
        tempMat = matrixPowTimesVector(m, ii, vec);
        std::transform(tempMat.begin(), tempMat.end(),tempMat.begin(), 
                       [&](double x){return x/factorial(ii);} );
        std::transform(sum.begin(), sum.end(),tempMat.begin(),sum.begin(), 
                       std::plus<double>() );
    }
    
    return sum;
}

std::vector<double> matrixPowTimesVector(std::vector<double> m, int power, std::vector<double> vec) {
    // tried recursion first -- crashes on power > 3
    std::vector<double> temp = vec;
    for (int ii = 0; ii != power; ++ii) {
        temp = tensorMixedInner(m, temp);
    }    
    return temp;
}

int factorial(int n) {
    if (n > 1)
        return factorial(n-1) * n;
    
    return 1;
}

void tensorDisplay(const std::vector<double> &tensor, const int rows, const int columns) {
    
    std::cout.precision(4);
    
    for (int ii = 0; ii != rows; ++ii) {
        for (int jj = 0; jj != columns; ++jj) {
            
            if (ii == 0 && jj == 0){
                std::cout << std::scientific
                          << tensor[ii*columns+jj] << " ";
            } else {
            
            std::cout << tensor[ii*columns+jj] << " ";
            
            }
            
            if (jj == columns-1)
                std::cout << "\n";
            
        }
    }
    
    
    std::cout << std:: endl;
    
}
