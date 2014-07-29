/*
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>

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
std::vector<double> tensorOuter(std::vector<double> &t1, std::vector<double> &t2) {
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
std::vector<double> tensorMixedInner(std::vector<double> &t1, std::vector<double> &t2) {
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
double tensorMagnitude(std::vector<double> &t1) {
    double mag = 0;
    for (auto &ii : t1)
        mag += ii*ii;
        
    return sqrt(mag);
}
