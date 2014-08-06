#include <iostream>
#include <vector>
#include "fevor_crystal.hh"
#include "test_crystal.hh"
#include "fevor_distribution.hh"
#include "vector_tensor_opperations.hh"

// Test the interface to fevor_crystal class
int main()
{
    std::cout << "\n"
              << "***************** \n"
              << "BEGIN TEST_ANGLES \n"
              << "***************** \n" << std::endl;
    test_angles();
    
    std::cout << "\n"
              << "*************** \n"
              << "BEGIN TEST_GROW \n"
              << "*************** \n" << std::endl;
    double K;
    K = test_grow();
    
    std::cout << "Growth Constant: " << K
              << " \n" << std::endl;
    
    std::cout << "\n"
              << "******************** \n"
              << "BEGIN TEST_DISLOCATE \n"
              << "******************** \n" << std::endl;
    test_dislocate( K );          
    
    std::cout << "\n"
              << "**************** \n"
              << "BEGIN TEST_MIGRE \n"
              << "**************** \n" << std::endl;
    test_migRe();
    
    std::cout << "\n"
              << "********************* \n"
              << "BEGIN TEST_POLYGONIZE \n"
              << "********************* \n" << std::endl;
    test_polygonize();
    
    std::cout << "\n"
              << "******************* \n"
              << "BEGIN TEST_RESOLVEM \n"
              << "******************* \n" << std::endl;
    test_resolveM();
    
    std::cout << "\n"
          << "*********************** \n"
          << "BEGIN TEST_DISTRIBUTION \n"
          << "*********************** \n" << std::endl;
    
    fevor_distribution d1({3,3,3}, 1.0); 
    
    d1.saveDistribution();
        
    
    std::cout << "\n"
          << "*************************** \n"
          << "BEGIN TEST_MATRIX_TRANSPOSE \n"
          << "*************************** \n" << std::endl;
        
        std::vector<double> matrix = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0};
        tensorDisplay(matrix, 3, 3);
        
        matrix = matrixTranspose(matrix, 3, 3);
        tensorDisplay(matrix, 3, 3);
        
    
    return 0;
}

