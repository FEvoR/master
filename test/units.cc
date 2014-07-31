#include <iostream>
#include <vector>
#include "fevor_crystal.hh"
#include "test_crystal.hh"
#include "fevor_distribution.hh"

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
        
        fevor_distribution d1(10, 1.0); 
        
    
    return 0;
}

