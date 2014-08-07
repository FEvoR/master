#include <iostream>
#include <vector>
#include "fevor_crystal.hh"
#include "test_crystal.hh"
#include "fevor_distribution.hh"
#include "test_distribution.hh"
#include "vector_tensor_opperations.hh"

// Test the interface to fevor_crystal class
int main()
{
    std::cout << "\n"
              << "****************** \n"
              << "BEGIN TEST_CRYSTAL \n"
              << "****************** \n" << std::endl;
              
        std::cout << "\n"
                  << "    ***************** \n"
                  << "    BEGIN TEST_ANGLES \n"
                  << "    ***************** \n" << std::endl;
        test_angles();
        
        std::cout << "\n"
                  << "    *************** \n"
                  << "    BEGIN TEST_GROW \n"
                  << "    *************** \n" << std::endl;
        double K;
        K = test_grow();
        
        std::cout << "    Growth Constant: " << K
                  << " \n" << std::endl;
        
        std::cout << "\n"
                  << "    ******************** \n"
                  << "    BEGIN TEST_DISLOCATE \n"
                  << "    ******************** \n" << std::endl;
        test_dislocate( K );          
        
        std::cout << "\n"
                  << "    **************** \n"
                  << "    BEGIN TEST_MIGRE \n"
                  << "    **************** \n" << std::endl;
        test_migRe();
        
        std::cout << "\n"
                  << "    ********************* \n"
                  << "    BEGIN TEST_POLYGONIZE \n"
                  << "    ********************* \n" << std::endl;
        test_polygonize();
        
        std::cout << "\n"
                  << "    ******************* \n"
                  << "    BEGIN TEST_RESOLVEM \n"
                  << "    ******************* \n" << std::endl;
        test_resolveM();
        
        std::cout << "\n"
                  << "    ***************** \n"
                  << "    BEGIN TEST_ROTATE \n"
                  << "    ***************** \n" << std::endl;
        test_rotate();
        
    std::cout << "\n"
              << "*********************** \n"
              << "BEGIN TEST_DISTRIBUTION \n"
              << "*********************** \n" << std::endl;
        
        std::cout << "\n"
                  << "    ********************* \n"
                  << "    BEGIN TEST_STEPINTIME \n"
                  << "    ********************* \n" << std::endl;
        test_stepInTime();
        
        std::cout << "\n"
                  << "    ********************** \n"
                  << "    BEGIN TEST_GETSOFTNESS \n"
                  << "    ********************** \n" << std::endl;
        test_getSoftness();
        
        std::cout << "\n"
                  << "    ********************** \n"
                  << "    BEGIN TEST_SETSOFTNESS \n"
                  << "    ********************** \n" << std::endl;    
        test_setSoftness();
            
        std::cout << "\n"
                  << "    *************************** \n"
                  << "    BEGIN TEST_SAVEDISTRIBUTION \n"
                  << "    *************************** \n" << std::endl;
        test_saveDistribution();
    
    std::cout << "\n"
              << "********* \n"
              << "END TESTS \n"
              << "********* \n" << std::endl;
    
    std::cout << "Your tests suck! Units out." << "\n" << std::endl;
    
    return 0;
}

