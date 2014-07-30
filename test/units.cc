#include <iostream>
#include <vector>
#include "fevor_crystal.hh"
#include "test_crystal.hh"

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
    test_grow();
    
    std::cout << "\n"
              << "******************** \n"
              << "BEGIN TEST_DISLOCATE \n"
              << "******************** \n" << std::endl;
    //~ test_dislocate();          
    
    std::cout << "\n"
              << "**************** \n"
              << "BEGIN TEST_MIGRE \n"
              << "**************** \n" << std::endl;
    //~ test_migRe();
    
    std::cout << "\n"
              << "********************* \n"
              << "BEGIN TEST_POLYGONIZE \n"
              << "********************* \n" << std::endl;
    //~ test_polygonize();
    
    std::cout << "\n"
              << "******************* \n"
              << "BEGIN TEST_RESOLVEM \n"
              << "******************* \n" << std::endl;
    //~ test_resolveM();
    
    return 0;
}

