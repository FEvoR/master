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
    
    double temperature = -10.0;
    std::vector<double> stress = { 10000,     0, 10000,
                                       0,     0,     0,
                                   10000,     0,-10000};
    double modelTime = 0.0;
    double timeStep = 1000.0*365.0*24.0*60.0*60.0;
    double nMigre, nPoly;
    nMigre = nPoly = 0;
    std::vector<double> bulkEdot= { 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0};
    
    std::cout.precision(4);
    std::cout << std::scientific 
              << "Model time: " << modelTime << "\n"
              << "Time Step:  " << timeStep << "\n"
              << "Temp:  " << temperature << "\n"
              << std::endl;
    
    std::cout << "Stress:" << std::endl;
    
    tensorDisplay(stress,3,3);
    
    d1.stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
    
    std::cout.precision(4);
    std::cout << std::scientific 
              << "Model time: " << modelTime << "\n"
              << "Number of MigRes: " << nMigre << "\n"
              << "Number of Polys:  " << nPoly << "\n"
              << std::endl;
    
    std::cout << "edot:" << std::endl;
        
    tensorDisplay(bulkEdot,3,3);
    
    return 0;
}

