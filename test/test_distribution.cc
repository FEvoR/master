/*
 */

#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include "fevor_distribution.hh"
#include "test_distribution.hh"
#include "vector_tensor_opperations.hh"


void test_stepInTime() {
    fevor_distribution d1({3,3,3}, 1.0); 
    
    //~ d1.saveDistribution();
    
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
}
        
void test_getSoftness() {
    //TODO: test fevor_distribution.getSoftness()
    std::cout << "Test?! Test?! We don't see no stinkin\' Test!" << std::endl;
}
    
void test_setSoftness() {
    //TODO: test fevor_distribution.setSoftnessParam()
    //TODO: rename  setSoftness to setSoftnessParam
    std::cout << "Test?! Test?! We don't see no stinkin\' Test!" << std::endl;
}
    
void test_saveDistribution() {
    //TODO: test fevor_distribution.saveDistribution()
    //TODO: have two overloaded methods: one with no input just prints to std::cout
            // the other takes in an input filename and prints to that file 
            // (options to append, overwrite, etc)
    std::cout << "Test?! Test?! We don't see no stinkin\' Test!" << std::endl;
}
