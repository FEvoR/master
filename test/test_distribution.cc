/**************************************************
 * FEvoR: Fabric Evolution with Recrystallization *
 **************************************************
 * Copyright (C) 2009-2014  Joseph H Kennedy
 *
 * This file is part of FEvoR.
 *
 * FEvoR is free software: you can redistribute it and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * FEvoR is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with FEvoR.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Additional permission under GNU GPL version 3 section 7
 *
 * If you modify FEvoR, or any covered work, to interface with
 * other modules (such as MATLAB code and MEX-files) available in a
 * MATLAB(R) or comparable environment containing parts covered
 * under other licensing terms, the licensors of FEvoR grant
 * you additional permission to convey the resulting work.
 */
/*
 */

#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include "fevor_distribution.hh"
#include "test_distribution.hh"
#include "vector_tensor_operations.hh"

namespace FEvoR {

void test_stepInTime() {
    FEvoR::Distribution d1({3,3,3}, 1.0); 
    
    //~ d1.saveDistribution();
    
    double temperature = 273.15-10.0;
    std::vector<double> stress = { 10000,     0, 10000,
                                       0,     0,     0,
                                   10000,     0,-10000};
    double modelTime = 0.0;
    double timeStep = 1000.0*365.0*24.0*60.0*60.0;
    unsigned int nMigre, nPoly;
    nMigre = nPoly = 0;
    std::vector<double> bulkEdot= { 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0};
    
    std::cout.precision(4);
    std::cout << std::scientific 
              << "Model time: " << modelTime << "s \n"
              << "Time Step:  " << timeStep << "s \n"
              << "Temp:  " << temperature << "K \n"
              << std::endl;
    
    std::cout << "Stress (Pa):" << std::endl;
    
    tensorDisplay(stress,3,3);
    
    d1.stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
    modelTime+=timeStep;
    
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
    //TODO: test FEvoR::distribution.getSoftness()
    std::cout << "Test?! Test?! We don't see no stinkin\' Test!" << std::endl;
}
    
void test_setSoftnessRatio() {
    //TODO: test FEvoR::distribution.setSoftnessParam()
    //TODO: rename  setSoftness to setSoftnessParam
    std::cout << "Test?! Test?! We don't see no stinkin\' Test!" << std::endl;
}
    
void test_saveDistribution() {
    //TODO: test FEvoR::distribution.saveDistribution()
    //TODO: have two overloaded methods: one with no input just prints to std::cout
            // the other takes in an input filename and prints to that file 
            // (options to append, overwrite, etc)
    std::cout << "Test?! Test?! We don't see no stinkin\' Test!" << std::endl;
}

} // end of namespace FEvoR
