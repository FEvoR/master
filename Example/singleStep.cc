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

/* This is an example of using FEvoR to preform a single time step.
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <string>
#include "fevor_distribution.hh"
#include "vector_tensor_operations.hh"

int main()
{
    // to find the execution time of a single time step
    clock_t start, end;
    double msecs;

    
    double temperature = -10.0;
    std::vector<double> stress = { 10000,     0, 10000,
                                       0,     0,     0,
                                   10000,     0,-10000};
    double modelTime = 0.0;
    double timeStep = 10000.0*365.0*24.0*60.0*60.0;
    
    double nMigre, nPoly;
    nMigre = nPoly = 0;
    std::vector<double> bulkEdot(9, 0.0);
    
    std::vector<unsigned int> packingDimensions = {3,3,3};
    
    // Create a random distribution using the Watson ODF
        double watsonK = 0.0;      // isotropic 
        //~ double watsonK = 10000.0;  // perfect girdle
        //~ double watsonK = -10000.0; // perfect bipolar
        //~ double watsonK = 2.0;      // girdle
        //~ double watsonK = -2.0;     // bipolar
    fevor_distribution d1(packingDimensions, watsonK);
    
    // load a distribution from a CSV file named FILENAME.csv
    //~ fevor_distribution d1(packingDimensions, "FILENAME.csv");
    
    std::cout << "\n" << "************************************ \n"
                      << "Example time step. \n"
                      << "************************************ \n"
                      << "\n"
                      << "Model Setup is: \n"
                      << "    " << "Softness Parameters: 1, 0 \n"
                      << "    " << "Temperature: " << temperature << " degrees C\n"
                      << "    " << "Stress (Pa): " << std::endl;
                    tensorDisplay(stress, 3, 3);
    std::cout << "\n" << "    " << "Initial time: " << modelTime << "s \n"
                      << "    " "Time Step: " << timeStep << "s \n"
                      << std::endl;
    
    std::string saveInitial = "./Example/distribution_initial.csv";
    d1.saveDistribution(saveInitial);
    
    std::cout << "\n" << "Initial distribution saved in: " 
                      << saveInitial << "\n"
                      << std::endl;
    
    start = clock();
    
    d1.stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
    
    end = clock();
    msecs = ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;
    
    std::string saveFinal = "./Example/distribution_final.csv";
    d1.saveDistribution(saveFinal);
    
    std::cout << "\n" << "Stepped Distribution! Saved in: " 
                      << saveFinal << "\n" 
                      << std::endl;
    
    std::cout << "Bulk Edot was:" << std::endl;
    tensorDisplay(bulkEdot,3,3);
    std::cout << "Number of Polygonizations:" << nPoly << std::endl;
    std::cout << "Number of Migration Recrystallizations:" << nMigre << std::endl;
    
    
    
        
    std::cout << "\n" << "Time to preform a timestep: " << msecs << " ms \n" << std::endl;
    
    
    return 0;
}
