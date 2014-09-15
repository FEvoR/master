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
#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>
#include "fevor_distribution.hh"
#include "vector_tensor_operations.hh"

// Perform a compairison time step between FEvoR and Thor
int main()
{
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
    
    std::vector<unsigned int> packingDimensions = {20,20,20};
    //~ double watsonK = 0.0;      // isotropic 
    //~ double watsonK = 10000.0;  // perfect girdle
    //~ double watsonK = -10000.0; // perfect bipolar
    //~ double watsonK = 2.0;      // girdle
    //~ double watsonK = -2.0;     // bipolar
    //~ "./util/compairisonDist.csv"
    fevor_distribution d1(packingDimensions, "./util/compairisonDist.csv");
    
    std::cout << "\n" << "************************************ \n"
                    << "Compairison timeStep: FEvoR vs Thor. \n"
                    << "************************************ \n"
                    << "\n"
                    << "Model Setup is: \n"
                    << "    " << "Softness Parameters: 1, 0 \n"
                    << "    " << "Temperature: " << temperature << " degrees C\n"
                    << "    " << "Stress (Pa): " << std::endl;
                    tensorDisplay(stress, 3, 3);
    std::cout << "\n" << "Initial time: " << modelTime << "s \n"
                    << "Time Step: " << timeStep << "s \n"
                    << std::endl;
    
    std::cout << "\n" << "Got initial distribution! \n" << std::endl;
    
    start = clock();
    
    d1.stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
    
    end = clock();
    msecs = ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;
    
    std::cout << "\n" << "Stepped Distribution! \n" << std::endl;
    d1.saveDistribution("./util/compairisonDist_stepped_cpp.csv");
    
    std::cout << "Bulk Edot was:" << std::endl;
    tensorDisplay(bulkEdot,3,3);
    std::cout << "Number of Polygonizations:" << nPoly << std::endl;
    std::cout << "Number of Migration Recrystallizations:" << nMigre << std::endl;
    
    
    
        
    std::cout << "\n" << "Time to preform a timestep: " << msecs << " ms \n" << std::endl;
    
    
    return 0;
}
