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

/* This is an example of using FEvoR to calculate an enhancement factor.
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <ctime>
#include <string>
#include "fevor_distribution.hh"
#include "vector_tensor_operations.hh"

int main(int argc, char *argv[])
{
    if (argc != 10) {
        std::cout << "usage: " << argv[0] << " w T Sxx(Pa) Sxy Sxz Syy Syz Szz dt(years)" << std::endl;
        std::cout << "example: " << argv[0] << "-3.0 -30.0 0.0 0.0 500 0.0 0.0 0.0 dt(years)" << std::endl;
        return 0;
    }
    
    
    double watsonK, temperature, s11, s12, s13, s22, s23, s33, dt;
    std::stringstream ss;
    ss << argv[1]; ss >> watsonK;
    ss.clear();
    ss << argv[2]; ss >> temperature;
    ss.clear();
    ss << argv[3]; ss >> s11;
    ss.clear();
    ss << argv[4]; ss >> s12;
    ss.clear();
    ss << argv[5]; ss >> s13;
    ss.clear();
    ss << argv[6]; ss >> s22;
    ss.clear();
    ss << argv[7]; ss >> s23;
    ss.clear();
    ss << argv[8]; ss >> s33;
    ss.clear();
    ss << argv[9]; ss >> dt;
    
    
    // to find the execution time of a single time step
    clock_t start, end;
    double msecs;

    std::vector<double> stress = { s11, s12, s13,
                                   s12, s22, s23,
                                   s13, s23, s33};
    
    double modelTime = 0.0;
    double modelTime_iso = 0.0;
    double timeStep = dt*365.0*24.0*60.0*60.0;
    
    unsigned int nMigre, nPoly;
    nMigre = nPoly = 0;
    std::vector<double> bulkEdot(9, 0.0);
    std::vector<double> bulkEdot_iso(9, 0.0);
    
    std::vector<unsigned int> packingDimensions = {10,10,10};
    
    // Create a random distribution using the Watson ODF
        //~ double watsonK = 0.0;      // isotropic 
        //~ double watsonK = 10000.0;  // perfect girdle
        //~ double watsonK = -10000.0; // perfect bipolar
        //~ double watsonK = 5.0;      // girdle
        //~ double watsonK = -2.0;     // bipolar
    FEvoR::Distribution d1(packingDimensions, watsonK);
    FEvoR::Distribution d_iso(packingDimensions,0.0);
    
    
    std::cout << "\n" << "************************************ \n"
                      << "Example time step. \n"
                      << "************************************ \n"
                      << "\n"
                      << "Model Setup is: \n"
                      << "    " << "Softness Parameters: 1, 0 \n"
                      << "    " << "Temperature: " << temperature << " degrees C\n"
                      << "    " << "Stress (Pa): " << std::endl;
                    FEvoR::tensorDisplay(stress, 3, 3);
    std::cout << "\n" << "    " << "Initial time: " << modelTime << "s \n"
                      << "    " "Time Step: " << timeStep << "s \n"
                      << std::endl;
    
    std::string saveInitial = "distribution_initial.csv";
    d1.saveDistribution(saveInitial);
    std::string saveInitial_iso = "distribution_initial_iso.csv";
    d1.saveDistribution(saveInitial_iso);
    
    std::cout << "\n" << "Initial distribution saved in: \n" 
                      << saveInitial << "\n"
                      << saveInitial_iso << "\n"
                      << std::endl;
    
    start = clock();
    
    d1.stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
    d_iso.stepInTime(temperature, stress, modelTime_iso, timeStep, bulkEdot_iso);
    
    
    double E = 0.0,
           M = 0.0,
       M_iso = 0.0;
    
    M = FEvoR::tensorMagnitude(bulkEdot);
    M_iso = FEvoR::tensorMagnitude(bulkEdot_iso);
    
    if (M_iso != 0.0) {
        E = M / M_iso;
      } else {
        E = 1.0;
      }
    
    end = clock();
    msecs = ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;
    
    std::string saveFinal = "distribution_final.csv";
    d1.saveDistribution(saveFinal);
    std::string saveFinal_iso = "distribution_final_iso.csv";
    d_iso.saveDistribution(saveFinal_iso);
    
    std::cout << "\n" << "Stepped Distribution! Saved in: \n" 
                      << saveFinal << "\n" 
                      << saveFinal_iso << "\n" 
                      << std::endl;
    
    std::cout << "Bulk Edot was:" << std::endl;
    FEvoR::tensorDisplay(bulkEdot,3,3);
    std::cout << "Bulk Edot_iso was:" << std::endl;
    FEvoR::tensorDisplay(bulkEdot_iso,3,3);

    std::cout << "Number of Polygonizations:" << nPoly << std::endl;
    std::cout << "Number of Migration Recrystallizations:" << nMigre << std::endl;
    
    std::cout << "M was:" << M << std::endl;
    std::cout << "M_iso was:" << M_iso << std::endl;
    std::cout << "Enhancement factor was:" << E << std::endl;
        
    std::cout << "\n" << "Time to preform a timestep: " << msecs << " ms \n" << std::endl;
    
    
    return 0;
}
