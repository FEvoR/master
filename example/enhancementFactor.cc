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
#include <cmath>
#include <numeric>
#include "fevor_distribution.hh"
#include "vector_tensor_operations.hh"

int main(int argc, char *argv[])
{
    if (argc != 12) {
        std::cout << "usage: " << argv[0] << " w T(K) Sxx(Pa) Sxy Sxz Syy Syz Szz dt(years) cc cn" << std::endl;
        std::cout << "example: " << argv[0] << " -3.0 263.15 0.0 0.0 500 0.0 0.0 0.0 1000.0 1.0 0.0" << std::endl;
        return 0;
    }
    
    
    double watsonK, temperature, sxx, sxy, sxz, syy, syz, szz, dt, cc, cn;
    std::stringstream ss;
    ss << argv[1]; ss >> watsonK;
    ss.clear();
    ss << argv[2]; ss >> temperature;
    ss.clear();
    ss << argv[3]; ss >> sxx;
    ss.clear();
    ss << argv[4]; ss >> sxy;
    ss.clear();
    ss << argv[5]; ss >> sxz;
    ss.clear();
    ss << argv[6]; ss >> syy;
    ss.clear();
    ss << argv[7]; ss >> syz;
    ss.clear();
    ss << argv[8]; ss >> szz;
    ss.clear();
    ss << argv[9]; ss >> dt;
    ss.clear();
    ss << argv[10]; ss >> cc;
    ss.clear();
    ss << argv[11]; ss >> cn;
    
    
    // to find the execution time of a single time step
    clock_t start, end;
    double msecs;

    std::vector<double> stress = { sxx, sxy, sxz,
                                   sxy, syy, syz,
                                   sxz, syz, szz};
    
    double timeStep = dt*365.0*24.0*60.0*60.0;
    
    unsigned int nMigre = 0, 
                 nPoly  = 0;
    
    std::vector<double> bulkEdot(9, 0.0);
    std::vector<double> bulkEdot_iso(9, 0.0);
    
    std::vector<unsigned int> packingDimensions = {5,5,5};
    

    FEvoR::Distribution d1(packingDimensions, watsonK);
    FEvoR::Distribution di(packingDimensions,0.0);
    d1.setSoftnessRatio(cc, cn);
    di.setSoftnessRatio(cc, cn);
    
    d1.saveDistribution("eDist_wk_initial.csv");
    di.saveDistribution("eDist_iso_initial.csv");
    
    std::cout << "\n" << "************************************ \n"
                      << "Example time step. \n"
                      << "************************************ \n"
                      << "\n"
                      << "Model Setup is: \n"
                      << "    " << "Watson K: " << watsonK << " \n"
                      << "    " << "Softness Parameters: " << cc << " , " << cn << " \n"
                      << "    " << "Temperature: " << temperature << " Kelvin\n"
                      << "    " << "Stress (Pa): " << std::endl;
    FEvoR::tensorDisplay(stress, 3, 3);
    std::cout         << "    " << "Time Step: " << timeStep << "s \n"
                      << std::endl;
    
    start = clock();
    
    double modelTime = 0.0;
    d1.stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
    di.stepInTime(temperature, stress, modelTime, timeStep, bulkEdot_iso);
    
    
    d1.saveDistribution("eDist_wk_final.csv");
    di.saveDistribution("eDist_iso_final.csv");
    
    double M = 0.0,
       M_iso = 0.0;
    M     = FEvoR::tensorMagnitude(bulkEdot);
    M_iso = FEvoR::tensorMagnitude(bulkEdot_iso);
    
    end = clock();
    msecs = ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;
    
    std::cout << "\n" << "Time to preform a timestep: " << msecs << " ms \n" << std::endl;
    
    
    std::cout << "Bulk Edot was:" << std::endl;
    FEvoR::tensorDisplay(bulkEdot,3,3);
    std::cout << "Bulk Edot_iso was:" << std::endl;
    FEvoR::tensorDisplay(bulkEdot_iso,3,3);
    std::cout << "M            was:" << M     << std::endl;
    std::cout << "M_iso        was:" << M_iso << "\n" << std::endl;
    
    double E_mag = M/M_iso;
    std::cout << "E_mag        was:" << std::fixed << E_mag << "\n" << std::endl;
    
    /* Indexing: {0, 1, 2,
     *            3, 4, 5,
     *            6, 7, 8}
     */
    double  E_33 = 0.0,
            E_13 = 0.0,
        E_weights= 0.0;
    
    if (bulkEdot_iso[8] != 0.0)
        E_33 = std::abs(bulkEdot[8]/bulkEdot_iso[8]);
    std::cout << "E_33         was:" << std::fixed << E_33 << "\n" << std::endl;
    
    if (bulkEdot_iso[2] != 0.0)
        E_13 = std::abs(bulkEdot[2]/bulkEdot_iso[2]);        
    std::cout << "E_13         was:" << std::fixed << E_13 << "\n" << std::endl;
    
    if (std::abs(stress[2])+std::abs(stress[8]) != 0.0 )
        E_weights = (E_13*std::abs(stress[2])+E_33*std::abs(stress[8]))/(std::abs(stress[2])+std::abs(stress[8]));
    std::cout << "E_weights     was:" << std::fixed << E_weights << "\n" << std::endl;
    
    return 0;
}
