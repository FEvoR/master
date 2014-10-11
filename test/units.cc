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
#include <vector>
#include "fevor_crystal.hh"
#include "test_crystal.hh"
#include "fevor_distribution.hh"
#include "test_distribution.hh"
#include "vector_tensor_operations.hh"

// Test the interface to FEvoR::crystal class
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
        FEvoR::test_angles();
        
        std::cout << "\n"
                  << "    *************** \n"
                  << "    BEGIN TEST_GROW \n"
                  << "    *************** \n" << std::endl;
        double K;
        K = FEvoR::test_grow();
        
        std::cout << "    Growth Constant: " << K
                  << " \n" << std::endl;
        
        std::cout << "\n"
                  << "    ******************** \n"
                  << "    BEGIN TEST_DISLOCATE \n"
                  << "    ******************** \n" << std::endl;
        FEvoR::test_dislocate( K );          
        
        std::cout << "\n"
                  << "    **************** \n"
                  << "    BEGIN TEST_MIGRE \n"
                  << "    **************** \n" << std::endl;
        FEvoR::test_migRe();
        
        std::cout << "\n"
                  << "    ********************* \n"
                  << "    BEGIN TEST_POLYGONIZE \n"
                  << "    ********************* \n" << std::endl;
        FEvoR::test_polygonize();
        
        std::cout << "\n"
                  << "    ******************* \n"
                  << "    BEGIN TEST_RESOLVEM \n"
                  << "    ******************* \n" << std::endl;
        FEvoR::test_resolveM();
        
        std::cout << "\n"
                  << "    ***************** \n"
                  << "    BEGIN TEST_ROTATE \n"
                  << "    ***************** \n" << std::endl;
        FEvoR::test_rotate();
        
    std::cout << "\n"
              << "*********************** \n"
              << "BEGIN TEST_DISTRIBUTION \n"
              << "*********************** \n" << std::endl;
        
        std::cout << "\n"
                  << "    ********************* \n"
                  << "    BEGIN TEST_STEPINTIME \n"
                  << "    ********************* \n" << std::endl;
        FEvoR::test_stepInTime();
        
        std::cout << "\n"
                  << "    ********************** \n"
                  << "    BEGIN TEST_GETSOFTNESS \n"
                  << "    ********************** \n" << std::endl;
        FEvoR::test_getSoftness();
        
        std::cout << "\n"
                  << "    *************************** \n"
                  << "    BEGIN TEST_SETSOFTNESSRATIO \n"
                  << "    *************************** \n" << std::endl;    
        FEvoR::test_setSoftnessRatio();
            
        std::cout << "\n"
                  << "    *************************** \n"
                  << "    BEGIN TEST_SAVEDISTRIBUTION \n"
                  << "    *************************** \n" << std::endl;
        FEvoR::test_saveDistribution();
    
    std::cout << "\n"
              << "********* \n"
              << "END TESTS \n"
              << "********* \n" << std::endl;
    
    std::cout << "Your tests suck! Units out." << "\n" << std::endl;
    
    return 0;
}

