/**************************************************
 * FEvoR: Fabric Evolution with Recrystallization *
 **************************************************
 * Copyright (C) 2009-2017  Joseph H. Kennedy
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

#include "catch.hpp"
#include "fevor_crystal.hh"
#include "test_crystal.hh"
#include "vector_tensor_operations.hh"

//TODO: develop a robust test suite that includes edge cases

namespace FEvoR {

void test_rotate(){
    
    //TODO: test FEvoR::Crystal.rotate
    std::cout << "Test?! Test?! We don't see no stinkin\' Test!" << std::endl;
    
}

void test_resolveM() {
    FEvoR::Crystal c1(std::vector<double> {0,0,1},0.01,1e11);
    
    double theta, phi;
    theta = M_PI/4.0;
    phi = M_PI/4.0;
    c1.setNewAxis(theta, phi);
    
    std::vector<double> stress = {10000.0,     0.0,     0.0,
                                      0.0,     0.0,     0.0,
                                      0.0,     0.0,-10000.0};
    double temperature = 273.15-10.0;
    double Mrss = 1.0;
    double Medot = 1.0;
    
    std::vector<double> bigM;
    bigM = c1.resolveM(temperature, stress, Mrss, Medot);
    
    std::cout << "bigM size should be: 81" << "\n"
              << "bigM size is: " << bigM.size() << std::endl;
              
    std::vector<double> vel;
    vel = tensorMixedInner(bigM, stress);
    
    std::cout << "\n" << "vel is:" << std::endl;
    tensorDisplay(vel, 9, 1);
    
    std::cout << "\n" << "bigM is:" << std::endl;
    tensorDisplay(bigM, 9, 9);
    
    std::cout << "\n" << "stress is:" << std::endl;
    tensorDisplay(stress, 9, 1);
    
    
}

TEST_CASE("Normal grain growth",  "[crystal]") {
    FEvoR::Crystal c1(std::vector<double> {1.0,0.0,0.0},0.01,1.0e10);
    
    double ca0, ca1, ca2, csz, cdd, ctlr, cslr; 
    
    double temperature, model_time, K;
    temperature = 273.15-11.0; // units: Kelvin
    model_time = 1000.0*365.0*24.0*60.0*60.0;   // units: m
    
    K = c1.grow(temperature, model_time);
    c1.getAll(ca0, ca1, ca2, csz, cdd, ctlr, cslr);

    //onlything changed should be the crystal size
    REQUIRE( ca0 == Approx( 1.0 ) );
    REQUIRE( ca1 == Approx( 0.0 ) );
    REQUIRE( ca2 == Approx( 0.0 ) );
    REQUIRE( cdd == Approx( 1.0e10 ) );
    REQUIRE( ctlr == Approx( 0.0 ) );
    REQUIRE( cslr == Approx( 0.01 ) );
    
    REQUIRE( csz == Approx( 1.0055e-2).epsilon( 1.e-6 ) );
    
    SECTION("Update dislocation density") {
        double Medot = sqrt(1.0/2.0);
        
        c1.dislocate(model_time, Medot, K);
        c1.getAll(ca0, ca1, ca2, csz, cdd, ctlr, cslr);

        //onlything changed should be the dislocation density
        REQUIRE( ca0 == Approx( 1.0 ) );
        REQUIRE( ca1 == Approx( 0.0 ) );
        REQUIRE( ca2 == Approx( 0.0 ) );
        REQUIRE( csz == Approx( 1.0055e-2).epsilon( 1.e-6 ) );
        REQUIRE( ctlr == Approx( 0.0 ) );
        REQUIRE( cslr == Approx( 0.01 ) );
        
        REQUIRE( cdd == Approx( 4.9282e21).epsilon( 1e17 ) );
    }
}


void test_migRe() {
    FEvoR::Crystal c1(std::vector<double> {1.0,0.0,0.0},0.01,1.0e11);
    
    double model_time, time_step;
    model_time = time_step = 1000.0*365.0*24.0*60.0*60.0;
    
    std::vector<double> stress = {10000.0,     0.0,     0.0,
                                      0.0,     0.0,     0.0,
                                      0.0,     0.0,-10000.0};
    
    c1.migRe(stress, model_time, time_step);

    c1.seeCrystal();
    
    std::cout << "New Grain Size should be: 4.6416e-6" << "\n"
              << "New Dislocation Density should be: 1.0000e10" << std::endl;
    
    double theta, phi;
    theta = phi = 0.0;
    c1.getAxisAngles(theta, phi);
    
    std::cout.precision(4);
    
    std::cout << "New Angles are: \n"
              << "    Theta: " << std::fixed << theta << "\n" 
              << "    Phi:   " << phi << "\n "
              << "New Angles should be: \n"
              << "     0.5236 <= theta <= 1.0472 \n"
              << "    -6.2832 <=  phi  <= 6.2832" << std::endl;
}

void test_polygonize() {
    FEvoR::Crystal c1(std::vector<double> {0.0,0.0,1.0},0.01,1.0e11);
    
    double model_time, time_step;
    model_time = time_step = 1000.0*365.0*24.0*60.0*60.0;
    
    std::vector<double> stress = {10000.0,     0.0,     0.0,
                                      0.0,     0.0,     0.0,
                                      0.0,     0.0,-10000.0};
    
    c1.polygonize(stress, 1.0, model_time, time_step);
    c1.seeCrystal();
    
    std::cout << "New Grain Size should be: 5.0000e-3" << "\n"
              << "New Dislocation Density should be: 4.6000e10" << std::endl;
    
    double theta, phi;
    theta = phi = 0.0;
    c1.getAxisAngles(theta, phi);
    
    std::cout.precision(4);
    
    std::cout << "New Angles are: \n"
              << "    Theta: " << std::fixed << theta << "\n" 
              << "    Phi:   " << phi << "\n "
              << "New Angles should be: \n"
              << "     theta == 0.0873 \n"
              << "    -6.2832 <=  phi  <= 6.2832" << std::endl;
    
    
}

void test_angles() {
    
    FEvoR::Crystal c1(std::vector<double> {1.0,0.0,0.0},1.0,10.0);
    c1.seeCrystal();
    
    double theta = 0.0, phi = 0.0;
    c1.getAxisAngles(theta, phi);
    
    std::cout.precision(4);
    std::cout << std::fixed << "Angles are: \n" << "Theta: " << theta << " Phi: " << phi << std::endl;
    
    std::random_device seed;
    std::mt19937 gen(seed());
    
    std::uniform_real_distribution<double> dPhi(0.0,2.0*M_PI);
    std::uniform_real_distribution<double> dTheta(0.0,M_PI/2.0);
    
    theta = dTheta(gen);
    phi = dPhi(gen);
    
    std::cout.precision(4);
    std::cout << std::fixed << "New angles are: \n" << "Theta: " << theta << " Phi: " << phi << std::endl;
    
    c1.setNewAxis(theta, phi);
    
    c1.seeCrystal();
    
    
}

} // end of namspace FEvoR
