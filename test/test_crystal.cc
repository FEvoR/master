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
#include "vector_tensor_operations.hh"

//TODO: develop a robust test suite that includes edge cases

namespace FEvoR {

TEST_CASE("Rotate", "[crystal]") {
    //TODO: test FEvoR::Crystal.rotate
}


TEST_CASE("Resolve M",  "[crystal]") {
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
   
    REQUIRE( bigM.size() == 81 );
    //TODO: actual tests for bigM
}


TEST_CASE("Normal grain growth",  "[crystal]") {
    FEvoR::Crystal c1(std::vector<double> {1.0,0.0,0.0},0.01,1.0e10);
    
    double ca0, ca1, ca2, csz, cdd, ctlr, cslr; 
    
    double temperature, model_time, K;
    temperature = 273.15-11.0; // units: Kelvin
    model_time = 1000.0*365.0*24.0*60.0*60.0;   // units: m
    
    K = c1.grow(temperature, model_time);
    c1.getAll(ca0, ca1, ca2, csz, cdd, ctlr, cslr);

    //only thing changed should be the crystal size
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

        //only thing changed should be the dislocation density
        REQUIRE( ca0 == Approx( 1.0 ) );
        REQUIRE( ca1 == Approx( 0.0 ) );
        REQUIRE( ca2 == Approx( 0.0 ) );
        REQUIRE( csz == Approx( 1.0055e-2).epsilon( 1.e-6 ) );
        REQUIRE( ctlr == Approx( 0.0 ) );
        REQUIRE( cslr == Approx( 0.01 ) );
        
        REQUIRE( cdd == Approx( 4.9282e21).epsilon( 1e17 ) );
    }
}


TEST_CASE("Migration recrystallization",  "[crystal]") {
    FEvoR::Crystal c1(std::vector<double> {1.0,0.0,0.0},0.01,1.0e11);
    
    double ca0, ca1, ca2, csz, cdd, ctlr, cslr; 
    
    double model_time, time_step;
    model_time = time_step = 1000.0*365.0*24.0*60.0*60.0;
    
    std::vector<double> stress = {10000.0,     0.0,     0.0,
                                      0.0,     0.0,     0.0,
                                      0.0,     0.0,-10000.0};
    
    c1.migRe(stress, model_time, time_step);
    c1.getAll(ca0, ca1, ca2, csz, cdd, ctlr, cslr);

    REQUIRE( csz == Approx(4.6416e-6).epsilon(1.0e-10) );
    REQUIRE( cdd == Approx(1.0e10) );
    
    double theta, phi;
    c1.getAxisAngles(theta, phi);
    
    REQUIRE( theta == Approx(M_PI/4.0).epsilon(M_PI/12.0));
    REQUIRE( phi == Approx(0.0).epsilon(M_PI) );
}


TEST_CASE("Polygonize",  "[crystal]") {
    FEvoR::Crystal c1(std::vector<double> {0.0,0.0,1.0},0.01,1.0e11);
    
    double ca0, ca1, ca2, csz, cdd, ctlr, cslr; 
    
    double model_time, time_step;
    model_time = time_step = 1000.0*365.0*24.0*60.0*60.0;
    
    std::vector<double> stress = {10000.0,     0.0,     0.0,
                                      0.0,     0.0,     0.0,
                                      0.0,     0.0,-10000.0};
    
    c1.polygonize(stress, 1.0, model_time, time_step);
    c1.getAll(ca0, ca1, ca2, csz, cdd, ctlr, cslr);
   
    REQUIRE( csz == Approx(5.0e-3).epsilon(1.0e-7) );
    REQUIRE( cdd == Approx(4.6e10).epsilon(1.0e6) );

    double theta, phi;
    c1.getAxisAngles(theta, phi);
    
    REQUIRE( theta == Approx(M_PI/36.0) ); // 5 degrees
    REQUIRE( phi == Approx(0.0).epsilon(M_PI) );
}


TEST_CASE("Crystal angles",  "[crystal]") {
    FEvoR::Crystal c1(std::vector<double> {1.0,0.0,0.0},1.0,10.0);
   
    double theta = 0.0, phi = 0.0;
    c1.getAxisAngles(theta, phi);

    REQUIRE( theta == Approx( M_PI/2.0 ) );
    REQUIRE( phi == Approx( 0.0 ) );

    c1.setNewAxis(M_PI/2.0, M_PI/2.0);
    
    double ca0, ca1, ca2; 
    c1.getAxis(ca0, ca1, ca2);

    REQUIRE( ca0 == Approx(0.0) );
    REQUIRE( ca1 == Approx(1.0) );
    REQUIRE( ca2 == Approx(0.0) );

    std::random_device seed;
    std::mt19937 gen(seed());
    
    std::uniform_real_distribution<double> dPhi(-M_PI,M_PI);
    std::uniform_real_distribution<double> dTheta(0.0,M_PI/2.0);
    
    theta = dTheta(gen);
    phi = dPhi(gen);
    
    c1.setNewAxis(theta, phi);
    
    double theta2 = 0.0, phi2 = 0.0;
    c1.getAxisAngles(theta2, phi2);

    REQUIRE( theta2 == Approx(theta) );
    REQUIRE( phi2 == Approx(phi) );
}


} // end of namspace FEvoR
