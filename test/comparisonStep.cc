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
#include <sstream>
#include <vector>
#include <ctime>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <random>
#include "fevor_distribution.hh"
#include "vector_tensor_operations.hh"


// Perform a compairison time step between FEvoR and Thor
int main(int argc, char *argv[])
{
    if (argc != 4) {
        std::cout << "  usage: " << argv[0] << " w T(K) dt(years)" << std::endl;
        std::cout << "example: " << argv[0] << "-3.0 263.15 1.0" << std::endl;
        return 0;
    }
    
    
    double watsonK, temperature, dt;
    std::stringstream ss;
    ss << argv[1]; ss >> watsonK;
    ss.clear();
    ss << argv[2]; ss >> temperature; // Kelvin
    ss.clear();
    ss << argv[3]; ss >> dt; // years
    ss.clear();
    
    double modelTime = 0.0;
    double timeStep = dt*365.0*24.0*60.0*60.0;
    
    // stress
    std::vector<double> stress = { 10000,     0, 10000,
                                       0,     0,     0,
                                   10000,     0,-10000};
    
    
    // Diagnostics
    unsigned int nMigre, nPoly;
    nMigre = nPoly = 0;
    std::vector<double> bulkEdot(9, 0.0);
    
    // softness parameters
    double cc = 1.0, cn = 0.0;
    std::vector<unsigned int> packingDimensions = {20,20,20};
    
    FEvoR::Distribution d1(packingDimensions, watsonK);
    d1.setSoftnessRatio(cc, cn);
    
    // calculate the flow factor
    double glenExp = 3.0;
    double R = 0.008314472; // units: kJ K^{-1} mol^{-1}
    double beta = 630.0; // from Thors 2001 paper (pg 510, above eqn 16)
    double Q = 0.0;
    double A = 0.0;
    
    Q = (temperature > 263.15 ? 115.0 : 60.0);
    A = 3.5e-25*exp(-(Q/R)*(1.0/temperature - 1.0/263.15));
    
    
    std::cout   << "\n" 
                << "************************************ \n"
                << "Compairison timeStep: FEvoR vs Thor. \n"
                << "************************************ \n"
                << "\n"
                << "Model Setup is: \n"
                << "    " << "Softness Parameters: " << cc << " , " << cn << " \n"
                << "    " << "Temperature: " << temperature << " Kelvin\n"
                << "    " << "Flow factor: " << A << " s^-1 Pa^-3\n"
                << "    " << "Stress (Pa): " << std::endl;
    FEvoR::tensorDisplay(stress, 3, 3);
    
    std::cout   << "\n" 
                << "Initial time: " << modelTime << " s \n"
                << "   Time Step: " << timeStep  << " s \n"
                << std::endl;
    
    double theta = M_PI/3,
             phi = M_PI/3; 
    
    FEvoR::Crystal c1({0.0,0.0,1.0}, 0.01, 1.0e10);
    c1.setNewAxis(theta, phi);
    
    double  Mrss = 1.0,
            Medot =0.0;
    std::vector<double> cM = c1.resolveM(temperature, stress, Mrss, Medot);
    
    // Burgers vector for each slip system
    std::vector<double> burger1, burger2, burger3;
    burger1 = burger2 = burger3 = {0.0,0.0,0.0};
        
        // sines and cosines so calculation only has to be preformed once
        double  st = std::sin(theta), 
                ct = std::cos(theta),
                sp = std::sin(phi),
                cp = std::cos(phi), 
                sq3= sqrt(3.0);

    burger1 = {ct*cp/3.0,             ct*sp/3.0,          -st/3.0};
    burger2 = {(-ct*cp - sq3*sp)/6.0, (-ct*sp+sq3*cp)/6.0, st/6.0};
    burger3 = {(-ct*cp + sq3*sp)/6.0, (-ct*sp-sq3*cp)/6.0, st/6.0};
    
    // calculate shmidt tensors
        // 1x9 vector containing the  3x3 shmidt tensor in ROW-MAJOR order. 
    std::vector<double> cAxis = {st*cp, st*sp, ct};
    std::vector<double> shmidt1, shmidt2,shmidt3;
    shmidt1 = FEvoR::tensorOuter(burger1,cAxis);
    shmidt2 = FEvoR::tensorOuter(burger2,cAxis);
    shmidt3 = FEvoR::tensorOuter(burger3,cAxis);
    
    double rss1, rss2, rss3;
    rss1 = std::inner_product(shmidt1.cbegin(),shmidt1.cend(), stress.cbegin(), 0.0);
    rss2 = std::inner_product(shmidt2.cbegin(),shmidt2.cend(), stress.cbegin(), 0.0);
    rss3 = std::inner_product(shmidt3.cbegin(),shmidt3.cend(), stress.cbegin(), 0.0);
    
    std::vector<double> vel;
    vel = FEvoR::tensorMixedInner(cM, stress);
    
    
    std::vector<double> velT, edot;
    edot.resize(vel.size());
    velT = FEvoR::matrixTranspose(vel, 3, 3);
    
    std::transform(vel.begin(), vel.end(), velT.begin(),edot.begin(), 
                   std::plus<double>());
    std::transform(edot.begin(),edot.end(),edot.begin(), 
                    [&](double x){return x/2.0;});
    
    std::cout   << "************************************ \n"
                << "   Stepping a crystal \n"
                << "************************************ \n"
                << "  shmidt1:" << std::endl;
    FEvoR::tensorDisplay(shmidt1,3,3);
    std::cout   << "  shmidt2:" << std::endl;
    FEvoR::tensorDisplay(shmidt2,3,3);
    std::cout   << "  shmidt3:" << std::endl;
    FEvoR::tensorDisplay(shmidt3,3,3);
    std::cout   << "     rss1:" << rss1 << "\n"
                << "     rss2:" << rss2 << "\n"
                << "     rss3:" << rss3 << "\n"
                << "     Mrss:" << Mrss << "\n"
                << "Crystal M:" << std::endl;
    FEvoR::tensorDisplay(cM, 9,9);
    std::cout   << "      vel:" << std::endl;
    FEvoR::tensorDisplay(vel, 3,3);
    std::cout   << "     velT:" << std::endl;
    FEvoR::tensorDisplay(velT, 3,3);
    std::cout   << "     edot:" << std::endl;
    FEvoR::tensorDisplay(edot, 3,3);
    std::cout   << "    Medot:" << Medot << std::endl;
    
    std::cout   << "************************************ \n"
                << "   Stepping a distribution \n"
                << "************************************ \n"
                << std::endl;
    
    d1.saveDistribution("comparisonDist_initial.csv");
    
    std::vector<double> bulkM;
    bulkM = d1.stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
    
    d1.saveDistribution("comparisonDist_final.csv");
    
    std::cout   << "   nMigRe:" << nMigre << "\n"
                << "    nPoly:" << nPoly << "\n"
                << "  dist. M:" << std::endl;
    FEvoR::tensorDisplay(bulkM, 9,9);
    std::cout   << "    bEdot:" << std::endl;
    FEvoR::tensorDisplay(bulkEdot, 3,3);
    
    return 0;
}
