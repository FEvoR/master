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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <random>
#include <cassert>

#include "fevor_crystal.hh"
#include "vector_tensor_operations.hh"
#include "fevor_distribution.hh"
#include "Faddeeva.hh"

namespace FEvoR {
// Each crystal has: the three components of the crystal axis, size,
// dislocation density, time of last recrystallization, and size at
// last recrystallization.
const unsigned int Distribution::numberParameters = 7;

Distribution::Distribution(std::vector<unsigned int> lwh): dimensions(lwh) {
    numberCrystals = dimensions[0]*dimensions[1]*dimensions[2];
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        crystals.push_back( Crystal ({0.0,0.0,1.0}, 0.01, 1.0e10) );
        softness.push_back(1.0);
        magRSS.push_back(1.0);
    }
    
    setSoftnessRatio(6.0, 1.0);
}

// construct a distribution from a saved distribution of crystals
Distribution::Distribution(std::vector<unsigned int> lwh, std::string fname): Distribution(lwh)  {
    loadDistribution(fname);
}

// construct a distribution using the Watson distribution for axis angles
Distribution::Distribution(std::vector<unsigned int> lwh, double wk): Distribution(lwh)  {
    generateWatsonAxes(wk);
}

// construct a distribution from a big vector of all distribution data
Distribution::Distribution(std::vector<unsigned int> lwh, std::vector<double> &data): dimensions(lwh) {
    numberCrystals = dimensions[0]*dimensions[1]*dimensions[2];
    const int n = numberParameters;
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        
        crystals.push_back( Crystal ({data[ii*n+0],data[ii*n+1],data[ii*n+2]},
                                           data[ii*n+3], data[ii*n+4],
                                           data[ii*n+5], data[ii*n+6]) );
        softness.push_back(1.0);
        magRSS.push_back(1.0);
    }
    
    setSoftnessRatio(6.0, 1.0);
}

// Define function members
std::vector<double> Distribution::stepInTime(const double &temperature, const std::vector<double> &stress, const double &modelTime, const double &timeStep, unsigned int &nMigre, unsigned int &nPoly, std::vector<double> &bulkEdot){
    
    std::vector<double> crystalMagEdot(numberCrystals, 0.0);
    std::vector<double> bulkM(81, 0.0);
    std::vector<std::vector<double>> crystalM;
    double crystalK=0.0;

    /* Ensure deviatoric stress 
     * Indexing: {0, 1, 2,
     *            3, 4, 5,
     *            6, 7, 8}
     */
    double hydrostatic = (stress[0] + stress[4] + stress[8])/3.0;
    std::vector<double> devStress = stress;
    devStress[0] -= hydrostatic;
    devStress[4] -= hydrostatic;
    devStress[8] -= hydrostatic;
    
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        crystalM.push_back( crystals[ii].resolveM(temperature, devStress, magRSS[ii], crystalMagEdot[ii]) );
    }
    
    getSoftness(crystalM, crystalMagEdot, bulkM, bulkEdot, devStress);
    
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {    
        crystalK = crystals[ii].grow(temperature, modelTime);

        crystals[ii].dislocate(timeStep, crystalMagEdot[ii], crystalK);

        nMigre += crystals[ii].migRe(devStress, modelTime, timeStep);

        nPoly  += crystals[ii].polygonize(devStress, magRSS[ii], modelTime, timeStep);
        
    }
    
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        
        crystals[ii].rotate(crystalM[ii], bulkEdot, devStress, timeStep);
    }
    
    return bulkM;
}

std::vector<double> Distribution::stepInTime(const double &temperature, const std::vector<double> &stress, const double &modelTime, const double &timeStep, std::vector<double> &bulkEdot){
    unsigned int nMigre = 0,
                 nPoly  = 0;
 
    return stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
}

void Distribution::getSoftness(std::vector<std::vector<double> > &crystalM, std::vector<double> &crystalMagEdot, std::vector<double> &bulkM, std::vector<double> &bulkEdot, const std::vector<double> &stress) {
    
    /* Ensure deviatoric stress 
     * Indexing: {0, 1, 2,
     *            3, 4, 5,
     *            6, 7, 8}
     */
    double hydrostatic = (stress[0] + stress[4] + stress[8])/3.0;
    std::vector<double> devStress = stress;
    devStress[0] -= hydrostatic;
    devStress[4] -= hydrostatic;
    devStress[8] -= hydrostatic;
    
    if (contribNeighbor != 0.0) {
        double glenExp = 3.0;
        
        unsigned int front, back, left, right, top, bottom;
        front = back = left = right = top = bottom = 0;
        double softy = 10.0;
        
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
            front = (ii + numberCrystals + 1) % numberCrystals; 
            back  = (ii + numberCrystals - 1) % numberCrystals; 
            left  = (ii + numberCrystals + dimensions[0]) % numberCrystals;
            right = (ii + numberCrystals - dimensions[0]) % numberCrystals;
            top   = right = (ii + numberCrystals + dimensions[0]*dimensions[1]) % numberCrystals;
            bottom= right = (ii + numberCrystals - dimensions[0]*dimensions[1]) % numberCrystals;
            
            if (magRSS[ii] != 0.0)
                softy = (magRSS[front] + magRSS[back] + magRSS[left] + magRSS[right] + magRSS[top] + magRSS[bottom])/magRSS[ii];
                
            // maximum softness
            if (softy > 10.0)
                softy = 10.0; 
            
            assert(contribCrystal != 0.0 || contribNeighbor !=0.0);
            softness[ii] = 1.0/(contribCrystal + 6.0*contribNeighbor)*(contribCrystal + contribNeighbor*softy);
            
            std::transform(crystalM[ii].begin(),crystalM[ii].end(),crystalM[ii].begin(), 
                        [&](double x){return x*std::pow(softness[ii],glenExp);});
            
            crystalMagEdot[ii] *= std::pow(softness[ii],glenExp);
                        
            std::transform(bulkM.begin(),bulkM.end(),crystalM[ii].begin(), bulkM.begin(),
                           std::plus<double>());
        }
    } else {
        for(unsigned int ii = 0; ii!= numberCrystals; ++ii) {
            
            std::transform(bulkM.begin(),bulkM.end(),crystalM[ii].begin(), bulkM.begin(),
                           std::plus<double>());
            
        }
    }
    
    assert(numberCrystals != 0);
    std::transform(bulkM.begin(),bulkM.end(),bulkM.begin(), 
                    [&](double x){return x/numberCrystals;});
    
    std::vector<double> bulkVel;
    bulkVel  = tensorMixedInner(bulkM, devStress);
    
    std::vector<double> bulkVelT;
    bulkVelT = matrixTranspose(bulkVel, 3, 3);
    
    std::transform(bulkVel.begin(), bulkVel.end(), bulkVelT.begin(),bulkEdot.begin(), 
                   std::plus<double>());
    std::transform(bulkEdot.begin(),bulkEdot.end(),bulkEdot.begin(), 
                    [&](double x){return x/2.0;});
}

void Distribution::setSoftnessRatio(double cc, double cn) {
    contribCrystal  = cc;
    contribNeighbor = cn;
}

void Distribution::saveDistribution() const {
    std::cout << "# "
              << "Crystal"                << ", "
              << "C-Axis (x)"             << ", "
              << "C-Axis (y)"             << ", "
              << "C-Axis (z)"             << ", "
              << "Size (m)"               << ", "
              << "Disl. dens. (1/m^2)"    << ", "
              << "Last recr. time (s)"    << ", "
              << "Size at last recr. (m)" << std::endl;
    
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        
        std::cout << std::setw(5) << ii << ", ";
                  
        crystals[ii].printCrystal();
    }
}
void Distribution::saveDistribution(std::string fname) const {
    std::ofstream file(fname);
    
    file << "# "
              << "Crystal"                << ", "
              << "C-Axis (x)"             << ", "
              << "C-Axis (y)"             << ", "
              << "C-Axis (z)"             << ", "
              << "Size (m)"               << ", "
              << "Disl. dens. (1/m^2)"    << ", "
              << "Last recr. time (s)"    << ", "
              << "Size at last recr. (m)" << std::endl;
    
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        
        file << std::setw(5) << ii << ", ";
                  
        crystals[ii].printCrystal(file);
    }
}
void Distribution::saveDistribution(std::vector<double> &data) const {
    data.resize(numberParameters*numberCrystals);
    const int n = numberParameters;
  
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        crystals[ii].getAll(data[ii*n+0], data[ii*n+1], data[ii*n+2],
                            data[ii*n+3], data[ii*n+4],
                            data[ii*n+5], data[ii*n+6]);
    }
}

void Distribution::loadDistribution( std::string fname ) {
    std::ifstream file(fname);
    std::string line;
    std::string field;
    double fieldVal;
    unsigned int ii = 0;
    std::vector<double> data;
    
    while (std::getline(file, line) && ii != numberCrystals) {
        if (line.find("#") == 0)
            continue;
            
        std::istringstream streamLine(line);
        while (std::getline(streamLine,field,',')) {
            std::istringstream streamField(field);
            
            streamField >> fieldVal;
            data.push_back(fieldVal);
        }
        
        ++ii; 
    }
    
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        // data[ii*8] is crystal number, skip this
        crystals[ii].setAll(data[ii*8+1],
                            data[ii*8+2],
                            data[ii*8+3],
                            data[ii*8+4],
                            data[ii*8+5],
                            data[ii*8+6],
                            data[ii*8+7]);
    }
}
void Distribution::loadDistribution( const std::vector<double> &data ) {
    // TODO: check size!
    // assert(data.size() == numberCrystals*7); 
  const int n = numberParameters;
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        crystals[ii].setAll(data[ii*n+0], data[ii*n+1], data[ii*n+2],
                            data[ii*n+3], data[ii*n+4], 
                            data[ii*n+5], data[ii*n+6]);
    }
}

void Distribution::generateWatsonAxes(const double &wk) {

    //std::random_device seed;
    //std::mt19937 gen(seed());
    double seed = 1337.0;
    std::mt19937 gen(seed);
    
    if (wk == 0) {
        // isotropic
        std::normal_distribution<double> rNumb(0.0,1.0);
        std::vector<double> axis;
        
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) { 
        
            axis = {rNumb(gen), rNumb(gen), rNumb(gen)};
            
            
            double axisMag = tensorMagnitude(axis);
            assert(axisMag != 0.0);
            if (axisMag != 1.0) {
                std::transform(axis.begin(), axis.end(), axis.begin(), 
                                [&](double x){return x/axisMag;} );
            }
        
            crystals[ii].setNewAxis(axis);
        }
        
    } else if (wk < -700.0) {
        // perfect bipolar (single maximum)
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) { 
                crystals[ii].setNewAxis({0.0,0.0,1.0});
        }
    } else if (wk < 0.0) {
        // bipolar (single maximum)
        std::vector<double> weights;
        std::vector<double> x;
        double abs_wk = std::abs(wk);
        double DawsonVal = Faddeeva::Dawson(abs_wk);
        
        assert(DawsonVal != 0.0);
        for (double ww = 0.0; ww < M_PI; ww+=0.001){
            x.push_back(ww);
            weights.push_back( 1.0/(4.0*M_PI)*std::sqrt(abs_wk)*std::exp(wk)/DawsonVal*std::exp(-wk*std::cos(ww)*std::cos(ww)) );
        }
        
        std::uniform_real_distribution<double> gPhi(0,2.0*M_PI);
        std::piecewise_linear_distribution<double> dTheta (x.begin(),x.end(),weights.begin());
        double theta, phi;
        
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) { 
                phi   = gPhi(gen);
                theta = dTheta(gen);
                
                crystals[ii].setNewAxis(theta, phi);
        }
        
        
    } else if (wk > 500.0) {
        // perfect equatorial girdle
        std::uniform_real_distribution<double> gPhi(0,2.0*M_PI);
        double theta, phi;
        theta = M_PI/2.0;
        
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) { 
                phi   = gPhi(gen);
                crystals[ii].setNewAxis(theta, phi);
        }
    } else { //(wk > 0.0)
        // equatorial girdle
        std::vector<double> weights;
        std::vector<double> x;
        double erfVal = Faddeeva::erf(wk);
        
        assert(erfVal != 0.0);
        for (double ww = 0.0; ww < M_PI; ww+=0.001){
            x.push_back(ww);
            weights.push_back( 1.0/(2.0*M_PI)*std::sqrt(wk/M_PI)*(1.0/erfVal)*std::exp(-wk*std::cos(ww)*std::cos(ww)) );
        }
        
        std::uniform_real_distribution<double> gPhi(0,2.0*M_PI);
        std::piecewise_linear_distribution<double> dTheta (x.begin(),x.end(),weights.begin());
        double theta, phi;
        
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) { 
                phi   = gPhi(gen);
                theta = dTheta(gen);
                
                crystals[ii].setNewAxis(theta, phi);
        }
    }
}

unsigned int Distribution::getNumberCrystals() const {
    return numberCrystals;
}

} // end of namespace FEvoR
