/*
 */

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <random>
#include "fevor_crystal.hh"
#include "vector_tensor_opperations.hh"
#include "fevor_distribution.hh"

// Define function members
std::vector<double> fevor_distribution::stepInTime(const double &temperature, const std::vector<double> &stress, double &modelTime, const double &timeStep, double &nMigre, double &nPoly, std::vector<double> &bulkEdot) {
    
    double crystalMagEdot;
    std::vector<double> bulkM(81, 0.0);
    std::vector<std::vector<double>> crystalM;
    double crystalK;
    
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        crystalM.push_back( crystals[ii].resolveM(temperature, stress, magRSS[ii], crystalMagEdot) );
        
        crystalK = crystals[ii].grow(temperature, modelTime);

        crystals[ii].dislocate(timeStep, crystalMagEdot, crystalK);

        nMigre += crystals[ii].migRe(stress, modelTime, timeStep);

        nPoly  += crystals[ii].polygonize(stress, magRSS[ii], modelTime, timeStep);
        
    }
    
    
    getSoftness(crystalM, bulkM, bulkEdot, stress);
    
    for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
        
        crystals[ii].rotate(crystalM[ii], bulkEdot, stress);
    }
    
    modelTime += timeStep;
    
    return bulkM;
}

void fevor_distribution::getSoftness(std::vector<std::vector<double>> &crystalM, std::vector<double> &bulkM, std::vector<double> &bulkEdot, const std::vector<double> &stress) {
    if (contribNeighbor != 0.0) {
        
        unsigned int front, back, left, right, top, bottom;
        front = back = left = right = top = bottom = 0;
        double softy;
        
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
            front = (ii + numberCrystals + 1) % numberCrystals; 
            back  = (ii + numberCrystals - 1) % numberCrystals; 
            left  = (ii + numberCrystals + dimensions[0]) % numberCrystals;
            right = (ii + numberCrystals - dimensions[0]) % numberCrystals;
            top   = right = (ii + numberCrystals + dimensions[0]*dimensions[1]) % numberCrystals;
            bottom= right = (ii + numberCrystals - dimensions[0]*dimensions[1]) % numberCrystals;
            
            softy = (magRSS[front] + magRSS[back] + magRSS[left] + magRSS[right] + magRSS[top] + magRSS[bottom])/magRSS[ii];
            
            softness[ii] = 1.0/(contribCrystal + 6*contribNeighbor)*(contribCrystal + contribNeighbor*softy);
            
            std::transform(crystalM[ii].begin(),crystalM[ii].end(),crystalM[ii].begin(), 
                        [&](double x){return x*softness[ii];});
                        
            std::transform(bulkM.begin(),bulkM.end(),crystalM[ii].begin(), bulkM.begin(),
                           std::plus<double>());
        }
    } else {
        for(unsigned int ii = 0; ii!= numberCrystals; ++ii) {
            
            std::transform(bulkM.begin(),bulkM.end(),crystalM[ii].begin(), bulkM.begin(),
                           std::plus<double>());
            
        }
    }
    
    std::vector<double> bulkMtrans;
    bulkMtrans = matrixTranspose(bulkM,9,9);
    std::transform(bulkM.begin(), bulkM.end(), bulkMtrans.begin(),bulkMtrans.begin(), 
                   std::plus<double>());
    std::transform(bulkMtrans.begin(),bulkMtrans.end(),bulkMtrans.begin(), 
                    [&](double x){return x/(2.0*numberCrystals);});
                    
    bulkEdot = tensorMixedInner(bulkMtrans, stress);
}

void fevor_distribution::setSoftness(double cc, double cn) {
    contribCrystal  = cc;
    contribNeighbor = cn;
}

void fevor_distribution::saveDistribution() {
    // FIXME: Make it actually save, now just prints to cout
    
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

void fevor_distribution::generateWatsonAxes(const double &wk) {
    const double PI  =3.141592653589793238463;
    //TODO: look up random_device vs default_random_engine and how to preserve
    // it across function calls
    std::random_device seed;
    
    if (wk == 0) {
        // isotropic
        std::normal_distribution<double> rNumb(-1.0,1.0);
        std::vector<double> axis;
        
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) { 
        
            axis = {rNumb(seed), rNumb(seed), rNumb(seed)};
            double axisMag = tensorMagnitude(axis);
            if (axisMag != 1.0) {
                std::transform(axis.begin(), axis.end(), axis.begin(), 
                                [&](double x){return x/sqrt(axisMag);} );
            }
        
            crystals[ii].getNewAxis(axis);
        }
        
    } else if (wk < -1000.5) {
        // perfect bipolar (single maximum)
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) { 
                crystals[ii].getNewAxis({0.0,0.0,1.0});
        }
    } else if (wk < 0) {
        // bipolar (single maximum)
        //TODO: find c++ equivelent of randomsample() in MATLAB
        
    } else if (wk > 500.0) {
        // perfect girdle
        std::uniform_real_distribution<double> gPhi(0,2.0*PI);
        double theta, phi;
        theta = PI/2.0;
        
        for (unsigned int ii = 0; ii!= numberCrystals; ++ii) { 
                phi   = gPhi(seed);
                crystals[ii].getNewAxis(theta, phi);
        }
    } else { //(wk > 0)
        // girdle
        //TODO: find c++ equivelent of erf() and erfinv() in MATLAB
            // or better way to calculate Q
    }
}
