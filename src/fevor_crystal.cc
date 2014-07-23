/*
 */

#include <vector>
#include <iostream>
#include <cmath>
#include "fevor_crystal.hh"

// Define function members

double fevor_crystal::grow(const double &tempature, const double &modelTime) {
    double K_0 = 8.2e-9; // units: m^2 s^{-1}
    double R = 0.008314472; // units: kJ K^{-1} mol^{-1}
    double Q = 0;
    if (tempature >= -10) // in degrees C
        Q = 0.7*115; // units: kJ mol^{-1}
    else
        Q = 0.7*60; // units: kJ mol^{-1}
        // From Kuffy + Patterson (4 ed.) pg. 40
    
    double K = 0;
    K  = K_0*exp(-Q/(R*(273.13+tempature))); // units: m^2 s^{-1}
    
    cSize = std::sqrt(K*(modelTime-cTimeLastRecrystal) + cSizeLastRecrystal*cSizeLastRecrystal);
    
    return K;

}


void fevor_crystal::dislocate(const double &timeStep, const double &Medot, const double &K) {
    double b = 3.69e-10; // units: m
    double alpha = 1; // units: -
        // Thor. 2002: constant grater than 1. However, everyone just uses 1:
            // Thor 2002, De La Chapelle 1998, Montagnant 2000
    double rhoDot = 0; // units: m^{-2} s^{-1}
    
    
    // Change in disloation density
    rhoDot = Medot/(b*cSize) - alpha*cDislDens*K/(cSize*cSize);  // units: m^{-2} s^{-1}
    
    cDislDens = cDislDens + rhoDot*timeStep; // units: m^{-2}
    
    // set boundary condition -- can't have a negative dislocation density
    if (cDislDens < 0)
        cDislDens = 0; // units: m^{-2}
}

//~ void fevor_crystal::MigRe();

//~ void fevor_crystal::polygonize();



void fevor_crystal::seeCrystal() { 
    
    std::cout << "The C-Axis orientation is:" << std::endl;
    for (auto i : cAxis)
        std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "The crystal size is:" << std::endl;
    std::cout << cSize << " m" << std::endl;

    std::cout << "The crystal dislocation density is:" << std::endl;
    std::cout << cDislDens << " m^{-2}" << std::endl;
}

// Define static data members

