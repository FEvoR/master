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
    double b = 4.5-10; // units: m
    double alpha = 1; // units: -
        // Thor. 2002: constant grater than 1. However, everyone just uses 1:
            // Thor 2002, De La Chapelle 1998, Montagnant 2000
    double rhoDot = 0; // units: m^{-2} s^{-1}
    
    
    // Change in disloation density
    rhoDot = Medot/(b*cSize) - alpha*cDislDens*K/(cSize*cSize);  // units: m^{-2} s^{-1}
    
    std::cout << "\n rhoDot: " << rhoDot << "\n" << std::endl;
    
    cDislDens = cDislDens + rhoDot*timeStep; // units: m^{-2}
    
    // set boundary condition -- can't have a negative dislocation density
    if (cDislDens < 0)
        cDislDens = 0; // units: m^{-2}
}

void fevor_crystal::migRe(const double &Estress, const double &modelTime, const double &timeStep) {
    double Ggb = 0.065; // units: J m^{-2}
    double G = 3.4e9; // units: Pa
    double b = 4.5e-10; //units: m
    double kappa = 0.35;
    /* adjustable parameter -- Thor. 2002 eqn. 19 -- value set [38] 
     * Thor. 2002 sets this to 0.35, which aligns with paper he cites 
     * (> 0.1; Mohamed2000 -- a paper based on cold rolling copper to .35 
     * strain), however, this does not give the correct results -- it will 
     * cause an undeformed crystal with the minimum dislocation density to
     * recrystallize. From Thor. 2002, initially crystals start with:
     *   dislDens: 4e10 % m^{-2}
     *   size:     4 mm
     * which gives:
     *   Egb  = 48.75 J m^{-3}
     *   Edis = 89.79 J m^{-3}
     * then in paragraph [38] he states " initially..., so Edis is very
     * small relative to the grain boundary energy." Which is not at all the
     * case using his equations. 
     *
     * It looks like Edis should include the log term as well --
     * log(1/(std::sqrt(cDislDens).*b)) is equal to about 10. Removing the
     * log term you would end up with 9.63 J m^{-3} which is actually much
     * less than 48.75 J m^{-3}. So, keeping the log term in the equation,
     * you should use kappa = 0.035 (?). Or, drop the log term out of the 
     * Edis equation below  and keep kappa at 0.35. For computational
     * efficiency, I am going to drop the log term. It can be added in by
     * copying it from above. Also, Mohamed 2000 states that the log term
     * can be a constant. Note: log() is the natural log.
     */
    
    double Egb = 0, Edis = 0;
    Egb = 3*Ggb/cSize; // units: J m^{-3}
    Edis = kappa*G*cDislDens*b*b; // units: J m^{-3}
    
    if (Edis <= Egb)
        return;

    cDislDens = 1e10; // units: m^{-2} 
    
    double PC = 1; // units Pa^{4/3} m 
        // Shimizu 2008
    cSize = PC*Estress; // units: m
    
    // TODO: select an orientation!    
    
    cTimeLastRecrystal = modelTime+timeStep;
    cSizeLastRecrystal = cSize;
    
}

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

