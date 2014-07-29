/*
 */

#include <vector>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <random>
#include "fevor_crystal.hh"
#include "vector_tensor_opperations.hh"

// Define function members

std::vector<double> fevor_crystal::resolveM(const double &temperature, const std::vector<double> &stress, double &Mrss) {
    double glenExp = 3;
    double R = 0.008314472; // units: kJ K^{-1} mol^{-1}
    double beta = 630; // from Thors 2001 paper (pg 510, above eqn 16)
    
    double A = 0;
    if (temperature > -10) {
        double Q = 115;
        A = 24e-25*beta*exp(-Q/(R*(273.13+temperature))); // units: s^{-1} Pa^{-n}
    } else {
        double Q = 60;
        A = 3.5e-25*beta*exp(-Q/(R*(273.13+temperature))); // units: s^{-1} Pa^{-n}
    }
    // From Cuffy + Patterson (4 ed.) pg. 73
    
    // Burgers vector for each slip system
    std::vector<double> burger1, burger2, burger3;
    if ((cAxis[1] == 0) && (cAxis[2] == 0)) {
        burger1 = {1/3,0,0};
        burger2 = {(1+sqrt(3))/6,-sqrt(3)/6,0};
        burger3 = {(1-sqrt(3))/6,sqrt(3)/6,0};
    } else {
        double xyline = sqrt(cAxis[0]*cAxis[0]+cAxis[1]*cAxis[1]);
        
        burger1 = {cAxis[0]*cAxis[2]/xyline/3,
                   cAxis[1]*cAxis[2]/xyline/3,
                   -xyline/3};
        
        burger2 = {burger1[0]/2+sqrt(3)*cAxis[1]/xyline/6,
                   burger1[1]/2-sqrt(3)*cAxis[0]/xyline/6,
                   -xyline/6};
        
        burger3 = {burger1[0]/2-sqrt(3)*cAxis[1]/xyline/6,
                   burger1[1]/2+sqrt(3)*cAxis[0]/xyline/6,
                   -xyline/6};
    }
    
    // calculate shmidt tensors
        // 1x6 vector containing the 6 independent elements of the 3x3 
        // shmidt tensor (upper triangle) in ROW-MAJOR order. 
    std::vector<double> shmidt1, shmidt2,shmidt3;
    shmidt1 = vectorOuter(burger1,cAxis);
    shmidt2 = vectorOuter(burger2,cAxis);
    shmidt3 = vectorOuter(burger3,cAxis);
    
    std::vector<double> Mbase1, Mbase2, Mbase3;
    Mbase1 = tensorOuter(shmidt1, shmidt1);
    Mbase2 = tensorOuter(shmidt2, shmidt2);
    Mbase3 = tensorOuter(shmidt3, shmidt3);
    
    double rss1, rss2,rss3;
    rss1 = std::inner_product(shmidt1.cbegin(),shmidt1.cend(), stress.cbegin(), 0);
    rss2 = std::inner_product(shmidt2.cbegin(),shmidt2.cend(), stress.cbegin(), 0);
    rss3 = std::inner_product(shmidt3.cbegin(),shmidt3.cend(), stress.cbegin(), 0);
    
    std::transform(burger1.begin(),burger1.end(),burger1.begin(), 
                    [&](double x){return x*rss1;});
    std::transform(burger2.begin(),burger2.end(),burger2.begin(), 
                    [&](double x){return x*rss2;});
    std::transform(burger3.begin(),burger3.end(),burger3.begin(), 
                    [&](double x){return x*rss3;});
    
    std::transform(burger1.begin(),burger1.end(),burger2.begin(), burger1.begin(),
                    std::plus<double>());
    std::transform(burger1.begin(),burger1.end(),burger3.begin(), burger1.begin(),
                    std::plus<double>());
    
    Mrss = tensorMagnitude(burger1);
    
    
    std::transform(Mbase1.begin(), Mbase1.end(), Mbase1.begin(), 
                   [&](double x){return A*pow(rss1,glenExp-1)*x;} );
    std::transform(Mbase2.begin(), Mbase2.end(), Mbase2.begin(), 
                   [&](double x){return A*pow(rss2,glenExp-1)*x;} );
    std::transform(Mbase3.begin(), Mbase3.end(), Mbase3.begin(), 
                   [&](double x){return A*pow(rss3,glenExp-1)*x;} );
                   
    
    std::transform(Mbase1.begin(), Mbase1.end(), Mbase2.begin(),Mbase1.begin(), 
                   std::plus<double>());
    std::transform(Mbase1.begin(), Mbase1.end(), Mbase3.begin(),Mbase1.begin(), 
                   std::plus<double>());
    
    return Mbase1;
    
}

double fevor_crystal::grow(const double &tempature, const double &modelTime) {
    double K_0 = 8.2e-9; // units: m^2 s^{-1}
    double R = 0.008314472; // units: kJ K^{-1} mol^{-1}
    double Q = 0;
    if (tempature >= -10) // in degrees C
        Q = 0.7*115; // units: kJ mol^{-1}
    else
        Q = 0.7*60; // units: kJ mol^{-1}
        // From Cuffy + Patterson (4 ed.) pg. 40
    
    double K = 0;
    K  = K_0*exp(-Q/(R*(273.13+tempature))); // units: m^2 s^{-1}
    
    cSize = std::sqrt(K*(modelTime-cTimeLastRecrystal) + cSizeLastRecrystal*cSizeLastRecrystal);
    
    return K;

}


void fevor_crystal::dislocate(const double &timeStep, const double &Medot, const double &K) {
    double b = 4.5e-10; // units: m
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

unsigned int fevor_crystal::migRe(const std::vector<double> &stress, const double &modelTime, const double &timeStep) {
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
        return 0;

    cDislDens = 1e10; // units: m^{-2} 
    
    double Mstress = 0;
    std::vector<double> c = {1, sqrt(2), sqrt(2), 1, sqrt(2), 1};
    std::transform(c.begin(), c.end(), stress.begin(),c.begin(), 
                   std::multiplies<double>());
    
    Mstress = tensorMagnitude(c);
    
    double PC = 1; // units Pa^{4/3} m 
        // Shimizu 2008
    cSize = PC*Mstress; // units: m
    
    // Select an orientation!
        // Should be close to max MRSS
    double theta = 0, phi = 0;
    getAxisAngles(theta, phi);
    const double PI  =3.141592653589793238463;
    
    int stressIndex1 = (stress[3] > stress[0] ? 3 : 0);
    int stressIndex2 = (stress[1] > stress[2] ? 1 : 2);
    stressIndex2     = (stress[4] > stress[stressIndex2] ? 4 : stressIndex2);
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> dPhi(0,2*PI);
    
     if (stress[stressIndex1] < stress[stressIndex2]) {
        // simple shear -- orientation will be near vertical
        std::uniform_real_distribution<double> dTheta(0,PI/6);
        theta = dTheta(generator);
        phi = dPhi(generator);
        
        
    } else {
        // uniaxial comp. or pure shear -- orientation will be near theta = 45 degrees
        std::uniform_real_distribution<double> dTheta(PI/6,PI/3);
        theta = dTheta(generator);
        phi = dPhi(generator);
    }
    
    getNewAxis(theta, phi);
    
    cTimeLastRecrystal = modelTime+timeStep;
    cSizeLastRecrystal = cSize;
    
    return 1;
}

unsigned int fevor_crystal::polygonize( const std::vector<double> &stress, const double &Mrss, const double &modelTime, const double &timeStep) {
    double del = 0.065; // units: - 
        // ratio threshold -- Thor. 2002 [26]
    double rhop = 5.4e10; // units: m^{-2} 
        // dislocation density needed to form a wall -- Thor. 2002 [26]
    double Mstress = 0;
    std::vector<double> c = {1, sqrt(2), sqrt(2), 1, sqrt(2), 1};
    std::transform(c.begin(), c.end(), stress.begin(),c.begin(), 
                   std::multiplies<double>());
    
    Mstress = tensorMagnitude(c);
    
    if (Mrss/Mstress >= del || cDislDens < rhop)
        return 0;
    
    cDislDens -= rhop;
    cSize /= 2;
    
    // Select an orientation!
        // Should be toward max MRSS --away from vertical in uni. comp. or 
        // pure shear, towards vertical if simple shear.
    double theta = 0, phi = 0;
    getAxisAngles(theta, phi);
    const double PI  =3.141592653589793238463;
    
    int stressIndex1 = (stress[3] > stress[0] ? 3 : 0);
    int stressIndex2 = (stress[1] > stress[2] ? 1 : 2);
    stressIndex2     = (stress[4] > stress[stressIndex2] ? 4 : stressIndex2);
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0,100);
    
    if (stress[stressIndex1] < stress[stressIndex2] && theta < PI/6) {
        // toward vertical
        theta -= PI/36;
        
    } else if (theta < PI/6) {
        // away from vetical
        theta += PI/36;
    } else {
        // randomly toward/away from vertical
        theta += (  distribution(generator) < 50 ? -PI/36 : PI/36);
    } 
    
    getNewAxis(theta, phi);
    
    cTimeLastRecrystal = modelTime+timeStep;
    cSizeLastRecrystal = cSize;
    
    return 1;
}



void fevor_crystal::seeCrystal() { 
    
    std::cout << "The C-Axis orientation is:" << std::endl;
    for (auto &i : cAxis)
        std::cout << i << " ";
    std::cout << std::endl;

    std::cout << "The crystal size is:" << std::endl;
    std::cout << cSize << " m" << std::endl;

    std::cout << "The crystal dislocation density is:" << std::endl;
    std::cout << cDislDens << " m^{-2}" << std::endl;
}

void fevor_crystal::getAxisAngles(double &theta, double &phi) {
    
    double HXY = sqrt(cAxis[0]*cAxis[0] + cAxis[1]*cAxis[1]);
    theta = atan2(HXY, cAxis[2]);
    phi   = atan2(cAxis[1], cAxis[0]);
    
}

void fevor_crystal::getNewAxis(double &theta, double &phi) {
    
    cAxis = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
    
}
