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

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <random>
#include <cassert>

#include "fevor_crystal.hh"
#include "vector_tensor_operations.hh"

namespace FEvoR {

Crystal::Crystal(std::vector<double> ca, double cs, double cdd):
    cAxis(ca), cSize(cs), cDislDens(cdd) {
    cTimeLastRecrystal = 0.0;
    cSizeLastRecrystal = cSize;
}
Crystal::Crystal(std::vector<double> ca, double cs, double cdd, double cto, double cso):
    cAxis(ca), cSize(cs), cDislDens(cdd), cTimeLastRecrystal(cto), cSizeLastRecrystal(cso) {
        /* nothing to do */
}
    
// Define function members

std::vector<double> Crystal::resolveM(const double &temperature, const std::vector<double> &stress, double &Mrss, double &Medot) {
    double glenExp = 3.0;
    double R = 0.008314472; // units: kJ K^{-1} mol^{-1}
    double beta = 630.0; // from Thors 2001 paper (pg 510, above eqn 16)
    double Q = 0.0;
    double A = 0.0;
    
    assert(temperature != 0.0);
    Q = (temperature > 263.15 ? 115.0 : 60.0);
    A = 3.5e-25*beta*std::exp(-(Q/R)*(1.0/temperature - 1.0/263.15)); // units: s^{-1} Pa^{-n}
    // From Cuffy + Patterson (4 ed.) pg. 73
    
    // Burgers vector for each slip system
    std::vector<double> burger1, burger2, burger3;
    burger1 = burger2 = burger3 = {0.0,0.0,0.0};
    
        // calculate using axis angles
        double theta, phi;
        angles(theta, phi);
        
        // sines and cosines so calculation only has to be preformed once
        double  st = std::sin(theta), 
                ct = std::cos(theta),
                sp = std::sin(phi),
                cp = std::cos(phi), 
                sq3= std::sqrt(3.0);

    burger1 = {ct*cp/3.0,             ct*sp/3.0,          -st/3.0};
    burger2 = {(-ct*cp - sq3*sp)/6.0, (-ct*sp+sq3*cp)/6.0, st/6.0};
    burger3 = {(-ct*cp + sq3*sp)/6.0, (-ct*sp-sq3*cp)/6.0, st/6.0};
    
    // calculate shmidt tensors
        // 1x9 vector containing the  3x3 shmidt tensor in ROW-MAJOR order. 
    std::vector<double> shmidt1, shmidt2,shmidt3;
    shmidt1 = tensorOuter(burger1,cAxis);
    shmidt2 = tensorOuter(burger2,cAxis);
    shmidt3 = tensorOuter(burger3,cAxis);
    
    // calculate M on each slip system
        // 1x81 vector containing the 9x9 matrix in ROW-MAJOR order. 
    std::vector<double> Mbase1, Mbase2, Mbase3;
    Mbase1 = tensorOuter(shmidt1, shmidt1);
    Mbase2 = tensorOuter(shmidt2, shmidt2);
    Mbase3 = tensorOuter(shmidt3, shmidt3);
    
    double rss1, rss2, rss3;
    rss1 = std::inner_product(shmidt1.cbegin(),shmidt1.cend(), stress.cbegin(), 0.0);
    rss2 = std::inner_product(shmidt2.cbegin(),shmidt2.cend(), stress.cbegin(), 0.0);
    rss3 = std::inner_product(shmidt3.cbegin(),shmidt3.cend(), stress.cbegin(), 0.0);
        
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
                   [&](double x){return A*std::pow(rss1,glenExp-1.0)*x;} );
    std::transform(Mbase2.begin(), Mbase2.end(), Mbase2.begin(), 
                   [&](double x){return A*std::pow(rss2,glenExp-1.0)*x;} );
    std::transform(Mbase3.begin(), Mbase3.end(), Mbase3.begin(), 
                   [&](double x){return A*std::pow(rss3,glenExp-1.0)*x;} );
           
    
    std::transform(Mbase1.begin(), Mbase1.end(), Mbase2.begin(),Mbase1.begin(), 
                   std::plus<double>());
    std::transform(Mbase1.begin(), Mbase1.end(), Mbase3.begin(),Mbase1.begin(), 
                   std::plus<double>());

    
    std::vector<double> vel;
    vel = FEvoR::tensorMixedInner(Mbase1, stress);
    
    std::vector<double> velT, edot;
    edot.resize(vel.size());
    velT = FEvoR::matrixTranspose(vel, 3, 3);
    
    std::transform(vel.begin(), vel.end(), velT.begin(),edot.begin(), 
                   std::plus<double>());
    std::transform(edot.begin(),edot.end(),edot.begin(), 
                    [&](double x){return x/2.0;});
                    
    
    Medot = tensorMagnitude(edot)*std::sqrt(1.0/2.0);
    
    return Mbase1;
    
}

double Crystal::grow(const double &temperature, const double &modelTime) {
    double K_0 = 8.2e-9; // units: m^2 s^{-1}
    double R = 0.008314472; // units: kJ K^{-1} mol^{-1}
    
    double Q = 0.7*60.0; // units: kJ mol^{-1}
    if (temperature >= 263.15) // units: Kelvin
        Q = 0.7*115.0; // units: kJ mol^{-1}
        // From Cuffy + Patterson (4 ed.) pg. 40
    
    double K = 0.0;
    assert(temperature != 0.0);
    K  = K_0*std::exp(-Q/(R*temperature)); // units: m^2 s^{-1}
    
    cSize = std::sqrt(K*(modelTime-cTimeLastRecrystal) + cSizeLastRecrystal*cSizeLastRecrystal);
    
    //~ std::cout << "Crystal K: " << K << std::endl;
    //~ std::cout << "Crystal Q: " << Q << std::endl;
    //~ std::cout << "Crystal K_0: " << K_0 << std::endl;
    //~ std::cout << "Temperature: " << temperature << std::endl;
    
    return K;

}

void Crystal::dislocate(const double &timeStep, const double &Medot, const double &K) {
    double b = 4.5e-10; // units: m
    double alpha = 1.0; // units: -
        // Thor. 2002: constant grater than 1. However, everyone just uses 1:
            // Thor 2002, De La Chapelle 1998, Montagnant 2000
    double rhoDot = 0.0; // units: m^{-2} s^{-1}
    
    
    assert(cSize != 0.0);
    // Change in dislocation density
    rhoDot = Medot/(b*cSize) - alpha*cDislDens*K/(cSize*cSize);  // units: m^{-2} s^{-1}
    
    cDislDens = cDislDens + rhoDot*timeStep; // units: m^{-2}
    
    // set boundary condition -- can't have a negative dislocation density
    if (cDislDens < 0.0)
        cDislDens = 0.0; // units: m^{-2}
}

unsigned int Crystal::migRe(const std::vector<double> &stress, const double &modelTime, const double &timeStep) {
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
    
    double Egb, Edis;
    assert(cSize != 0.0);
    Egb = 3.0*Ggb/cSize; // units: J m^{-3}
    Edis = kappa*G*cDislDens*b*b; // units: J m^{-3}
    
    if (Edis <= Egb)
        return 0;

    cDislDens = 1.0e10; // units: m^{-2} 
    
    double Mstress = 0.0;
    Mstress = tensorMagnitude(stress)*std::sqrt(1.0/2.0); 
        // effective stress: Thor. 2002, 3-4, [29]
    
    double PC = 1.0; // units Pa^{4/3} m 
        // Shimizu 2008
    cSize = PC*std::pow(Mstress, -4.0/3.0); // units: m
    
    // Select an orientation!
        // Should be close to max MRSS
    double theta = 0.0, phi = 0.0;
    angles(theta, phi);
    
    /* Indexing: {0, 1, 2,
     *            3, 4, 5,
     *            6, 7, 8}
     */
    int stressIndex1 = (std::abs(stress[4]) > std::abs(stress[8]) ? 4 : 8);
    int stressIndex2 = (std::abs(stress[1]) > std::abs(stress[2]) ? 1 : 2);
    stressIndex2     = (std::abs(stress[5]) > std::abs(stress[stressIndex2]) 
                        ? 5 : stressIndex2);
    
    std::random_device seed;
    std::mt19937 gen(seed());
    
    std::uniform_real_distribution<double> dPhi(0.0,2.0*M_PI);
    double u_pSix   = (1.0-std::cos(M_PI/6.0))/2.0,
           u_pThree = (1.0-std::cos(M_PI/3.0))/2.0;
    
    if ( std::abs(stress[stressIndex1]) < std::abs(stress[stressIndex2]) ) {
        // simple shear -- orientation will be near vertical
        std::uniform_real_distribution<double> du(0.0,u_pSix);
        double u   = du(gen);
        
        theta = std::acos(2.0*u-1.0); 
        phi   = dPhi(gen);
        
        
    } else {
        // uniaxial comp. or pure shear -- orientation will be near theta = 45 degrees
        std::uniform_real_distribution<double> du(u_pSix,u_pThree);
        double u   = du(gen);
        
        theta = std::acos(2.0*u-1.0);
        phi   = dPhi(gen);
    }
    
    (void) axis(theta, phi);
    
    cTimeLastRecrystal = modelTime+timeStep;
    cSizeLastRecrystal = cSize;
    
    return 1;
}

unsigned int Crystal::polygonize( const std::vector<double> &stress, const double &Mrss, const double &modelTime, const double &timeStep) {
    double del = 0.065; // units: - 
        // ratio threshold -- Thor. 2002 [26]
    double rhop = 5.4e10; // units: m^{-2} 
        // dislocation density needed to form a wall -- Thor. 2002 [26]
    double Mstress = 0.0;
    std::vector<double> c = {1.0, std::sqrt(2.0), std::sqrt(2.0), 1.0, std::sqrt(2.0), 1.0};
    std::transform(c.begin(), c.end(), stress.begin(),c.begin(), 
                   std::multiplies<double>());
    
    Mstress = tensorMagnitude(c)*std::sqrt(1.0/2.0);
    
    if (Mstress == 0.0 || Mrss/Mstress >= del || cDislDens < rhop)
        return 0;
    
    cDislDens -= rhop;
    cSize /= 2.0;
    
    // Select an orientation!
        // Should be toward max MRSS --away from vertical in uni. comp. or 
        // pure shear, towards vertical if simple shear.
    double theta = 0.0, phi = 0.0;
    angles(theta, phi);
    
    int stressIndex1 = (stress[4] > stress[0] ? 4 : 0);
    int stressIndex2 = (stress[1] > stress[2] ? 1 : 2);
    stressIndex2     = (stress[5] > stress[stressIndex2] ? 5 : stressIndex2);
    
    std::random_device seed;
    std::uniform_real_distribution<double> distribution(0.0,100.0);
    
    if (stress[stressIndex1] < stress[stressIndex2] && theta < M_PI/6.0) {
        // toward vertical
        theta -= M_PI/36.0;
        
    } else if (theta < M_PI/6.0) {
        // away from vertical
        theta += M_PI/36.0;
    } else {
        // randomly toward/away from vertical
        theta += (  distribution(seed) < 50.0 ? -M_PI/36.0 : M_PI/36.0);
    } 
    
    (void) axis(theta, phi);
    
    cTimeLastRecrystal = modelTime+timeStep;
    cSizeLastRecrystal = cSize;
    
    return 1;
}

void Crystal::rotate(const std::vector<double> &bigM, const std::vector<double> &bulkEdot, const std::vector<double> &stress, const double &timeStep) {
    
    std::vector<double> rdot, rdotTrans;
    rdot = tensorMixedInner(bigM, stress);
    rdotTrans = matrixTranspose(rdot, 3, 3);
    std::transform(rdot.begin(), rdot.end(), rdotTrans.begin(),rdot.begin(), 
                   std::minus<double>());
    std::transform(rdot.begin(),rdot.end(),rdot.begin(), 
                    [&](double x){return x/2.0;});
    
    std::vector<double> r = { 0, 1, 1,
                             -1, 0, 1,
                             -1,-1, 0};
    
    
    std::transform(r.begin(),r.end(), bulkEdot.begin(), r.begin(), 
                    std::multiplies<double>() );
    std::transform(r.begin(),r.end(), rdot.begin(), r.begin(), 
                    std::minus<double>() );
                    
    std::transform(r.begin(),r.end(),r.begin(), 
                    [&](double x){return x*timeStep;});
    // rotate crystals C' = (I+A+A^2/2+A^3/6+...)*C
    cAxis = vectorTimesExpm(cAxis,r, 3);
        // two terms seem to work fine. Three to be safe.
    
    
    //make unit vectors!
    double magCaxis;
    magCaxis = tensorMagnitude(cAxis);
    assert(magCaxis != 0.0);
    std::transform(cAxis.begin(),cAxis.end(),cAxis.begin(), 
                    [&](double x){return x/magCaxis;});
}

void Crystal::view() { 
    
    
    std::cout << "The C-Axis orientation is:" << std::endl;
    
    std::cout.precision(4);
    
    for (auto &i : cAxis)
        std::cout << std::fixed << i << " ";
    std::cout << std::endl;

    std::cout << "The crystal size is:" << std::endl;
    std::cout << std::scientific << cSize << " m" << std::endl;

    std::cout << "The crystal dislocation density is:" << std::endl;
    std::cout << std::scientific << cDislDens << " m^{-2}" << std::endl;
    
    std::cout.precision(0);
}

void Crystal::print() const {
     
    std::cout.precision(4);
    
    std::cout << std::fixed
              << std::setw(8) << cAxis[0]           << ", "
              << std::setw(8) << cAxis[1]           << ", "
              << std::setw(8) << cAxis[2]           << ", "
              << std::setw(11) << std::scientific
              << std::setw(11) << cSize              << ", "
              << std::setw(11) << cDislDens          << ", "
              << std::setw(11) << cTimeLastRecrystal << ", "
              << std::setw(11) << cSizeLastRecrystal << std::endl;
    
}
void Crystal::print(std::ofstream &file) const {
     
    file.precision(4);
    
    file << std::fixed
         << std::setw(8) << cAxis[0]           << ", "
         << std::setw(8) << cAxis[1]           << ", "
         << std::setw(8) << cAxis[2]           << ", "
         << std::setw(11) << std::scientific
         << std::setw(11) << cSize              << ", "
         << std::setw(11) << cDislDens          << ", "
         << std::setw(11) << cTimeLastRecrystal << ", "
         << std::setw(11) << cSizeLastRecrystal << std::endl;
    
}

void Crystal::angles(double &theta, double &phi) {
    
    double HXY = std::sqrt(cAxis[0]*cAxis[0] + cAxis[1]*cAxis[1]);
    theta = std::atan2(HXY, cAxis[2]);
    phi   = std::atan2(cAxis[1], cAxis[0]);
    
}

std::vector<double> Crystal::axis() {
    return cAxis;
}

std::vector<double> Crystal::axis(const double &theta, const double &phi) {
    
    cAxis = {sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
    
    double cAxisMag = tensorMagnitude(cAxis);
    
    assert(cAxisMag != 0.0);
    if (cAxisMag != 1.0) {
    std::transform(cAxis.begin(), cAxis.end(), cAxis.begin(), 
                   [&](double x){return x/cAxisMag;} );
    }

    return cAxis;
}

std::vector<double> Crystal::axis(const std::vector<double> &ax) {
    assert(ax.size() == 3);
    cAxis = ax;
    return cAxis;
}


std::vector<double> Crystal::all() const {
    return {cAxis[0], cAxis[1], cAxis[2], cSize, cDislDens, 
            cTimeLastRecrystal, cSizeLastRecrystal};
}

std::vector<double> Crystal::all(const std::vector<double> &al) {
    
    assert(al.size() == 7);
    
    (void) axis({al[0], al[1], al[2]});
    
    cSize     = al[3];
    cDislDens = al[4];
    
    cTimeLastRecrystal = al[5];
    cSizeLastRecrystal = al[6];
}

} // end of namespace FEvoR
