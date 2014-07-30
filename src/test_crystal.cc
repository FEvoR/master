/*
 */

#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include "fevor_crystal.hh"
#include "test_crystal.hh"

void test_resolveM() {
    
    
    
    
}

void test_grow(){
    
    
    
}

void test_dislocate() {
    
    
    
}

void test_migRe() {
    
    
    
}

void test_polygonize() {
    
    
    
}

void test_angles() {
    
    fevor_crystal c1(std::vector<double> {1,0,0},1,10);
    c1.seeCrystal();
    
    double theta = 0, phi = 0;
    c1.getAxisAngles(theta, phi);
    
    std::cout.precision(4);
    std::cout << std::fixed << "Angles are: \n" << "Theta: " << theta << " Phi: " << phi << std::endl;
    
    const double PI  =3.141592653589793238463;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<double> dPhi(0.0,2*PI);
    std::uniform_real_distribution<double> dTheta(0.0,PI/2);
    
    theta = dTheta(generator);
    phi = dPhi(generator);
    
    std::cout.precision(4);
    std::cout << std::fixed << "New angles are: \n" << "Theta: " << theta << " Phi: " << phi << std::endl;
    
    c1.getNewAxis(theta, phi);
    
    c1.seeCrystal();
    
    
}
