/*
 */

#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <cmath>
#include "fevor_crystal.hh"
#include "test_crystal.hh"
#include "vector_tensor_opperations.hh"


void test_resolveM() {
    fevor_crystal c1(std::vector<double> {0,0,1},0.01,1e11);
    
    const double PI  =3.141592653589793238463;
    double theta, phi;
    theta = PI/4.0;
    phi = PI/4.0;
    c1.getNewAxis(theta, phi);
    
    std::vector<double> stress = {10000.0,     0.0,     0.0,
                                      0.0,     0.0,     0.0,
                                      0.0,     0.0,-10000.0};
    double temperature = -10.0;
    double Mrss = 1.0;
    double Medot = 1.0;
    
    std::vector<double> bigM;
    bigM = c1.resolveM(temperature, stress, Mrss, Medot);
    
    std::cout << "bigM size should be: 81" << "\n"
              << "bigM size is: " << bigM.size() << std::endl;
              
    std::vector<double> edot = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    edot = tensorMixedInner(bigM, stress);
    
    std::cout << "\n" << "edot is:" << std::endl;
    tensorDisplay(edot, 9, 1);
    
    //~ std::cout << "\n" << "bigM is:" << std::endl;
    //~ tensorDisplay(bigM, 9, 9);
    
    //~ std::cout << "\n" << "stress is:" << std::endl;
    //~ tensorDisplay(stress, 9, 1);
    
    
}

double test_grow(){
    fevor_crystal c1(std::vector<double> {1.0,0.0,0.0},0.01,10.0);
    //~ c1.seeCrystal();
    
    double temperature, model_time, K;
    temperature = -11.0; // units: degrees C
    model_time = 1000.0*365.0*24.0*60.0*60.0;   // units: m
    
    K = c1.grow(temperature, model_time);
    
    c1.seeCrystal();
    
    std::cout << "New Crystal Size should be: 1.0055e-2" << std::endl;
    
    return K;
}

void test_dislocate(const double &K) {
    fevor_crystal c1(std::vector<double> {1.0,0.0,0.0},1.0055e-2,1.0e10);
    //~ c1.seeCrystal();
    
    double timeStep = 1000.0*365.0*24.0*60.0*60.0;
    double Medot = sqrt(1.0/2.0); //TODO: Decide if Medot needs to be devided by 2!
    
    c1.dislocate(timeStep, Medot, K);
    
    c1.seeCrystal();
    
    std::cout << "Growth Constant was:" << K << "\n"
              << "New Dislocation Density should be: 4.9282e21" << std::endl;
}

void test_migRe() {
    fevor_crystal c1(std::vector<double> {1.0,0.0,0.0},0.01,1.0e11);
    
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
    fevor_crystal c1(std::vector<double> {0.0,0.0,1.0},0.01,1.0e11);
    
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
    
    fevor_crystal c1(std::vector<double> {1.0,0.0,0.0},1.0,10.0);
    c1.seeCrystal();
    
    double theta = 0.0, phi = 0.0;
    c1.getAxisAngles(theta, phi);
    
    std::cout.precision(4);
    std::cout << std::fixed << "Angles are: \n" << "Theta: " << theta << " Phi: " << phi << std::endl;
    
    const double PI  =3.141592653589793238463;
    std::random_device seed;
    std::uniform_real_distribution<double> dPhi(0.0,2.0*PI);
    std::uniform_real_distribution<double> dTheta(0.0,PI/2.0);
    
    theta = dTheta(seed);
    phi = dPhi(seed);
    
    std::cout.precision(4);
    std::cout << std::fixed << "New angles are: \n" << "Theta: " << theta << " Phi: " << phi << std::endl;
    
    c1.getNewAxis(theta, phi);
    
    c1.seeCrystal();
    
    
}


