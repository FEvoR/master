/*
 */

#include <vector>
#include <iostream>
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
    
    c1.getNewAxis(theta, phi);
    
    c1.seeCrystal();
    
}
