/*
 */

#include <vector>
#include <iostream>
#include "fevor_crystal.hh"
 
// Define function members
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
 
