/*
 */

#ifndef FEVOR_CRYSTAL
#define FEVOR_CRYSTAL

#include <vector>
#include <iostream>

class fevor_crystal {
    public:
        // constructors
        fevor_crystal() = default;
        fevor_crystal(std::vector<double> ca, double cs, long double cdd):
            cAxis(ca), cSize(cs), cDislDens(cdd) {/* TODO: test validity of inputs */}
        
        // functions
        void seeCrystal() { 
            
            std::cout << "The C-Axis orientation is:" << std::endl;
            for (auto i : cAxis)
                std::cout << i << " ";
            std::cout << std::endl;
            
            std::cout << "The crystal size is:" << std::endl;
            std::cout << cSize << " m" << std::endl;
            
            std::cout << "The crystal dislocation density is:" << std::endl;
            std::cout << cDislDens << " m^{-2}" << std::endl;
        }
        
        
    private:
        // holds the crystals c-axis orientation vector in cartesian coordinates
        std::vector<double> cAxis; // unit vector
        
        // holds the size of the crystal
        double cSize; // unit: m
        
        // holds the dislocation density of the crystal
        long double cDislDens; // unit: #/m^2
};



#endif
