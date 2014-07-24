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
        fevor_crystal(std::vector<double> ca, double cs, double cdd):
            cAxis(ca), cSize(cs), cDislDens(cdd) {
                cTimeLastRecrystal = 0;
                cSizeLastRecrystal = cSize;
            }
        
    // functions
        void resolveScalar();
        
        // grow the crystal
        double grow(const double &Tempature, const double &stepTime);
        // get new dislocation density
        void dislocate(const double &timeStep, const double &Medot, const double &K);
        // migration recrystallize if favorable to do so
        void migRe(const double &Estress, const double &modelTime, const double &timeStep);
        // polygonize if favorable to do so
        void polygonize();
        
        void resolveM();
        
        void seeCrystal();
        
        
    private:
        // holds the crystals c-axis orientation vector in cartesian coordinates
        std::vector<double> cAxis; // unit vector
        // holds the size of the crystal
        double cSize; // units: m
        // holds the dislocation density of the crystal
        double cDislDens; // units: #/m^2
        // holds the time of last recrystallization
        double cTimeLastRecrystal; // units: #/m^2
        // holds the cSize at last recrystallization
        double cSizeLastRecrystal; // units: m
        
        // Static members that are needed to compute member functions
        
        
        
        // Default values for static members
};



#endif
