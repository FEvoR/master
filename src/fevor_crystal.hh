/*
 */

#ifndef FEVOR_CRYSTAL
#define FEVOR_CRYSTAL

#include <vector>

class fevor_crystal {
    public:
    // constructors
        fevor_crystal(std::vector<double> ca, double cs, double cdd):
            cAxis(ca), cSize(cs), cDislDens(cdd) {
                cTimeLastRecrystal = 0.0;
                cSizeLastRecrystal = cSize;
            }
        
    // functions
        // get the strain-stress relation tensor (4th order)
        std::vector<double> resolveM(const double &temperature, const std::vector<double> &stress, double &Mrss, double &Medot);
        // grow the crystal
        double grow(const double &Tempature, const double &modelTime);
        // get new dislocation density
        void dislocate(const double &timeStep, const double &Medot, const double &K);
        // migration recrystallize if favorable to do so
        unsigned int migRe(const std::vector<double> &stress, const double &modelTime, const double &timeStep);
        // polygonize if favorable to do so
        unsigned int polygonize( const std::vector<double> &stress, const double &Mrss, const double &modelTime, const double &timeStep);
        // rotate the crystals
        void rotate(const std::vector<double> &bigM, const std::vector<double> &bulkEdot, const std::vector<double> &stress);
        
        void getAxisAngles(double &theta, double &phi);
        
        void getNewAxis(double &theta, double &phi);
        
        void seeCrystal();
        void printCrystal();
        
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
};



#endif
