/*
 */

#ifndef FEVOR_DISTRIBUTION
#define FEVOR_DISTRIBUTION

#include <vector>
#include "fevor_crystal.hh"

class fevor_distribution {
    public:
    // constructors
        // construct a distribution from a saved distribution
        
        // construct a distribution using the Watson distribution for axis angles
        fevor_distribution(unsigned int nc, double wk): numberCrystals(nc), watsonK(wk) {
            for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
                
                crystals.push_back( fevor_crystal ({0.0,0.0,1.0}, 0.01, 1.0e10) );
                
                crystals[ii].seeCrystal();
                
            }
        }
        
        
    // functions
        
        // save distribution to disk
        
        // calculate distribution eigenvalues
        
        // calculate distribution bulk M
        
        // calculate NNI softness parameter
        
        
        
    private:
    unsigned int  numberCrystals;

    double watsonK;

    std::vector<fevor_crystal> crystals;
    
    std::vector<double> softness;
    
    std::vector<fevor_crystal> bulkM;
    
    std::vector<fevor_crystal> bulkEdot;
    

    

    
};

#endif 
