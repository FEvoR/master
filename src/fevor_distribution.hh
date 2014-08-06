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
        fevor_distribution(std::vector<unsigned int> lwh, double wk): dimensions(lwh), watsonK(wk)  {
            // TODO: check input of dimensions -- make sure vector is a length, width, and height
            
            numberCrystals = dimensions[0]*dimensions[1]*dimensions[2];
            for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
                
                crystals.push_back( fevor_crystal ({0.0,0.0,1.0}, 0.01, 1.0e10) );
                softness.push_back(1.0);
                magRSS.push_back(1.0);
                
                contribCrystal  = 1.0;
                contribNeighbor = 0.0;
                
                
            }
        }
        
        
    // functions
        // preform a time step
        std::vector<double> stepInTime(const double &temperature, const std::vector<double> &stress, const double &modelTime, const double &timeStep, double &nMigre, double &nPoly, std::vector<double> &bulkEdot);
        
        // calculate NNI softness parameter
        void getSoftness();
        
        // set softness ratio
        void setSoftness(double cc, double cn);
        
        // save distribution to disk
        void saveDistribution();
        
        
        
    private:
        std::vector<unsigned int> dimensions;
        
        double watsonK;

        unsigned int  numberCrystals;
        
        std::vector<fevor_crystal> crystals;
        
        std::vector<double> softness;
        std::vector<double> magRSS;
        
        double contribCrystal;
        double contribNeighbor;
    
    
};

#endif 
