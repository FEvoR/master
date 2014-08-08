/*
 */

#ifndef FEVOR_DISTRIBUTION
#define FEVOR_DISTRIBUTION

#include <vector>
#include "fevor_crystal.hh"

class fevor_distribution {
    public:
    // constructors
        // TODO: check input of dimensions -- make sure vector is a length, width, and height
        fevor_distribution(std::vector<unsigned int> lwh): dimensions(lwh) {
            numberCrystals = dimensions[0]*dimensions[1]*dimensions[2];
            for (unsigned int ii = 0; ii!= numberCrystals; ++ii) {
                crystals.push_back( fevor_crystal ({0.0,0.0,1.0}, 0.01, 1.0e10) );
                softness.push_back(1.0);
                magRSS.push_back(1.0);
                
                contribCrystal  = 1.0;
                contribNeighbor = 0.0;
            }
        }
        // construct a distribution from a saved distribution of crystals
        fevor_distribution(std::vector<unsigned int> lwh, std::string fname): fevor_distribution(lwh)  {
            loadDistribution(fname);
        }
        // construct a distribution using the Watson distribution for axis angles
        fevor_distribution(std::vector<unsigned int> lwh, double wk): fevor_distribution(lwh)  {
            generateWatsonAxes(wk);
        }
        
        
    // functions
        // preform a time step
        std::vector<double> stepInTime(const double &temperature, const std::vector<double> &stress, double &modelTime, const double &timeStep, double &nMigre, double &nPoly, std::vector<double> &bulkEdot);
        
        // calculate NNI softness parameter
        void getSoftness(std::vector<std::vector<double>> &crystalM, std::vector<double> &bulkM, std::vector<double> &bulkEdot, const std::vector<double> &stress);
        
        // set softness ratio
        void setSoftness(double cc, double cn);
        
        // save distribution to disk
        void saveDistribution();
        void saveDistribution(std::string fname);
        
        // save distribution to disk
        void loadDistribution( std::string fname );
        
        //TODO: generate watson cAxis for distrobution
        void generateWatsonAxes(const double &wk);
        
    private:
        std::vector<unsigned int> dimensions;

        unsigned int  numberCrystals;
        
        std::vector<fevor_crystal> crystals;
        
        std::vector<double> softness;
        std::vector<double> magRSS;
        
        double contribCrystal;
        double contribNeighbor;
    
    
};

#endif 
