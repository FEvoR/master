/**************************************************
 * FEvoR: Fabric Evolution with Recrystallization *
 **************************************************
 * Copyright (C) 2009-2014  Joseph H Kennedy
 *
 * This file is part of FEvoR.
 *
 * FEvoR is free software: you can redistribute it and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version.
 *
 * FEvoR is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with FEvoR.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Additional permission under GNU GPL version 3 section 7
 *
 * If you modify FEvoR, or any covered work, to interface with
 * other modules (such as MATLAB code and MEX-files) available in a
 * MATLAB(R) or comparable environment containing parts covered
 * under other licensing terms, the licensors of FEvoR grant
 * you additional permission to convey the resulting work.
 */

/*
 */

#ifndef FEVOR_DISTRIBUTION
#define FEVOR_DISTRIBUTION

#include <vector>
#include "fevor_crystal.hh"

namespace FEvoR {

class Distribution {
    public:
    // constructors
    // TODO: check input of dimensions -- make sure vector is a length, width, and height
    Distribution(std::vector<unsigned int> lwh);
    // construct a distribution from a saved distribution of crystals
    Distribution(std::vector<unsigned int> lwh, std::string fname);
    // construct a distribution using the Watson distribution for axis angles
    Distribution(std::vector<unsigned int> lwh, double wk);
    // construct a distribution from a big vector of all distribution data
    Distribution(std::vector<unsigned int> lwh, std::vector<double> &data);
        
    // functions
        // preform a time step
        std::vector<double> stepInTime(const double &temperature, const std::vector<double> &stress, const double &modelTime, const double &timeStep, unsigned int &nMigre, unsigned int &nPoly, std::vector<double> &bulkEdot);
        // preform a time step without tracking recrystallization
        std::vector<double> stepInTime(const double &temperature, const std::vector<double> &stress, const double &modelTime, const double &timeStep, std::vector<double> &bulkEdot);
        
        // calculate NNI softness parameter
        void getSoftness(std::vector<std::vector<double> > &crystalM, std::vector<double> &crystalMagEdot, std::vector<double> &bulkM, std::vector<double> &bulkEdot, const std::vector<double> &stress);
        
        // set softness ratio
        void setSoftnessRatio(double cc, double cn);
        
        // save distribution to disk
        void saveDistribution() const;
        void saveDistribution(std::string fname) const;
        void saveDistribution(std::vector<double> &data) const;
        
        // Load distribution from disk
        void loadDistribution( std::string fname );
        void loadDistribution( const std::vector<double> &data );
        //generate watson cAxis for distribution
        void generateWatsonAxes(const double &wk);
        
        unsigned int getNumberCrystals() const;
        
  static const unsigned int numberParameters;
    private:
        std::vector<unsigned int> dimensions;

        unsigned int  numberCrystals;
        
        std::vector<Crystal> crystals;
        
        std::vector<double> softness;
        std::vector<double> magRSS;
        
        double contribCrystal;
        double contribNeighbor;
    
    
};

} // end of namespace FEvoR

#endif 
