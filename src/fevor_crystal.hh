/**************************************************
 * FEvoR: Fabric Evolution with Recrystallization *
 **************************************************
 * Copyright (C) 2009-2017  Joseph H. Kennedy
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

#ifndef FEVOR_CRYSTAL
#define FEVOR_CRYSTAL

#include <vector>
#include <fstream>

namespace FEvoR {

class Crystal {
    public:
    // constructors
    Crystal(std::vector<double> ca, double cs, double cdd);
    Crystal(std::vector<double> ca, double cs, double cdd, double cto, double cso);

    // functions
        // get the strain-stress relation tensor (4th order)
        std::vector<double> resolveM(const double &temperature, const std::vector<double> &stress, double &Mrss, double &Medot);
        // grow the crystal
        double grow(const double &temperature, const double &modelTime);
        // get new dislocation density
        void dislocate(const double &timeStep, const double &Medot, const double &K);
        // migration recrystallize if favorable to do so
        unsigned int migRe(const std::vector<double> &stress, const double &modelTime, const double &timeStep);
        // polygonize if favorable to do so
        unsigned int polygonize( const std::vector<double> &stress, const double &Mrss, const double &modelTime, const double &timeStep);
        // rotate the crystals
        void rotate(const std::vector<double> &bigM, const std::vector<double> &bulkEdot, const std::vector<double> &stress, const double &timeStep);

        void seeCrystal();

        void printCrystal() const;
        void printCrystal(std::ofstream &file) const;

        void getAxis(double &ca0, double &ca1, double &ca2);
        void getAxisAngles(double &theta, double &phi);

        void setNewAxis(const double &theta, const double &phi);
        void setNewAxis(const std::vector<double> &ax);

        void setAll(const double &ca1, const double &ca2, const double &ca3,
                    const double &csz, const double &cdd,
                    const double &ctlr, const double &cslr);

        void getAll(double &ca0, double &ca1, double &ca2,
                    double &csz, double &cdd, double &ctlr, double &cslr) const;
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

} // end of namespace FEvoR

#endif
