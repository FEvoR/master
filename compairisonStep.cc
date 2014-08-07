#include <iostream>
#include <iomanip>
#include <vector>
#include "fevor_distribution.hh"


// Perform a compairison time step between FEvoR and Thor

// TODO: figure out how to time the step
int main()
{
    double temperature = -10.0;
    std::vector<double> stress = { 10000,     0, 10000,
                                       0,     0,     0,
                                   10000,     0,-10000};
    double modelTime = 0.0;
    double timeStep = 1000.0*365.0*24.0*60.0*60.0;
    double nMigre, nPoly;
    nMigre = nPoly = 0;
    std::vector<double> bulkEdot(9, 0.0);
    
    std::vector<unsigned int> packingDimensions(3, 20);
    double watsonK = 0.0;
    
    fevor_distribution d1(packingDimensions, watsonK);
    
    std::cout << "\n" << "# Initial Distribution: \n" << std::endl;
    d1.saveDistribution();
    
    d1.stepInTime(temperature, stress, modelTime, timeStep, nMigre, nPoly, bulkEdot);
    
    std::cout << "\n" << "# Stepped Distribution: \n" << std::endl;
    d1.saveDistribution();
    
    return 0;
}
