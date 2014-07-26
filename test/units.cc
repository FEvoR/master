#include <iostream>
#include <vector>
#include "../src/fevor_crystal.hh"
    // FIXME: bad practice! Fix includes.
    
int main()
{
    // Test the interface to fevor_crystal class
        fevor_crystal c1(std::vector<double> {1,0,0},1,10);
        c1.seeCrystal();
        double T = -9;
        double t = 1e25;
        double k = c1.grow(T, t);
        std::cout << "\n Growth constant: " << k << "\n " << std::endl;
        
        c1.dislocate(1, 1, k);
        
        c1.seeCrystal();
        
        unsigned int numberMigRe = 0;
        numberMigRe += c1.migRe(1, t, t);
        std::cout << "\n Times MigRed: " << numberMigRe << "\n " << std::endl;
        
        unsigned int numberPoly = 0;
        numberPoly += c1.polygonize(1, 1, t, t);
        std::cout << "\n Times Polygonized: " << numberPoly << "\n " << std::endl;
        
        c1.seeCrystal();
        
    return 0;
}

