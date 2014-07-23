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
        std::cout << "Growth constant: " << k << std::endl;
        
        c1.dislocate(t, 1, k);
        
        c1.seeCrystal();
        
    return 0;
}

