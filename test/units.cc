#include <iostream>
#include <vector>
#include "../src/fevor_crystal.hh"
    // FIXME: bad practice! Fix includes.
    
int main()
{
    // Test the interface to fevor_crystal class
        std::vector<double> ca1 = {1,0,0};
        double cs = 1;
        long double cdd = 10;
        
        fevor_crystal c1(ca1,cs,cdd);
        c1.seeCrystal();
    
    return 0;
}

