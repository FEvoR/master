===============================================
 FEvoR: Fabric Evolution with Recrystallization
===============================================
Copyright (C) 2009-2017  Joseph H. Kennedy

This file is part of FEvoR

FEvoR is hosted at https://github.com/FEvoR

FEvoR is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

FEvoR is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
FEvoR.  If not, see <http://www.gnu.org/licenses/>.

Additional permission under GNU GPL version 3 section 7

If you modify FEvoR, or any covered work, to interface with other modules (such
as MATLAB code and MEX-files) available in a MATLAB(R) or comparable environment
containing parts covered under other licensing terms, the licensors of FEvoR
grant you additional permission to convey the resulting work.

=============
 About FEvoR
=============

FEvoR is a statistical fabric evolution model that includes 
recrystallization processes such as:
    Grain Growth
    Polygonization
    Migration Recrystallization

FEvoR was developed at the University of Alaska Fairbanks (UAF) by 
Joseph H. Kennedy. For a detailed description of this model, see:
 
    Kennedy, J. H., and E. C. Pettit (2015). The response of climate induced
    fabric variations to simple shear and migration recrystallization. Journal
    of Glaciology, 61 (227), 537--550. doi:10.3189/2015JoG14J156. 
    
    Kennedy, J. H., Pettit, E. C., & Di Prinzio, C. L. (2013). 
    The evolution of crystal fabric in ice sheets and its link to climate 
    history. Journal of Glaciology, 59(214), 357–373. 
    doi:10.3189/2013JoG12J159

FEvoR is based on the work by Dr. Thorsteinsson, see:

    Thorsteinsson, T. (2001). 
    An analytical approach to deformation of anisotropic ice-cryrstal 
    aggregates. Journal of Glaciology, 47(158), 507–516.
    
    Thorsteinsson, T. (2002). 
    Fabric development with nearest-neighbor interaction and dynamic 
    recrystallization. Journal of Geophysical Research Solid Earth, 
    107(B1,2014), 1–13. doi:10.1029/2001JB000244

====================
INSTALL INSTRUCTIONS
====================

FEvoR only requires the C++ standard library and C++11 compliant compiler.
FEvoR comes packaged with the [Faddeeva package (Dec.
2012)](http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package), and
[Catch v1.7.2](://github.com/philsquared/Catch) for unit testing.

To install:

Clone the FEvoR repository from github. 

```bash
git clone git@github.com:FEvoR/master.git fevor
```

Then you can make and install FEvoR like so:

```bash
cd fevor
mkdir build && cd build
cmake ..
make install
```
There is no uninstall routine within the make file. However, you can uninstall
the FEvoR library by running

```bash
xargs rm < install_manifest.txt
```

from within the build directory. 
