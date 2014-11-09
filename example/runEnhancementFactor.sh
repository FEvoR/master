#!/bin/bash

# enhancementFactor watsonK temperature Sxx Sxy Sxz Syy Syz Szz dt cc cn 
#                      -      Kelvin   Pa  Pa  Pa  Pa  Pa  Pa  yr -  -

# path to executable
path=${1%/}
# Watson concentration parameter
k=${2}
# temperature
T=266.98
####### stress ######
#                   #
# | Sxx, Sxy, Sxz | #
# | Sxy, Syy, Syz | #
# | Sxz, Syz, Szz | #
#                   #
#####################
sxx=${4}
sxy=0.0
sxz=${3}
syy=${4}
syz=0.0
szz=${4}
# time step
dt=10
# nearest neighbor interaction
cc=1.0 # contribution of crystal
cn=0.0 # contribution of neighbors


${path}/enhancementFactor ${k} ${T} ${sxx} ${sxy} ${sxz} ${syy} ${syz} ${szz} ${dt} ${cc} ${cn}
