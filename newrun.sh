#!/bin/bash

clear

rm *.txt
rm fort.*
rm SpontNu
#rm ntrnlifedata

# Problem Information (Fission, Parasitic Absorption, Leakage)
chains=50000
lfission=1.0
Pleakage=0.0

sed -i.bak "118s/.*/	integer, parameter :: loop = $chains/" NtrnGammaInit.f90
sed -i.bak "97s/.*/	real, parameter :: ntrnlfission = $lfission/" NtrnGammaInit.f90
sed -i.bak "103s/.*/	real, parameter :: ntrnPleakage = $Pleakage/" NtrnGammaInit.f90

#Spontaneous fission source or neutron present initial condition flag
#If ICNtrnFlag = 0, N(0) = 0, else N(0) = 1
ICNtrnFlag=1
sed -i.bak "129s/.*/	integer, parameter :: ICNtrnFlag = $ICNtrnFlag/" NtrnGammaInit.f90

fissflag=1
sed -i.bak "130s/.*/	integer, parameter :: fissflag = $fissflag/" NtrnGammaInit.f90

branchlens=1000

sed -i.bak "160s/.*/	integer, parameter :: branchlens = $branchlens/" NtrnGammaInit.f90

timeint=20

sed -i.bak "119s/.*/	integer, parameter :: N = $timeint/" NtrnGammaInit.f90

binaryfission=0

sed -i.bak "133s/.*/	integer, parameter :: binaryfissionflag = $binaryfission/" NtrnGammaInit.f90

if [ $fissflag = 1 ]; then

	gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

else

	gfortran -fopenmp -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

fi	

time ./NtrnGammaInit.out

rm *.mod
rm *.bak

tput bel

./case.sh

tput bel