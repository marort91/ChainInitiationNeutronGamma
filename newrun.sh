#!/bin/bash

clear

rm *.txt
#rm ntrnlifedata

# Problem Information (Fission, Parasitic Absorption, Leakage)
chains=10000
lfission=0.01
Pleakage=0.0

sed -i.bak "118s/.*/	integer, parameter :: loop = $chains/" NtrnGammaInit.f90
sed -i.bak "97s/.*/	real, parameter :: ntrnlfission = $lfission/" NtrnGammaInit.f90
sed -i.bak "103s/.*/	real, parameter :: ntrnPleakage = $Pleakage/" NtrnGammaInit.f90

#Spontaneous fission source or neutron present initial condition flag
#If ICNtrnFlag = 0, N(0) = 0, else N(0) = 1
ICNtrnFlag=0
sed -i.bak "127s/.*/	integer, parameter :: ICNtrnFlag = $ICNtrnFlag/" NtrnGammaInit.f90

fissflag=1
sed -i.bak "128s/.*/	integer, parameter :: fissflag = $fissflag/" NtrnGammaInit.f90

branchlens=1000

sed -i.bak "152s/.*/	integer, parameter :: branchlens = $branchlens/" NtrnGammaInit.f90
#sed -i.bak "7s/.*/	INTEGER, PARAMETER :: ntrnlens = $branchlens/" ntrngammadataread.f90

timeint=25

sed -i.bak "119s/.*/	integer, parameter :: N = $timeint/" NtrnGammaInit.f90

if [ $fissflag = 1 ]; then

	gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

else

	gfortran -fopenmp -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

fi	

time ./NtrnGammaInit.out

#sed -i.bak "8s/.*/    INTEGER, PARAMETER :: N = $timeint/" Png.f90
#sed -i.bak "9s/.*/    INTEGER, PARAMETER :: chains = $chains/" Png.f90
#sed -i.bak "12s/.*/	REAL, PARAMETER :: chain = $chains/" PngTest.f90

sed -i.bak "8s/.*/    INTEGER, PARAMETER :: N = $timeint/" PngTest.f90
sed -i.bak "9s/.*/    INTEGER, PARAMETER :: chains = $chains/" PngTest.f90
sed -i.bak "12s/.*/	 REAL, PARAMETER :: chain = $chains/" PngTest.f90

#gfortran -o PNG.out Png.f90
gfortran -o PNG.out PngTest.f90
./PNG.out

rm *.mod
rm *.bak

tput bel

./case.sh

tput bel