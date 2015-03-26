#!/bin/bash

clear

rm *.txt

# Problem Information (Fission, Parasitic Absorption, Leakage)
chains=1000
lfission=0.5
Pleakage=0.0

sed -i.bak "121s/.*/	integer, parameter :: loop = $chains/" NtrnGammaInit.f90
sed -i.bak "97s/.*/	real, parameter :: ntrnlfission = $lfission/" NtrnGammaInit.f90
sed -i.bak "103s/.*/	real, parameter :: ntrnPleakage = $Pleakage/" NtrnGammaInit.f90

#Spontaneous fission source or neutron present initial condition flag
#If ICNtrnFlag = 0, N(0) = 0, else N(0) = 1
ICNtrnFlag=1
sed -i.bak "173s/.*/	ICNtrnFlag = $ICNtrnFlag/" NtrnGammaInit.f90

fissflag=0
sed -i.bak "174s/.*/	fissflag = $fissflag/" NtrnGammaInit.f90

branchlens=1000

sed -i.bak "143s/.*/	integer, parameter :: branchlens = $branchlens/" NtrnGammaInit.f90
#sed -i.bak "7s/.*/	INTEGER, PARAMETER :: ntrnlens = $branchlens/" ntrngammadataread.f90

gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90
./NtrnGammaInit.out

timeint=20

sed -i.bak "8s/.*/    INTEGER, PARAMETER :: N = $timeint/" Png.f90
sed -i.bak "9s/.*/    INTEGER, PARAMETER :: chains = $chains/" Png.f90

gfortran -o PNG.out Png.f90
./PNG.out

rm *.mod
rm *.bak

tput bel