#!/bin/bash

clear

rm *.txt

# Problem Information (Fission, Parasitic Absorption, Leakage)
chains=10000
lfission=0.0
Pleakage=0.0

sed -i.bak "122s/.*/	integer, parameter :: loop = $chains/" NtrnGammaInit.f90
sed -i.bak "97s/.*/	real, parameter :: ntrnlfission = $lfission/" NtrnGammaInit.f90
sed -i.bak "103s/.*/	real, parameter :: ntrnPleakage = $Pleakage/" NtrnGammaInit.f90

#Spontaneous fission source or neutron present initial condition flag
#If ICNtrnFlag = 0, N(0) = 0, else N(0) = 1
ICNtrnFlag=1
sed -i.bak "215s/.*/	ICNtrnFlag = $ICNtrnFlag/" NtrnGammaInit.f90

fissflag=0
sed -i.bak "216s/.*/	fissflag = $fissflag/" NtrnGammaInit.f90

branchlens=1000

sed -i.bak "151s/.*/	integer, parameter :: branchlens = $branchlens/" NtrnGammaInit.f90
#sed -i.bak "7s/.*/	INTEGER, PARAMETER :: ntrnlens = $branchlens/" ntrngammadataread.f90

if [ $fissflag = 1 ]; then

	gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

else

	gfortran -fopenmp -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

fi	
#gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90
#./NtrnGammaInit.out

timeint=20

#let loop=10000

#for i in $(seq 1 $loop);
#do
(
	./NtrnGammaInit.out
	##echo $i
	#cat ntrnlifedata >> ntrnlife.txt
	#cat gammalifedata >> gammalife.txt
	#cat ntrnfissiondata >> ntrnfission.txt
	#cat gammafissiondata >> gammafission.txt
	##cat ntrntalarraytest.txt >> ntrntal.txt
	##cat gammatalarraytest.txt >> gammatal.txt

	#if [ $fissflag = 1 ]; then

		#cat SpontNuEmission.txt >> SpontNuEmiss.txt
		#cat SpontNuEmit.txt >> SpontNuNumber.txt
		#cat SpontMuEmit.txt >> SpontMuNumber.txt
		#rm SpontNuEmission.txt
		#rm SpontNuEmit.txt
		#rm SpontMuEmit.txt

	#fi
	
	#rm ntrnlifedata
	#rm gammalifedata
	#rm ntrnfissiondata
	#rm gammafissiondata 
	##rm ntrntalarraytest.txt
	##rm gammatalarraytest.txt
	) #&

#done

sed -i.bak "8s/.*/    INTEGER, PARAMETER :: N = $timeint/" Png.f90
sed -i.bak "9s/.*/    INTEGER, PARAMETER :: chains = $chains/" Png.f90

gfortran -o PNG.out Png.f90
./PNG.out

rm *.mod
rm *.bak

tput bel