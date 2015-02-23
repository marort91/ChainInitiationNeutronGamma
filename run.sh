#!/bin/bash

# rm test.txt

clear

rm ntrnlife.txt
rm gammalife.txt
rm ProbGStats.txt
rm ProbNStats.txt
rm NtrnGammaInit.out
rm ntrnfission.txt
rm gammafission.txt

# Problem Information (Fission, Parasitic Absorption)
lfission=0.25

sed -i.bak "70s/.*/	real, parameter :: lfission = $lfission/" NtrnGammaInit.f90

gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

# touch test.txt

let loop=1000
let chain=100000
let idx=chain/loop

timeint=50

for j in $(seq 1 $idx);
do

	rm ntrnlife.txt
	rm gammalife.txt

for i in $(seq 1 $loop);
do
	./NtrnGammaInit.out
	echo $i
	cat ntrnlifedata >> ntrnlife.txt
	cat gammalifedata >> gammalife.txt
	cat ntrnfissiondata >> ntrnfission.txt
	cat gammafissiondata >> gammafission.txt
	rm ntrnlifedata
	rm gammalifedata
	rm ntrnfissiondata
	rm gammafissiondata
done

echo $j

sed -i.bak "6s/.*/	INTEGER, PARAMETER :: N = $timeint/" ntrngammadataread.f90
sed -i.bak "9s/.*/	INTEGER, PARAMETER :: chains = $loop/" ntrngammadataread.f90
sed -i.bak "8s/.*/    INTEGER, PARAMETER :: N = $timeint/" Png.f90
sed -i.bak "9s/.*/    INTEGER, PARAMETER :: chains = $loop/" Png.f90

gfortran -o ntrngammadatread ntrngammadataread.f90
./ntrngammadatread

gfortran -o PNG.out Png.f90
./PNG.out

cat ProbN.txt >> ProbNStats.txt
cat ProbG.txt >> ProbGStats.txt

done

rm ntrngammadatread
#rm mat_params.mod
#rm mcnp_params.mod
#rm mcnp_random.mod
#rm ntrngammadataread.f90.bak

rm *.mod
rm *.bak

sed -i.bak "5s/.*/lf = $lfission;/" PnStats.m
sed -i.bak "26s/.*/N = $timeint;/" PnStats.m
sed -i.bak "25s/.*/chains = $chain;/" PnStats.m

sed -i.bak "8s/.*/chains = $chain;/" PnPlotter.m
sed -i.bak "10s/.*/lf = $lfission;/" PnPlotter.m
sed -i.bak "13s/.*/N = $timeint;/" PnPlotter.m