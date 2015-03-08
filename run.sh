#!/bin/bash

clear

rm *.txt

# Problem Information (Fission, Parasitic Absorption)
lfission=0.0

sed -i.bak "77s/.*/	real, parameter :: lfission = $lfission/" NtrnGammaInit.f90

gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

# touch test.txt

let loop=1000
let chain=10000
let idx=chain/loop

timeint=50

for j in $(seq 1 $idx);
do

	#rm ntrnlife.txt
	#rm gammalife.txt

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

#Output moment data here using some sort of reading of data within Png.f90

cat PnMmnt0.txt >> ProbN0Stats.txt
cat PnMmnt1.txt >> ProbN1Stats.txt
cat PnMmnt2.txt >> ProbN2Stats.txt
cat PnMmnt3.txt >> ProbN3Stats.txt
cat PnMmnt4.txt >> ProbN4Stats.txt
#cat ProbG.txt >> ProbGStats.txt
cat PgMmnt0.txt >> ProbG0Stats.txt
cat PgMmnt1.txt >> ProbG1Stats.txt
cat PgMmnt2.txt >> ProbG2Stats.txt
cat PgMmnt3.txt >> ProbG3Stats.txt
cat PgMmnt4.txt >> ProbG4Stats.txt

done

rm ntrngammadatread
#rm mat_params.mod
#rm mcnp_params.mod
#rm mcnp_random.mod
#rm ntrngammadataread.f90.bak

rm *.mod
rm *.bak
#rm fort.*

sed -i.bak "5s/.*/lf = $lfission;/" PnStats.m
sed -i.bak "26s/.*/N = $timeint;/" PnStats.m
sed -i.bak "25s/.*/chains = $chain;/" PnStats.m

sed -i.bak "8s/.*/chains = $chain;/" PnPlotter.m
sed -i.bak "10s/.*/lf = $lfission;/" PnPlotter.m
sed -i.bak "13s/.*/N = $timeint;/" PnPlotter.m