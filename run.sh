#!/bin/bash

# rm test.txt

clear

rm ntrnlife.txt
rm gammalife.txt
rm ProbNStats.txt
rm NtrnGammaInit.out
rm ntrnfission.txt
rm gammafission.txt

gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

# touch test.txt

loop=5000
chain=10000
#idx = expr $chain/$loop
#echo $idx
timeint=50

for j in $(seq 1 10);
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
sed -i.bak "8s/.*/  INTEGER, PARAMETER :: N = $timeint/" Png.f90
sed -i.bak "9s/.*/  INTEGER, PARAMETER :: chains = $loop/" Png.f90

gfortran -o ntrngammadatread ntrngammadataread.f90

./ntrngammadatread

#sed -i.bak "59s/.*/  open( unit = 3, file = \"ProbN_$j.txt\" )" Png.f90

gfortran -o PNG.out Png.f90

./PNG.out

cat ProbN.txt >> ProbNStats.txt
cat ProbG.txt >> ProbGStats.txt

#rm ProbN.txt

done

rm ntrngammadatread
rm mat_params.mod
rm mcnp_params.mod
rm mcnp_random.mod
rm ntrngammadataread.f90.bak