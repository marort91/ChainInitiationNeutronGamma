#!/bin/bash

# rm test.txt

clear

rm ntrnlife.txt
rm gammalife.txt
rm NtrnGammaInit.out

gfortran -o NtrnGammaInit.out mcnp_random.f90 NtrnGammaInit.f90

# touch test.txt

chain=1e3

for i in $(seq 1 $chain);
do
	./NtrnGammaInit.out
	echo $i
	cat ntrnlifedata >> ntrnlife.txt
	cat gammalifedata >> gammalife.txt
	rm ntrnlifedata
	rm gammalifedata
done

sed -i.bak "9s/.*/	INTEGER, PARAMETER :: chains = $chain/" ntrngammadataread.f90

gfortran -o datread ntrngammadataread.f90

./datread

rm datread
rm mat_params.mod
rm mcnp_params.mod
rm mcnp_random.mod
#rm ntrngammadataread.f90.bak