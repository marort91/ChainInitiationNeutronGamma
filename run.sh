#!/bin/bash

# rm test.txt
rm ntrnlife.txt
rm FissionNuTwo

gfortran -o FissionNuTwo mcnp_random.f90 meeting212.f90

# touch test.txt

chain=5e4

for i in $(seq 1 $chain);
do
	./FissionNuTwo
	echo $i
	cat lifedata >> ntrnlife.txt
	rm lifedata
done

sed -i.bak "8s/.*/	INTEGER, PARAMETER :: chains = $chain/" dataread.f90

gfortran -o datread dataread.f90

./datread

rm datread
rm mat_params.mod
rm mcnp_params.mod
rm mcnp_random.mod
rm dataread.f90.bak