!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																											         !!
!! 	NtrnGammaInit.f90																						         !!
!! 	Description: Neutron and photon stochastic population tracking										             !!
!! 	Assumptions: Zero-dimension, one speed neutron and gamma ray creation									         !!
!! 	Created by: Mario I. Ortega																				         !!
!! 	Created: January 26th, 2015																					     !!
!!																			                                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																													 !!
!! Modules																											 !!
!! 																													 !!
!! mcnp_random		--> MCNP Random Number Generator (External)														 !!
!! mcnp_params		--> MCNP Fortran Options																		 !!
!! mat_params		--> Problem material neutron cross sections, multiplicity, reaction rate frequncy				 !!
!! timestmp(REMOVED)--> Initializes time stamp data for particular branch											 !!
!!																													 !!
!! Subroutines																										 !!
!!																													 !!
!! Functions																										 !!
!!																													 !!
!! TimeInteract     -->	Interaction time for neutron after birth													 !!
!! NtrnMult			--> Neutron multiplicity PDF and sampling														 !!
!! GammaMult		--> Gamma ray multiplicity PDF and sampling														 !!
!! SpontaneousSrc   --> Spontaneous neutron source emission															 !!					
!!																													 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!																													 !!
!! Revision History																									 !!
!!																													 !!
!! 01/29/2015																										 !!
!! - Fission function made not implemented. Instead function code resides in main program. 							 !!
!! - timestmp module removed.																						 !!
!!																													 !!
!! 02/04/2015																										 !!
!! - Problem is in the nu term in this code. Code works when there is no multiplication. Indexing???				 !!
!! - Implemented workaroud preventing segmentation faults. Currently code allows for any nu and prevents code from   !!
!!   accessing non-allocated array indices. The code then interacts all created neutrons and records time.           !!
!!																													 !!
!! 02/09/2015																										 !!
!! - Implemented neutron multiplicity function. Currently uses hard coded probability density function. Improvements !!
!!   could include reading in a PDF from an external text file directly into the function for different problems.    !!
!!																												     !!
!! 02/10/2015																										 !!
!! - Implemented gamma photon multiplicity function. Currently uses hard coded PDF. Improvements could include 		 !!
!!   reading PDF from external text file directly into the function.												 !!
!! - Implemened spontaneous neutron source emission function. Need to verify it works as intended					 !!
!!																													 !!
!! 02/13/2015																										 !!
!! - Implemented mubar and nubar (average number of neutrons and photons emitted per fission) in the GammaMult       !!
!!   and NtrnMult functions. Code verified for single source neutron present at time zero.							 !!
!!																													 !!
!! 02/16/2015																										 !!
!! - Fixed neutron and gamma multiplicity functions. Verified no fission possible case. No more segmentation 		 !!
!!   faults after fixing indexing. Nubar, mubar, and cumulative distribution functions implemented correctly!		 !!
!!																													 !!
!! 02/23/2015																										 !!
!! - Implemented spontaneuous fission source (currently emitting one neutron per second). Seems to work with 		 !!
!!   fission though it remains to be seen what happens to population after increasing branch length.				 !!
!!																													 !!
!! 02/26/2015																										 !!
!! -Implemented spontaneuous and initial neutron flag for different case analysis. Moment calculation through time   !!
!!  to be implemented where moments are defined as sum (n**i)*Pn where i is the particular moment to be calculated   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mcnp_params

	integer, parameter :: R8 = selected_real_kind(15,307)
    integer, parameter :: I8 = selected_int_kind(18)

    integer(I8) :: seed

END MODULE mcnp_params
	
MODULE mat_params

	real, parameter :: lfission = 0.0
    real, parameter :: lcapture = 1.0 - lfission
    real, parameter :: ltot = lfission + lcapture

    real, parameter :: Pfiss = lfission / ltot
    real, parameter :: Pcap = lcapture / ltot
    integer :: nu, mu, spontNu

    !Binary Fission Assumption
    !integer, parameter :: nu = 2

END MODULE mat_params

PROGRAM neutrongammachain

	use mcnp_random
	use mcnp_params
	use mat_params

	IMPLICIT NONE

	!INTERFACE
		
	!	SUBROUTINE SpotaneousSrc(rnm,t0,time,SpontNu)

	!	real, intent(IN) :: rnm, t0
	!	real, intent(OUT) :: time
	!	integer, intent(OUT) :: SpontNu

	!	END SUBROUTINE SpotaneousSrc

	!END INTERFACE	

	real(R8) :: dt, rnmN, rnmG, rnmSp, rnmSpNu
	real(R8) :: t0, tf, TimeInteract, SpotaneousSrcTime, tsp
	integer(I8) :: i, j, NtrnMult, GammaMult, SpotaneousNu
	integer :: nidx, gidx, fissflag

	integer, parameter :: branchlens = 1000
	real(R8), dimension(1,2,branchlens) :: ntrntime = 0
	real(R8), dimension(1,branchlens*2) :: gammatime = 0

	!Initialize MCNP Random Number Generator
	call system_clock(seed)
	call RN_init_problem( 2, seed, 0_I8, 0_I8, 1 )

	open( unit = 1, file = "ntrnlifedata" )
	open( unit = 2, file = "gammalifedata" )
	open( unit = 3, file = "ntrnfissiondata")
	open( unit = 4, file = "gammafissiondata")

	!Initialize neutron and gamma indices
	nidx = 0
	gidx = 1

	!Birth time of first neutron ( Random source to be implemented )
	!Chain reaction start time
	t0 = 0

	!rnmSp = rang()
	!ntrntime(1,1,1) = SpotaneousSrc(rnmSp)

	!Spontaneous fission ( fissflag = 1 ) or neutron present at start time flag ( fissflag = 0 )

	fissflag = 0

	do i=1,branchlens

		if ( fissflag .eq. 1 ) then

			rnmSp = rang()
			tsp = SpotaneousSrcTime(t0,rnmSp)

			rnmSpNu = rang()
			spontNu = SpotaneousNu(rnmSpNu)

			ntrntime(1,1,nidx+1:nidx+spontNu) = tsp
		
			nidx = nidx + spontNu
			t0 = tsp

		endif

		if ( ( ntrntime(1,1,i).eq.0_R8) .and. ( i.gt.1 ) ) exit
		
		ntrntime(1,2,i) = ntrntime(1,1,i) + TimeInteract(ntrntime(1,1,i))

		if (rang() .le. Pfiss) then

			rnmN = rang()
			rnmG = rang()

			nu = NtrnMult(rnmN)
			mu = GammaMult(rnmG)

			write( 3, * ), nu
			write( 4, * ), mu
			!print *, rnmN
			!print *
			!print *, "Time of Fission: ", ntrntime(1,2,i)
			!print *
			!print *, "Fission neutrons created: ", nu
			!print *, "Gammas generated:         ", mu
			!print *, nubar

			ntrntime(1,1,nidx+1:nidx+nu) = ntrntime(1,2,i)
			gammatime(1,gidx:gidx+mu-1) = ntrntime(1,2,i)
			
			nidx = nidx + nu
			gidx = gidx + mu

			if ( nu .eq. 0 ) then

				continue

			else

			!print *, "Gamma index: 		   ", gidx
			!print *

			endif
			!print *, idx

		endif

		if ( nidx + 1 .gt. branchlens) exit

		!rnmSp = rang()
		!nidx = nidx + 1


	enddo

	do i=1,branchlens

		if ( (ntrntime(1,1,i) .ne. 0) .and. (ntrntime(1,2,i) .eq. 0) ) then

			!time(1,2,i) = time(1,1,i) - (1/ltot)*log( rang() )
			ntrntime(1,2,i) = ntrntime(1,1,i) + TimeInteract(ntrntime(1,1,i))

		endif

	enddo

	do i=1,branchlens
		
		write( 1, * ), ntrntime(1,1,i), ntrntime(1,2,i)

	enddo

	do i=1,2*branchlens

		write( 2, * ), gammatime(1,i)

	enddo

END PROGRAM neutrongammachain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initialize

	!use timestmp
	use mcnp_params
	use mcnp_random

	!call system_clock(seed)
	!call RN_init_problem( 2, seed, 0_I8, 0_I8, 1 )

END SUBROUTINE initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION TimeInteract(t0)

	use mcnp_random
	use mcnp_params
	use mat_params

	IMPLICIT NONE

	!integer, save :: savestate

	real(R8) :: t0, TimeInteract

	if ( nu .eq. 1) then

	TimeInteract = - (1 / ltot) * log ( rang() )

	else

	!call system_clock(seed)
	!call RN_init_problem( 1, seed, 0_I8, 0_I8, 0 )

	TimeInteract = - (1 / ltot) * log ( rang() )

	endif

END FUNCTION TimeInteract

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

FUNCTION NtrnMult(rnmN)

	use mcnp_random
	use mcnp_params

	IMPLICIT NONE

	integer, parameter :: NtrnMax = 8
	integer :: i

	real(R8) :: rnm, rnmN
	real :: sum = 0
	real :: nubar = 0
	real, dimension( NtrnMax ) :: ProbNu
	real, dimension( NtrnMax ) :: CumNu = 0
	integer(I8) :: NtrnMult

	call system_clock( seed )
	call RN_init_problem( 1, seed, 0_I8, 0_I8, 0 )

	!ProbNu = (/ 0.2, 0.3, 0.2, 0.1, 0.1, 0.05, 0.05 /)

	open( unit = 10, file = "ntrn.mult" )

	read( unit = 10, FMT = * ) ProbNu

	close( unit = 10 )

	do i=1,NtrnMax

		CumNu(i) = sum + ProbNu(i)
		sum = CumNu(i)

	enddo

	!print *, "CumNu: ", CumNu

	do i=1,NtrnMax

		nubar = nubar + ProbNu(i)*(i-1)

	enddo

	!print *, nubar

	!print *, rnm

	do i=1,NtrnMax

		if ( rnmN .lt. CumNu(1) ) then

			NtrnMult = 0

		else if ( ( rnmN .gt. CumNu(i) ) .and. ( rnmN .le. CumNu(i+1) ) ) then

			NtrnMult = i 

		else if ( rnmN .gt. CumNu(NtrnMax) ) then 

			NtrnMult = NtrnMax - 1

		else

			continue

		endif

	enddo

	sum = 0
	nubar = 0

	NtrnMult = 2

END FUNCTION NtrnMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION GammaMult(rnm)

	use mcnp_random
	use mcnp_params

	IMPLICIT NONE

	integer, parameter :: GammaMax = 24
	integer :: i

	real(R8) :: rnm
	real(R8) :: sum = 0
	real(R8) :: mubar = 0
	real(R8), dimension( GammaMax ) :: ProbMu
	real(R8), dimension( GammaMax ) :: CumMu
	integer(I8) :: GammaMult

	!call system_clock( seed )
	!call RN_init_problem( 1, seed, 0_I8, 0_I8, 0 )

	!ProbMu = (/ 0.2, 0.3, 0.2, 0.1, 0.1, 0.05, 0.05 /)

	open( unit = 11, file = "gamma.mult" )

	read( unit = 11, FMT = * ) ProbMu

	close( unit = 11 )

	do i=1,GammaMax

		CumMu(i) = sum + ProbMu(i)
		sum = CumMu(i)

	enddo

	print *, "CumMu: ", CumMu

	do i=1,GammaMax

		mubar = mubar + ProbMu(i)*(i-1)

	enddo

	!rnm = rang()

	do i=1,GammaMax

		if ( rnm .lt. CumMu(1) ) then

			GammaMult = 0

		else if ( ( rnm .gt. CumMu(i) ) .and. ( rnm .le. CumMu(i+1) ) ) then

			GammaMult = i 

		else if ( rnm .gt. CumMu(GammaMax) ) then 

			GammaMult = GammaMax - 1

		else

			continue

		endif

	enddo

	sum = 0
	mubar = 0

END FUNCTION GammaMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION SpotaneousSrcTime(t0,rnm)

	use mcnp_random
	use mcnp_params

	IMPLICIT NONE

	integer, parameter :: S = 1
	real(R8), parameter :: tmax = 1

	real(R8) :: rnm, t0, SpotaneousSrcTime

	!call system_clock( seed )
	!call RN_init_problem( 1, seed, 0_I8, 0_I8, 0 )

	!rnm = rang()

	SpotaneousSrcTime = t0 + rnm * S * tmax

END FUNCTION SpotaneousSrcTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION SpotaneousNu(rnm)

	use mcnp_random
	use mcnp_params

	IMPLICIT NONE

	integer, parameter :: SpontNtrnMax = 8
	integer :: i

	real(R8) :: rnm, sum 
	real(R8) :: spontnubar = 0
	real(R8), dimension( SpontNtrnMax ) :: ProbSpontNu
	real(R8), dimension( SpontNtrnMax ) :: CumSpontNu
	integer(I8) :: SpotaneousNu

	open( unit = 12, file = "spontntrn.mult")

	read( unit = 12, FMT = * ) ProbSpontNu

	close( unit = 12 )

	do i=1,SpontNtrnMax

		CumSpontNu(i) = sum + ProbSpontNu(i)
		sum = CumSpontNu(i)

	enddo

	!print *, CumSpontNu

	do i=1,SpontNtrnMax

		spontnubar = spontnubar + ProbSpontNu(i)*(i-1)

	enddo

	!print *, spontnubar

	do i=1,SpontNtrnMax

		if ( rnm .lt. CumSpontNu(1) ) then

			SpotaneousNu = 0

		else if ( ( rnm .gt. CumSpontNu(i) ) .and. ( rnm .le. CumSpontNu(i+1) ) ) then

			SpotaneousNu = i 

		else if ( rnm .gt. CumSpontNu(SpontNtrnMax) ) then 

			SpotaneousNu = SpontNtrnMax - 1

		else

			continue

		endif

	enddo

	sum = 0
	spontnubar = 0

	!print *, rnm
	!print *, SpotaneousNu
	
	!SpotaneousNu = 1

END FUNCTION SpotaneousNu

!SUBROUTINE SpontNu(rnm,t0,time,SpontNu)

!	use mcnp_random
!	use mcnp_params

!	IMPLICIT NONE

!	real, intent(IN) :: rnm, t0
!	real, intent(OUT) :: time
!	integer, intent(OUT) :: SpontNu

!	integer, parameter :: S = 1
!	real, parameter :: tmax = 1

!	SpontNu = 1

!	time = t0 + rnm * S * tmax

!END SUBROUTINE SpotaneousSrc
