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
!! 																													 !!
!! 03/10/2015																										 !!
!! -Multiplicity data read currently being made into own subroutine. Prevents need to read data over branch loops.   !!
!!  Currently having trouble reading all values on one line. Might have to put each value on own line. To be tested. !!
!!  Also, should try to implement sed command in bash script that is independent of line numbers.                    !!
!!																													 !!
!! 03/12/2015																										 !!
!! -Implemented subroutine that reads neutron and gamma multiplicity data outside neutron chain reaction loop.       !!
!! -Need to use MCNP parameters in subroutines and functions to keep variable types consistent. "use mcnp_params"    !!
!!  should be used inside subroutine and variable initializations previously implemented used.                       !!
!! -Calculations of nubar and mubar should be done by subroutine in order to check for k-effective value of problem. !!
!! -Sed editing of files not dependent on line number not yet implemented, be careful!!!!!                           !!
!!																												     !!
!! 03/13/2015																										 !!
!! -Finally stamped out segmentation faults for close to critical problems. Gamma and neutron indices are checked    !!
!!  and logic statements prevent accessing non-allocated array locations.											 !!
!! -Noticed spontaneous fission driven problems overestimated number of neutrons present. Checked average emission   !!
!!  time and found that emission time was half than expected. Average emission time must be doubled I believe. Will  !!
!!  check for different source strengths.																			 !!
!!																													 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE mcnp_params

	integer, parameter :: R8 = selected_real_kind(15,307)
    integer, parameter :: I8 = selected_int_kind(18)

    integer(I8) :: seed = 1

END MODULE mcnp_params
	
MODULE mat_params

	real, parameter :: ntrnlfission = 0.1
    real, parameter :: ntrnlcapture = 1.0 - ntrnlfission
    real, parameter :: ntrnltot = ntrnlfission + ntrnlcapture

    real, parameter :: ntrnPcap = ntrnlcapture / ntrnltot
    real, parameter :: ntrnPfiss = ntrnlfission / ntrnltot
	real, parameter :: ntrnPleakage = 0.0

	real, parameter :: gammalcapture = 0.0
	real, parameter :: gammaleakage = 1.0 - gammalcapture
	real, parameter :: gammaltot = gammalcapture + gammaleakage

	real, parameter :: gammaPcap = gammalcapture / gammaltot
	real, parameter :: gammaPleak = gammaleakage / gammaltot

    integer :: nu, mu, spNu

    !Binary Fission Assumption
    !integer, parameter :: nu = 2

END MODULE mat_params

PROGRAM neutrongammachain

	use mcnp_random
	use mcnp_params
	use mat_params

	real(R8) :: dt, rnmN, rnmG, rnmSp, rnmSpNu, rnmLeakNtrn

	real(R8) :: t0, tf, NtrnTimeInteract, SpotaneousSrcTime, GammaTimeInteract,tsp
	integer(I8) :: i, j, NtrnMult, GammaMult, SpontNu
	integer(I8) :: ICNtrnFlag, fissflag, nidx, gidx

	integer, parameter :: branchlens = 1000
	real(R8), dimension(1,2,branchlens) :: ntrntime = 0
	real(R8), dimension(1,2,branchlens*5) :: gammatime = 0

	!Experimental implementation of new data reader subroutine and functions
	integer(I8) :: ntrnl, gammal, spntl
	real(R8), dimension(:), allocatable :: CumNu, CumMu, CumSpontNu

	call datalens(ntrnl, gammal, spntl)

	allocate(CumNu(ntrnl))
	allocate(CumMu(gammal))
	allocate(CumSpontNu(spntl))

	call CumDistFcnGen(ntrnl, gammal, spntl, CumNu, CumMu, CumSpontNu)

	!Initialize MCNP Random Number Generator
	call system_clock(seed)
	call RN_init_problem( 2, seed, 0_I8, 0_I8, 1 )
	
	open( unit = 1, file = "ntrnlifedata", position = "append", action = "write")
	open( unit = 2, file = "gammalifedata", position = "append", action = "write")
	open( unit = 3, file = "ntrnfissiondata", position = "append", action = "write")
	open( unit = 4, file = "gammafissiondata", position = "append", action = "write")
	
	t0 = 0
	
	ICNtrnFlag = 1
	fissflag = 0

	if ( ICNtrnFlag .eq. 1 ) then

		nidx = 1
		gidx = 1

		else

		nidx = 0
		gidx = 1

	endif

	do i=1,branchlens

		if ( fissflag .eq. 1 ) then

			rnmSp = rang()
			tsp = SpotaneousSrcTime(t0,rnmSp)

			rnmSpNu = rang()
			spNu = SpontNu(spntl,CumSpontNu,rnmSpNu)

			ntrntime(1,1,nidx+1:nidx+spNu) = tsp
		
			nidx = nidx + spNu 
			t0 = tsp

		endif

		if ( ( nidx .eq. 0 ) .and. ( fissflag .ne. 1 ) ) exit

		if ( ( ntrntime(1,1,i).eq.0_R8) .and. ( i.gt.1 ) ) exit
		
		ntrntime(1,2,i) = ntrntime(1,1,i) + NtrnTimeInteract(ntrntime(1,1,i))

		rnmLeakNtrn = rang()

		if ( rnmLeakNtrn .le. ntrnPleakage ) then

			cycle

		else

			if (rang() .le. ntrnPfiss) then

				rnmN = rang()
				rnmG = rang()

				!nu = NtrnMultNew(rnmN)
				!mu = GammaMultNew(rnmG)

				nu = NtrnMult(ntrnl,CumNu,rnmN)
				mu = GammaMult(gammal,CumMu,rnmG)

				write( 3, * ), nu
				write( 4, * ), mu

				ntrntime(1,1,nidx+1:nidx+nu) = ntrntime(1,2,i)
				gammatime(1,1,gidx:gidx+mu-1) = ntrntime(1,2,i)
			
				nidx = nidx + nu
				gidx = gidx + mu

				if ( nu .eq. 0 ) then

					continue

				else

				endif

			endif

		endif

		if ( nidx + 1 .gt. branchlens ) exit
		if ( gidx + 1 .gt. 10*branchlens ) exit

	enddo

	do i=1,branchlens

		if ( (ntrntime(1,1,i) .ne. 0) .and. (ntrntime(1,2,i) .eq. 0) ) then

			ntrntime(1,2,i) = ntrntime(1,1,i) + NtrnTimeInteract(ntrntime(1,1,i))

		endif

	enddo

	do i=1,5*branchlens

		if ( (gammatime(1,1,i) .ne. 0 ) .and. ( gammatime(1,2,i) .eq. 0 ) ) then

			gammatime(1,2,i) = gammatime(1,1,i) + GammaTimeInteract(gammatime(1,1,i))

		endif

	enddo	

	do i=1,branchlens
		
		write( 1, * ), ntrntime(1,1,i), ntrntime(1,2,i)

	enddo

	do i=1,5*branchlens

		write( 2, * ), gammatime(1,1,i), gammatime(1,2,i)

	enddo

	!open( unit = 17, file = 'ntrnlife.txt')

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

FUNCTION NtrnTimeInteract(t0)

	use mcnp_random
	use mcnp_params
	use mat_params

	IMPLICIT NONE

	real(R8) :: t0, NtrnTimeInteract

	NtrnTimeInteract = - (1 / ntrnltot ) * log ( rang() )

END FUNCTION NtrnTimeInteract

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION GammaTimeInteract(t0)

	use mcnp_random
	use mcnp_params
	use mat_params

	IMPLICIT NONE

	real(R8) :: t0, GammaTimeInteract

	GammaTimeInteract = - (1 / gammaltot ) * log ( rang() )

END FUNCTION GammaTimeInteract

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION SpotaneousSrcTime(t0,rnm)

	use mcnp_random
	use mcnp_params

	IMPLICIT NONE

	real(R8), parameter :: S = 1
	real(R8), parameter :: tmax = 2

	real(R8) :: rnm, t0, SpotaneousSrcTime

	!call system_clock( seed )
	!call RN_init_problem( 1, seed, 0_I8, 0_I8, 0 )

	!rnm = rang()

	open( unit = 15, file = 'SpontNuEmission.txt' )

	!write(15,*) tmax/2 !rnm/S*tmax !*2

	SpotaneousSrcTime = t0 + rnm / S * tmax ! *2
	!SpotaneousSrcTime = t0 + tmax

END FUNCTION SpotaneousSrcTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE datalens( nlines, glines, spntlines )

	use mcnp_params

	IMPLICIT NONE

	integer :: io
	integer(I8), intent(out) :: nlines, glines, spntlines

	nlines = 0
	glines = 0
	spntlines = 0

	open( unit = 97, file = 'ntrn.mult', status = 'old', action = 'read' )

	do
		read( 97, *, iostat = io )

		if ( io .ne. 0 ) exit

		nlines = nlines + 1

	enddo

	rewind(97)

	open( unit = 98, file = 'gamma.mult', status = 'old', action = 'read' )

	do

		read( 98, *, iostat = io )

		if ( io .ne. 0 ) exit

		glines = glines + 1

	enddo

	rewind(98)

	open( unit = 99, file = 'spontntrn.mult', status = 'old', action = 'read' )

	do

		read( 99, *, iostat = io )

		if ( io .ne. 0 ) exit

		spntlines = spntlines + 1

	enddo

	rewind(99)

	close( 97 )
	close( 98 )
	close( 99 )

END SUBROUTINE datalens

SUBROUTINE CumDistFcnGen( nlines, glines, spntlines, CumNu, CumMu, CumSpontNu )

	use mcnp_params

	IMPLICIT NONE

	integer(I8), intent(in) :: nlines, glines, spntlines
	real(R8), dimension(nlines), intent(out) :: CumNu
	real(R8), dimension(glines), intent(out) :: CumMu
	real(R8), dimension(spntlines), intent(out) :: CumSpontNu

	integer :: i
	real :: sumn = 0, sumg = 0, sumspont = 0
	real, dimension(nlines) :: ProbNu
	real, dimension(glines) :: ProbMu
	real, dimension(spntlines) :: ProbSpontNu

	open( unit = 97, file = 'ntrn.mult', status = 'old', action = 'read' )
	open( unit = 98, file = 'gamma.mult', status = 'old', action = 'read' )
	open( unit = 99, file = 'spontntrn.mult', status = 'old', action = 'read' )


	read( 97, FMT = * ) ProbNu
	read( 98, FMT = * ) ProbMu
	read( 99, FMT = * ) ProbSpontNu

	do i = 1,nlines

		CumNu(i) = sumn + ProbNu(i)
		sumn = CumNu(i)

	enddo

	do i = 1, glines

		CumMu(i) = sumg + ProbMu(i)
		sumg = CumMu(i)

	enddo

	do i = 1, spntlines

		CumSpontNu(i) = sumspont + ProbSpontNu(i)
		sumspont = CumSpontNu(i)

	enddo

	close( 97 )
	close( 98 )
	close( 99 )

END SUBROUTINE CumDistFcnGen

FUNCTION NtrnMult(N,CumNu,rnm)

	use mcnp_params

	IMPLICIT NONE

	integer(I8) :: N, NtrnMult, i
	real(R8), dimension(N) :: CumNu
	real(R8) :: rnm

	do i=1,N

		if ( rnm .lt. CumNu(1) ) then

			NtrnMult = 0

		else if ( ( rnm .gt. CumNu(i) ) .and. ( rnm .le. CumNu(i+1) ) ) then

			NtrnMult = i 

		else if ( rnm .gt. CumNu(N) ) then 

			NtrnMult = N - 1

		else

			continue

		endif

	enddo

	!NtrnMult = 2

END FUNCTION NtrnMult

FUNCTION GammaMult(G,CumMu,rnm)

	use mcnp_params

	IMPLICIT NONE

	integer(I8) :: G, GammaMult, i
	real(R8), dimension(G) :: CumMu
	real(R8) :: rnm

	do i=1,G

		if ( rnm .lt. CumMu(1) ) then

			GammaMult = 0

		else if ( ( rnm .gt. CumMu(i) ) .and. ( rnm .le. CumMu(i+1) ) ) then

			GammaMult = i 

		else if ( rnm .gt. CumMu(G) ) then 

			GammaMult = G - 1

		else

			continue

		endif

	enddo

	!GammaMult = 5

END FUNCTION GammaMult

FUNCTION SpontNu(S,CumSpontNu,rnm)

	use mcnp_params

	IMPLICIT NONE

	integer(I8) :: S, SpontNu, i
	real(R8), dimension(S) :: CumSpontNu
	real(R8) :: rnm

	do i=1,S

		if ( rnm .lt. CumSpontNu(1) ) then

			SpontNu = 0

		else if ( ( rnm .gt. CumSpontNu(i) ) .and. ( rnm .le. CumSpontNu(i+1) ) ) then

			SpontNu = i 

		else if ( rnm .gt. CumSpontNu(S) ) then 

			SpontNu = S - 1

		else

			continue

		endif

	enddo

	!SpontNu = 1

	open( unit = 16, file = 'SpontNuEmit.txt' )

	write(16,*) SpontNu

END FUNCTION SpontNu

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

!SUBROUTINE NtrnGammaDataMult(CumNu,CumMu)

!	IMPLICIT NONE

!	real, dimension(:), intent(out), allocatable :: CumNu, CumMu
!	integer :: lensntrn, lensgamma

!	open( unit = 98, file = 'ntrn.mult', status = 'old', action = 'read' )
!	open( unit = 99, file = 'gamma.mult', status = 'old', action = 'read' )

!	read( 98, * ), lensntrn
!	read( 99, * ), lensgamma

!	allocate( CumNu (lensntrn) ) 
!	allocate( CumMu (lensgamma) )

!END SUBROUTINE

