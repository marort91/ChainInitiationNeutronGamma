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

	real, parameter :: ntrnlfission = 0.01
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

    integer :: nu, mu

END MODULE mat_params

MODULE time_params

	integer, parameter :: loop = 10000
	integer, parameter :: N = 25
	real, parameter :: ti = 0, tf = 50, dt = ( tf - ti ) / N 
	real, dimension(N+1) :: time

END MODULE time_params

MODULE prob_flags

	integer, parameter :: ICNtrnFlag = 0
	integer, parameter :: fissflag = 1
	integer, parameter :: spontmuflag = 0
	integer, parameter :: expvaloutcome = 0 

END MODULE prob_flags

PROGRAM neutrongammachain

	use mcnp_random
	use mcnp_params
	use mat_params
	use time_params
	use prob_flags
	use omp_lib

	real(R8) :: rnmN, rnmG, rnmLeakNtrn
	real(R8), dimension(loop) :: rnmSp, rnmSpNu, rnmSpMu
	integer(I8), dimension(loop) :: spNu, spMu

	real(R8) :: NtrnTimeInteract, SpotaneousSrcTime, GammaTimeInteract
	real(R8), dimension(:), allocatable :: t0, tsp
	integer(I8) :: i, j, k, NtrnMult, GammaMult, SpontNu, SpontMu
	!integer(I8) :: ICNtrnFlag, fissflag, spontmuflag
	integer(I8), dimension(:), allocatable :: nidx, gidx
	integer, parameter :: branchlens = 1000
	
	real(R8), dimension(:,:,:,:), allocatable :: ntrntime, gammatime

	integer, dimension(loop,N+1) :: ntrntal = 0
	integer, dimension(loop,N+1) :: gammatal = 0

	integer(I8) :: ntrnl, gammal, spntl, spntgl
	real(R8), dimension(:), allocatable :: CumNu, CumMu, CumSpontNu, CumSpontMu
	real(R8) :: NuBar, MuBar, SpontNuBar, SpontMuBar

	call datalens(ntrnl, gammal, spntl, spntgl)

	allocate(CumNu(ntrnl))
	allocate(CumMu(gammal))
	allocate(CumSpontNu(spntl))
	allocate(CumSpontMu(spntgl))

	allocate(ntrntime(1,2,branchlens,loop))
	allocate(gammatime(1,2,10*branchlens,loop))

	allocate(nidx(loop))
	allocate(gidx(loop))

	allocate(t0(loop))
	allocate(tsp(loop))

	ntrntime(:,:,:,:) = 0
	gammatime(:,:,:,:) = 0

	nidx(:) = 0
	gidx(:) = 0

	t0(:) = 0
	tsp(:) = 0

	call CumDistFcnGen(ntrnl, gammal, spntl, spntgl, CumNu, CumMu, CumSpontNu, CumSpontMu, NuBar, MuBar, SpontNuBar, SpontMuBar )

	!call omp_set_num_threads(8)

	!print *, NuBar, MuBar, SpontNuBar, SpontMuBar

	!call sleep(3600)
	
	open( unit = 1, file = "ntrnlifedata", position = "append", action = "write")
	open( unit = 2, file = "gammalifedata", position = "append", action = "write")
	open( unit = 3, file = "ntrnfissiondata", position = "append", action = "write")
	open( unit = 4, file = "gammafissiondata", position = "append", action = "write")

	!dt = (tf - t0)/N

	do i = 1,N+1

		time(i) = ti + (i-1)*dt
		!print *, time(i)

	enddo

	open( unit = 7, file = 'ntrntal.txt' )
	open( unit = 8, file = 'gammatal.txt')

	call initialize

	if ( ICNtrnFlag .eq. 1 ) then

		nidx(:) = 1

	else

		nidx(:) = 0

	endif

	if ( spontmuflag .eq. 0 ) then

		gidx(:) = 1

	else
		
		gidx(:) = 0

	endif
	
	!$OMP PARALLEL DO
	do k = 1, loop

	print *, k 

	

	do i=1,branchlens

			if ( fissflag .eq. 1 ) then

			!t0(k) = t0(k) + 1

			!if ( i .eq. 1 ) then

			!rnmSpNu(k) = rang()
			!spNu(k) = SpontNu(spntl,CumSpontNu,rnmSpNu(k))

			!ntrntime(1,1,nidx(k)+1:nidx(k)+spNu(k),k) = t0(k)+0.6321/1000
			!nidx(k) = nidx(k) + spNu(k)

			!else

			rnmSp(k) = rang()
			tsp(k) = SpotaneousSrcTime(t0(k),rnmSp(k))

			rnmSpNu(k) = rang()
			spNu(k) = SpontNu(spntl,CumSpontNu,rnmSpNu(k),SpontNuBar)

			ntrntime(1,1,nidx(k)+1:nidx(k)+spNu(k),k) = tsp(k)
		
			nidx(k) = nidx(k) + spNu(k)

			!endif

			if ( spontmuflag .eq. 1 ) then
			
				rnmSpMu(k) = rang()
				spMu(k) = SpontMu(spntgl,CumSpontMu,rnmSpMu(k),SpontMuBar)

				gammatime(1,1,gidx(k)+1:gidx(k)+spMu(k),k) = tsp(k)

				gidx(k) = gidx(k) + spMu(k)

			endif

			t0(k) = tsp(k)
			!print *, t0(k)

		endif

		if ( ( nidx(k) .eq. 0 ) .and. ( fissflag .ne. 1 ) ) exit

		if ( ( ntrntime(1,1,i,k).eq.0_R8) .and. ( i.gt.1 ) ) exit
		
		ntrntime(1,2,i,k) = ntrntime(1,1,i,k) + NtrnTimeInteract(ntrntime(1,1,i,k))

		rnmLeakNtrn = rang()

		if ( rnmLeakNtrn .le. ntrnPleakage ) then

			cycle

		else

			if (rang() .le. ntrnPfiss) then

				rnmN = rang()
				rnmG = rang()

				nu = NtrnMult(ntrnl,CumNu,rnmN,NuBar)
				mu = GammaMult(gammal,CumMu,rnmG,MuBar)

				write( 3, * ), nu
				write( 4, * ), mu

				ntrntime(1,1,nidx(k)+1:nidx(k)+nu,k) = ntrntime(1,2,i,k)
				gammatime(1,1,gidx(k):gidx(k)+mu-1,k) = ntrntime(1,2,i,k)
			
				nidx(k) = nidx(k) + nu
				gidx(k) = gidx(k) + mu

				if ( nu .eq. 0 ) then

					continue

				else

				endif

			endif

		endif

		if ( nidx(k) + 1 .gt. branchlens ) exit
		if ( gidx(k) + 1 .gt. 10*branchlens ) exit

	enddo

	do i=1,branchlens

		if ( (ntrntime(1,1,i,k) .ne. 0) .and. (ntrntime(1,2,i,k) .eq. 0) ) then

			ntrntime(1,2,i,k) = ntrntime(1,1,i,k) + NtrnTimeInteract(ntrntime(1,1,i,k))

		endif

	enddo

	!do k = 1, loop
	do i = 1, branchlens

		write(1,*) ntrntime(1,:,i,k)

	enddo
	!enddo

	do i=1,10*branchlens

		if ( (gammatime(1,1,i,k) .ne. 0 ) .and. ( gammatime(1,2,i,k) .eq. 0 ) ) then

			gammatime(1,2,i,k) = gammatime(1,1,i,k) + GammaTimeInteract(gammatime(1,1,i,k))

		endif

	enddo

	do i = 1, branchlens

		do j = 1, N+1

			if ( ntrntime(1,1,i,k) .eq. ntrntime(1,2,i,k) ) then

					cycle

			endif

			if ( ( ntrntime(1,1,i,k) .le. time(j) ) .and. ( ntrntime(1,2,i,k) .gt. time(j) ) ) then

					ntrntal(k,j) = ntrntal(k,j) + 1

			endif

		enddo

	enddo

	do i = 1, branchlens*10

		do j = 1, N+1

				if ( gammatime(i,i,1,k) .eq. 0 ) then

					cycle

				endif

				if ( gammatime(1,1,i,k) .le. time(j) ) then ! .and. ( gammaarr(i,2) .gt. time(j) ) ) then

					gammatal(k,j) = gammatal(k,j) +  1

				endif

			enddo

		enddo	

	enddo
	!$OMP END PARALLEL DO

	do k = 1, loop

		write( 7, * ) ntrntal(k,:)
		write( 8, * ) gammatal(k,:)

	enddo

END PROGRAM neutrongammachain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initialize

	use mcnp_params
	use mcnp_random

	IMPLICIT NONE

	call system_clock(seed)
	call RN_init_problem( 2, seed, 0_I8, 0_I8, 1 )

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

	real(R8), parameter :: S = 1 !1
	real(R8), parameter :: tmax = 2

	real(R8) :: rnm, t0, SpotaneousSrcTime

	open( unit = 15, file = 'SpontNuEmission.txt' )

	SpotaneousSrcTime = t0 + rnm / S * tmax !/ S * tmax

	!write(15,*) rnm / S * tmax

END FUNCTION SpotaneousSrcTime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE datalens( nlines, glines, spntlines, spntglines )

	use mcnp_params

	IMPLICIT NONE

	integer :: io
	integer(I8), intent(out) :: nlines, glines, spntlines, spntglines

	nlines = 0
	glines = 0
	spntlines = 0
	spntglines = 0

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

	open( unit = 96, file = 'spontgamma.mult', status = 'old', action = 'read' )

	do

		read( 96, *, iostat = io )

		if ( io .ne. 0 ) exit

		spntglines = spntglines + 1

	enddo

	rewind(96)	

	close( 96 )
	close( 97 )
	close( 98 )
	close( 99 )

END SUBROUTINE datalens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CumDistFcnGen( nlines, glines, spntlines, spntglines, CumNu, CumMu, CumSpontNu, CumSpontMu, NuBar, MuBar, &
						  SpontNuBar, SpontMuBar)

	use mcnp_params

	IMPLICIT NONE

	integer(I8), intent(in) :: nlines, glines, spntlines, spntglines
	real(R8), dimension(nlines), intent(out) :: CumNu
	real(R8), dimension(glines), intent(out) :: CumMu
	real(R8), dimension(spntlines), intent(out) :: CumSpontNu
	real(R8), dimension(spntglines), intent(out) :: CumSpontMu
	real(R8), intent(out) :: NuBar, MuBar, SpontNuBar, SpontMuBar 

	integer :: i
	real :: sumn = 0, sumg = 0, sumspont = 0, sumspontg = 0
	real, dimension(nlines) :: ProbNu
	real, dimension(glines) :: ProbMu
	real, dimension(spntlines) :: ProbSpontNu
	real, dimension(spntglines) :: ProbSpontMu

	open( unit = 96, file = 'spontgamma.mult', status = 'old', action = 'read')
	open( unit = 97, file = 'ntrn.mult', status = 'old', action = 'read' )
	open( unit = 98, file = 'gamma.mult', status = 'old', action = 'read' )
	open( unit = 99, file = 'spontntrn.mult', status = 'old', action = 'read' )

	read( 96, FMT = * ) ProbSpontMu
	read( 97, FMT = * ) ProbNu
	read( 98, FMT = * ) ProbMu
	read( 99, FMT = * ) ProbSpontNu

	NuBar = 0
	MuBar = 0
	SpontNuBar = 0
	SpontMuBar = 0

	do i = 1,spntglines

		CumSpontMu(i) = sumspontg + ProbSpontMu(i)
		sumspontg = CumSpontMu(i)
		SpontMuBar = SpontMuBar + ProbSpontMu(i)*( i - 1 )

	enddo	

	do i = 1,nlines

		CumNu(i) = sumn + ProbNu(i)
		sumn = CumNu(i)
		NuBar = NuBar + ProbNu(i)*( i - 1 )

	enddo

	do i = 1, glines

		CumMu(i) = sumg + ProbMu(i)
		sumg = CumMu(i)
		MuBar = MuBar + ProbMu(i)*( i - 1 )

	enddo

	do i = 1, spntlines

		CumSpontNu(i) = sumspont + ProbSpontNu(i)
		sumspont = CumSpontNu(i)
		SpontNuBar = SpontNuBar + ProbSpontNu(i)*( i - 1 )

	enddo

	close( 96 )
	close( 97 )
	close( 98 )
	close( 99 )

END SUBROUTINE CumDistFcnGen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION NtrnMult(N,CumNu,rnm,NuBar)

	use mcnp_params
	use prob_flags

	IMPLICIT NONE

	integer(I8) :: N, NtrnMult, i
	real(R8), dimension(N) :: CumNu
	real(R8) :: rnm, NuBar

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

	NtrnMult = 2

END FUNCTION NtrnMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION GammaMult(G,CumMu,rnm,MuBar)

	use mcnp_params
	use prob_flags

	IMPLICIT NONE

	integer(I8) :: G, GammaMult, i
	real(R8), dimension(G) :: CumMu
	real(R8) :: rnm,MuBar

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

	!GammaMult = 2

END FUNCTION GammaMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION SpontNu(S,CumSpontNu,rnm,SpontNuBar)

	use mcnp_params
	use prob_flags

	IMPLICIT NONE

	integer(I8) :: S, SpontNu, i
	real(R8), dimension(S) :: CumSpontNu
	real(R8) :: rnm, SpontNuBar

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

	SpontNu = 1

	open( unit = 16, file = 'SpontNuEmit.txt' )

	write(16,*) SpontNu

END FUNCTION SpontNu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION SpontMu(S,CumSpontMu,rnm,SpontMuBar)

	use mcnp_params
	use prob_flags

	IMPLICIT NONE

	integer(I8) :: S, SpontMu, i
	real(R8), dimension(S) :: CumSpontMu
	real(R8) :: rnm, SpontMuBar

	do i=1,S

		if ( rnm .lt. CumSpontMu(1) ) then

			SpontMu = 0

		else if ( ( rnm .gt. CumSpontMu(i) ) .and. ( rnm .le. CumSpontMu(i+1) ) ) then

			SpontMu = i 

		else if ( rnm .gt. CumSpontMu(S) ) then 

			SpontMu = S - 1

		else

			continue

		endif

	enddo

	open( unit = 13, file = 'SpontMuEmit.txt' )
	write(13,*) SpontMu

END FUNCTION SpontMu	
