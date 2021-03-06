PROGRAM ntrngammadataread

	IMPLICIT NONE

	INTEGER :: i, j, k
	INTEGER, PARAMETER :: N = 20
	INTEGER, PARAMETER :: ntrnlens = 1000
	INTEGER, PARAMETER :: gammalens = 10 * ntrnlens
	INTEGER, PARAMETER :: chains = 10

	REAL, DIMENSION(ntrnlens*chains,2) :: ntrnarr
	REAL, DIMENSION(gammalens*chains,2) :: gammaarr
	REAL, DIMENSION(N+1) :: time
	INTEGER, DIMENSION(chains,N+1) :: ntrntal = 0
	INTEGER, DIMENSION(chains,N+1) :: gammatal = 0

	REAL :: t0, tf, dt

	! Reading neutron and gamma lifetime text files
	open( unit = 1, file = "ntrnlife.txt" )
	!open( unit = 1, file = "ntrnlifedata" )
	
	!$OMP PARALLEL DO
	do i=1,ntrnlens*chains

		read( unit = 1, FMT =  * ) ntrnarr(i,:)

	enddo
	!$OMP END PARALLEL DO

	close( unit = 1 )

	open( unit = 2, file = "gammalife.txt" )
	!open( unit = 2, file = "gammalifedata" )

	!$OMP PARALLEL DO
	do i=1,gammalens*chains

		read( unit = 2, FMT = * ) gammaarr(i,:)

	enddo
	!$OMP END PARALLEL DO

	close( unit = 2 )

	!Time Tally Bins
	t0 = 0
	tf = 20

	dt = (tf - t0)/N

	do i = 1,N+1

		time(i) = t0 + (i-1)*dt

	enddo

	!print *, time

	!Tallying of neutrons

	do k = 1, chains

		do i = 1+(k-1)*ntrnlens, ntrnlens*k

			do j = 1, N+1

				if ( ntrnarr(i,1) .eq. ntrnarr(i,2) ) then

					cycle

				endif

				if ( ( ntrnarr(i,1) .le. time(j) ) .and. ( ntrnarr(i,2) .gt. time(j) ) ) then

					ntrntal(k,j) = ntrntal(k,j) + 1

				endif

			enddo

		enddo

	enddo

	!Tallying of photons

	do k = 1, chains

		do i = 1+(k-1)*gammalens, gammalens*k

			do j = 1, N+1

				if ( gammaarr(i,1) .eq. 0 ) then

					cycle

				endif

				if ( gammaarr(i,1) .le. time(j) ) then ! .and. ( gammaarr(i,2) .gt. time(j) ) ) then

					gammatal(k,j) = gammatal(k,j) +  1

				endif

			enddo

		enddo
		
	enddo		

	
	!open( unit = 2, file = "data.txt")
	open ( unit = 11, file = "ntrnanalysis.txt")
	open ( unit = 12, file = "gammanalysis.txt")

	!print *, tal

	do k = 1, chains

		write(11,*), ntrntal(k,:)

	enddo

	do k = 1, chains

		write(12,*), gammatal(k,:)

	enddo

END PROGRAM ntrngammadataread
