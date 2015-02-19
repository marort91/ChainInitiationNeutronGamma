PROGRAM ntrndataread

	IMPLICIT NONE

	INTEGER :: i, j, k
	INTEGER, PARAMETER :: N = 1000
	INTEGER, PARAMETER :: lens = 50
	INTEGER, PARAMETER :: chains = 1e3

	REAL, DIMENSION(lens*chains,2) :: ntrnarr
	REAL, DIMENSION(N+1) :: time
	INTEGER, DIMENSION(chains,N+1) :: ntrntal = 0

	REAL :: t0, tf, dt

	! Reading neutron and gamma lifetime text files
	open( unit = 1, file = "ntrnlife.txt")

	do i=1,lens*chains

		read( unit = 1, FMT =  * ) ntrnarr(i,:)

	enddo

	!Time Tally Bins

	t0 = 0
	tf = 20

	dt = (tf - t0)/N

	do i = 1,N+1

		time(i) = t0 + (i-1)*dt

	enddo

	!Tallying of neutrons

	do k = 1, chains

		do i = 1+(k-1)*lens, lens*k

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

	
	!open( unit = 2, file = "data.txt")
	open ( unit = 2, file = "ntrnanalysis.txt")

	!print *, tal

	do k =1, chains

		write(2,*), ntrntal(k,:)

	enddo



END PROGRAM ntrndataread
