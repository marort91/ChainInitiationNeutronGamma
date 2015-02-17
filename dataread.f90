PROGRAM dataread

	IMPLICIT NONE

	INTEGER :: i, j, k
	INTEGER, PARAMETER :: N = 1000
	INTEGER, PARAMETER :: lens = 50
	INTEGER, PARAMETER :: chains = 2e4

	REAL, DIMENSION(lens*chains,2) :: arr
	REAL, DIMENSION(N+1) :: time
	INTEGER, DIMENSION(chains,N+1) :: tal = 0

	REAL :: t0, tf, dt

	!open( unit = 1, file = "test.txt")
	!open( unit = 1, file = "mario.txt")
    !open( unit = 1, file = "cut")
    open( unit = 1, file = "ntrnlife.txt")
    !open( unit = 1, file = "lifedata" )

	do i=1,lens*chains

		read( unit = 1, FMT =  * ) arr(i,:)

	enddo

	t0 = 0
	tf = 20

	dt = (tf - t0)/N

	do i = 1,N+1

		time(i) = t0 + (i-1)*dt

	enddo

	do k = 1, chains

		do i = 1+(k-1)*lens, lens*k

			do j = 1, N+1

				if ( arr(i,1) .eq. arr(i,2) ) then

					cycle

				endif

				if ( (arr(i,1) .le. time(j) ) .and. ( arr(i,2) .gt. time(j) ) ) then

					tal(k,j) = tal(k,j) + 1

				endif

			enddo

		enddo

	enddo

	do j = 1,N+1

		!print *, tal(j)

	enddo

	!open( unit = 2, file = "data.txt")
	open ( unit = 2, file = "fortran.txt")

	!print *, tal

	do k =1, chains

		write(2,*), tal(k,:)

	enddo



END PROGRAM dataread
