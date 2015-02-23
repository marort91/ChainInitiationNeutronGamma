PROGRAM Png

	IMPLICIT NONE

	INTEGER :: i, j, k
	INTEGER, PARAMETER :: ntrnlens = 50
	INTEGER, PARAMETER :: gammalens = 100
    INTEGER, PARAMETER :: N = 50
    INTEGER, PARAMETER :: chains = 1000
	INTEGER, PARAMETER :: neut = 30
	INTEGER, PARAMETER :: gama = 99

	INTEGER, DIMENSION(chains,N+1) :: PnData = 0
	INTEGER, DIMENSION(chains,N+1) :: PgData = 0

	REAL, DIMENSION(neut+1,N+1) :: Pn = 0
	REAL, DIMENSION(gama+1,N+1) :: Pg = 0

	open( unit = 1, file = "ntrnanalysis.txt" )

	do i = 1,chains

		read( unit = 1, FMT = * ) PnData(i,:)

	enddo

	close( unit = 1 )

	open( unit = 2, file = "gammanalysis.txt")

	do i = 1,chains

		read( unit = 2, FMT = * ) PgData(i,:)

	enddo

	close( unit = 2 )

	!print *, PnData
	!print *, PgData

	do i = 1,neut+1

		do j = 1,N+1

			do k = 1, chains

				if ( PnData(k,j) .eq. ( i - 1 ) ) then

					Pn(i,j) = Pn(i,j) + 1

				endif

			enddo

		enddo

	enddo

	open( unit = 3, file = "ProbN.txt" )

	do i = 1,neut+1

		write(3,*), Pn(i,:)/chains

	enddo

	close( unit = 3 )

	do i = 1,gama+1

		do j =1,N+1

			do k = 1,chains

				if ( PgData(k,j) .eq. ( i - 1) ) then

					Pg(i,j) = Pg(i,j) + 1

				endif

			enddo

		enddo

	enddo

	open( unit = 4, file = "ProbG.txt" )

	do i = 1,gama+1

		write(4,*), Pg(i,:)/chains

	enddo

END PROGRAM Png
