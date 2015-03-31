PROGRAM Png

	IMPLICIT NONE

	INTEGER :: i, j, k, mntidx, l, m, p, o
	INTEGER, PARAMETER :: ntrnlens = 100
	INTEGER, PARAMETER :: gammalens = 100
    INTEGER, PARAMETER :: N = 20
    INTEGER, PARAMETER :: chains = 10000
	INTEGER, PARAMETER :: neut = 30
	INTEGER, PARAMETER :: gama = 199

	INTEGER, PARAMETER :: batch = chains/100

	INTEGER, DIMENSION(chains,N+1) :: PnData = 0
	INTEGER, DIMENSION(chains,N+1) :: PgData = 0

	INTEGER, PARAMETER :: mmnt = 4

	!REAL, DIMENSION(neut+1,N+1,mmnt+1) :: Pn = 0
	REAL, DIMENSION(neut+1,N+1,batch) :: Pn = 0
	!REAL, DIMENSION(gama+1,N+1,mmnt+1) :: Pg = 0
	REAL, DIMENSION(chains/100,N+1) :: mean = 0
	REAL, DIMENSION(chains*neut,N+1) :: PnMmntData, PgMmntData

	REAL :: t0, tf, dt
	REAL, DIMENSION(N+1) :: time

	CHARACTER(LEN = 6), DIMENSION(mmnt+1) :: fid = 'PnMmnt'
	CHARACTER(LEN = 6), DIMENSION(mmnt+1) :: fidg = 'PgMmnt'
	CHARACTER(LEN = 25), DIMENSION(mmnt+1) :: filename, filenamegamma
	CHARACTER(LEN = 1) :: filenum

	!Initializing file names

	do i = 1,mmnt+1

		write(filenum,'(i1)') i-1
		filename(i) = fid(i)//filenum//'.txt'
		filenamegamma(i) = fidg(i)//filenum//'.txt'

	enddo

	tf = 20
	t0 = 0

	dt = (tf-t0)/N

	!print *, 0**0

	do i = 1,N+1

		time(i) = t0 + (i-1)*dt

	enddo

	open( unit = 1, file = "ntrntal.txt" )

	do i = 1,chains

		read( unit = 1, FMT = * ) PnData(i,:)

	enddo

	close( unit = 1 )

	open( unit = 2, file = "gammatal.txt" )

	do i = 1,chains

		read( unit = 2, FMT = * ) PgData(i,:)

	enddo

	close( unit = 2 )

	open( unit = 3, file = 'ntrnstat.test' )

	do i = 1, batch

		do j = 1, neut

			do k = 1, N+1

				do l = 1+100*(i-1), 100*i

					if ( PnData(l,k) .eq. (j-1) ) then

						Pn(j,k,i) = Pn(j,k,i) + 1

					endif

				enddo

			enddo

		enddo

		do p = 1, neut

			write(3,*) Pn(p,:,i)/100

		enddo 

	enddo

	close( unit = 3 )

	open( unit = 4, file = 'ntrnstat.test' )

	do i = 1, batch*neut

		read( unit = 4, FMT = * ) PnMmntData(i,:)
		!print *, PnMmntData(i,:)

	enddo

	close( unit = 4 )

	for i = 1, mmnt

		for j = 1, neut

			k = 

	!print *, batch*neut

	!Implement moments here, save to file names as done before.	
	
END PROGRAM Png
