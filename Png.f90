PROGRAM Png

	IMPLICIT NONE

	INTEGER :: i, j, k, mntidx
	INTEGER, PARAMETER :: ntrnlens = 100
	INTEGER, PARAMETER :: gammalens = 100
    INTEGER, PARAMETER :: N = 50
    INTEGER, PARAMETER :: chains = 1000
	INTEGER, PARAMETER :: neut = 30
	INTEGER, PARAMETER :: gama = 99

	INTEGER, DIMENSION(chains,N+1) :: PnData = 0
	INTEGER, DIMENSION(chains,N+1) :: PgData = 0

	INTEGER, PARAMETER :: mmnt = 4

	!REAL, DIMENSION(neut+1,N+1) :: Pn = 0
	!REAL, DIMENSION(gama+1,N+1) :: Pg = 0

	REAL, DIMENSION(neut+1,N+1,mmnt+1) :: Pn = 0
	REAL, DIMENSION(gama+1,N+1,mmnt+1) :: Pg = 0

	REAL :: t0, tf, dt
	REAL, DIMENSION(N+1) :: time

	!CHARACTER(LEN = 10) :: fid = "PnMmnt"
	!CHARACTER(LEN = 10) :: file

	CHARACTER(LEN = 6), DIMENSION(mmnt+1) :: fid = 'PnMmnt'
	CHARACTER(LEN = 6), DIMENSION(mmnt+1) :: fidg = 'PgMmnt'
	CHARACTER(LEN = 25), DIMENSION(mmnt+1) :: filename, filenamegamma
	CHARACTER(LEN = 1) :: filenum

	!Initializing file names

	do i = 1,mmnt+1

		write(filenum,'(i1)') i-1
		filename(i) = fid(i)//filenum//'.txt'
		filenamegamma(i) = fidg(i)//filenum//'.txt'
		!print *, filename
		print *, filenamegamma(i)

	enddo

	tf = 20
	t0 = 0

	dt = (tf-t0)/N

	!print *, 0**0

	do i = 1,N+1

		time(i) = t0 + (i-1)*dt

	enddo

	!print *, time

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

	!do mntidx = 1, mmnt+1

	do i = 1, neut+1

		do j = 1, N+1

			do k = 1, chains

				if ( PnData(k,j) .eq. ( i - 1 ) ) then

					Pn(i,j,:) = (Pn(i,j,:) + 1)!*((i-1)**(mntidx-1))
					!Pn(i,j,mntidx) = Pn(i,j,mntidx) + Pn(i,j,mntidx)*(mntidx**(mntidx-1))

				endif

			enddo

		enddo

	enddo

	Pn(:,:,:) = Pn(:,:,:)/chains

	do mntidx = 1,mmnt+1

		do i = 1, neut+1

		Pn(i,:,mntidx) = Pn(i,:,mntidx)*(i-1)**(mntidx-1)

		enddo

	enddo	

	!enddo

	!write(fid,'(I3)') chains
	!file = 'out'//fid//
	!file = trim(adjustl(file))

	!print *, file

	!open( unit = 3, file = "ProbN.txt" )

	do i = 1,neut+1

		!write(3,*), Pn(i,:)/chains

	enddo

	do j = 1,mmnt+1

		open( unit = 10+j, file = filename(j) )

		do i = 1,neut+1

			write(10+j,*), Pn(i,:,j)!/chains!/(chains*neut**(j-1))  !/(chains*(i**(j-1)))

		enddo

		close( unit = 10+j ) 

	enddo

	!close( unit = 3 )

	do i = 1,gama+1

		do j =1,N+1

			do k = 1,chains

				if ( PgData(k,j) .eq. ( i - 1) ) then

					Pg(i,j,:) = Pg(i,j,:) + 1

				endif

			enddo

		enddo

	enddo

	Pg(:,:,:) = Pg(:,:,:)/chains

	do mntidx = 1,mmnt+1

		do i = 1, gama+1

		Pg(i,:,mntidx) = Pg(i,:,mntidx)*(i-1)**(mntidx-1)

		enddo

	enddo

	do j = 1,mmnt+1

		open( unit = 20+j, file = filenamegamma(j) )

		do i = 1,gama+1

			write(20+j,*), Pg(i,:,j)!/chains!/(chains*neut**(j-1))  !/(chains*(i**(j-1)))

		enddo

		close( unit = 20+j ) 

	enddo	

	!open( unit = 4, file = "ProbG.txt" )

	!do i = 1,gama+1

	!	write(4,*), Pg(i,:)/chains

	!enddo

END PROGRAM Png
