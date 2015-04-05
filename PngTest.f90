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

	INTEGER :: batchidx = 0

	REAL, DIMENSION(neut+1,N+1,batch) :: Pn = 0
	REAL, DIMENSION(gama+1,N+1,batch) :: Pg = 0
	!REAL, DIMENSION(chains/100,N+1) :: mean = 0
	REAL, DIMENSION(chains*neut,N+1) :: PnMmntData, PgMmntData
	REAL, DIMENSION(neut,N+1,mmnt+1) :: PnMeans = 0!, PgMeans = 0
	REAL, DIMENSION(gama,N+1,mmnt+1) :: PgMeans = 0
	REAL, DIMENSION(batch*neut,N+1,mmnt+1,neut) :: PnMmntMatrix
	REAL, DIMENSION(batch*gama,N+1,mmnt+1,gama) :: PgMmntMatrix

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

	open( unit = 3, file = 'PnMmnt0.txt' )

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

	open( unit = 4, file = 'PnMmnt0.txt' )

	do i = 1, batch*neut

		read( unit = 4, FMT = * ) PnMmntData(i,:)

	enddo

	close( unit = 4 )

	!j = 0

	do i = 1, mmnt+1

		!write(filenum,'(i1)') i
		!filename(i) = fid(i)//filenum//'.txt'

		!open( unit = i+20, file = filename(i) )

	do j = 1, neut

		do k = 1, batch

				!write(i+20,*) PnMmntData(j+30*(k-1),:)*(j-1)**(i)
				PnMmntMatrix(k,:,i,j) = PnMmntData(j+neut*(k-1),:)*(j-1)**(i-1)

			enddo

		enddo

		!close(unit = i+20)

	enddo

	open( unit = 120, file = 'test.file' )
				
	do j = 1, neut

		do k = 1, batch

			write(120,*) PnMmntMatrix(k,:,1,j)

			enddo

		enddo

	!enddo	

	close(unit = 120)



	open( unit = 5, file = 'PgMmnt0.txt' )

	do i = 1, batch

		do j = 1, gama

			do k = 1, N+1

				do l = 1+100*(i-1), 100*i

					if ( PgData(l,k) .eq. (j-1) ) then

						Pg(j,k,i) = Pg(j,k,i) + 1

					endif

				enddo

			enddo

		enddo

		do p = 1, gama

			write(5,*) Pg(p,:,i)/100

		enddo 

	enddo

	close( unit = 5 )

	open( unit = 4, file = 'PgMmnt0.txt')

	do i = 1, batch*gama

		read( unit = 4, FMT = * ) PgMmntData(i,:)

	enddo

	close( unit = 4 )

	! Maybe I can create an 4-D array that will prevent having to read text files back and forth


	do i = 1, mmnt+1

		write(filenum,'(i1)') i
		filename(i) = fidg(i)//filenum//'.txt'

		open( unit = i+20, file = filename(i) )

		do k = 1, batch

			do j = 1, gama

				!write(i+20,*) PgMmntData(j+30*(k-1),:)*(j-1)**(i)
				PgMmntMatrix(k,:,i,j) = PgMmntData(j+gama*(k-1),:)*(j-1)**(i-1)

			enddo

		enddo

		close(unit = i+20)

	enddo

	open( unit = 121, file = 'testgamma.file' )

	do k = 1, batch

		do j = 1, gama

			write(121,*) PgMmntMatrix(k,:,2,j)

		enddo

	enddo

	!print *, batch

	!Perhaps make 4D to 3D array????s

	do i = 1, 1 !mmnt

		do j = 1, 1 !neut

			do k = j,batch*neut,neut !batch

				!print *, k
				PnMeans(j,:,i) = PnMeans(j,:,i) +  PnMmntMatrix(j+(k-1)*neut,:,i,j)
				!print *, PnMeans(j,:,i)
				!PnMeans(j,:,i) = PnMeans(j,:,i) +  PnMmntMatrix(j+(k-1)*neut,:,i,j)

			enddo

		enddo
		
	enddo	

	!print *, k	

	!print *, j+(k-1)*neut
	
	!print *, PnMeans(1,:,1)

	!print *, PnMmntMatrix(11,:,1,1)	

END PROGRAM Png
