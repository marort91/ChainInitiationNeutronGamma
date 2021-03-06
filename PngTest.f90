PROGRAM Png

	IMPLICIT NONE

	INTEGER :: i, j, k, mntidx, l, m, p, o
	INTEGER, PARAMETER :: ntrnlens = 100
	INTEGER, PARAMETER :: gammalens = 100
    INTEGER, PARAMETER :: N = 20
    INTEGER, PARAMETER :: chains = 20000
	INTEGER, PARAMETER :: neut = 60
	INTEGER, PARAMETER :: gama = 199
	INTEGER, PARAMETER :: batchint = 100
	!REAL, PARAMETER :: chain = 20000

	INTEGER, PARAMETER :: batch = chains/batchint

	INTEGER, DIMENSION(:,:), ALLOCATABLE :: PnData
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: PgData

	INTEGER, PARAMETER :: mmnt = 4

	INTEGER :: batchidx = 0, idx = 0

	REAL, DIMENSION(neut,N+1,batch) :: Pn = 0
	REAL, DIMENSION(gama,N+1,batch) :: Pg = 0
	REAL, DIMENSION(:,:), ALLOCATABLE :: PnMmntData, PgMmntData
	REAL, DIMENSION(neut,N+1,mmnt+1) :: PnMeans = 0 !, PgMeans = 0
	REAL, DIMENSION(gama,N+1,mmnt+1) :: PgMeans = 0
	!REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PnMmntMatrix
	REAL, DIMENSION(:,:,:), ALLOCATABLE :: PnMmntMatrix
	REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: PgMmntMatrix
	REAL, DIMENSION(:,:), ALLOCATABLE :: NeutData, GamaData
	REAL, DIMENSION(neut,N+1,mmnt+1) :: varneut = 0
	REAL, DIMENSION(gama,N+1,mmnt+1) :: vargama = 0
	REAL, DIMENSION(neut,2*(N+1),mmnt+1) :: PnMeanVar = 0
	REAL, DIMENSION(gama,2*(N+1),mmnt+1) :: PgMeanVar = 0

	REAL :: t0, tf, dt
	REAL, DIMENSION(N+1) :: time

	CHARACTER(LEN = 6), DIMENSION(mmnt+1) :: fid = 'PnMmnt'
	CHARACTER(LEN = 6), DIMENSION(mmnt+1) :: fidg = 'PgMmnt'
	CHARACTER(LEN = 25), DIMENSION(mmnt+1) :: filename, filenamegamma
	CHARACTER(LEN = 1) :: filenum

	print *, "Batch Number: ", batch

	allocate(PnData(chains,N+1))
	allocate(PgData(chains,N+1))
	!allocate(PnMmntMatrix(batch*neut,N+1,mmnt+1,neut))

	allocate(PnMmntMatrix(batch*neut,N+1,mmnt+1))

	allocate(PgMmntMatrix(batch*gama,N+1,mmnt+1,gama))
	allocate(PnMmntData(neut*chains,N+1))
	allocate(PgMmntData(gama*chains,N+1))
	allocate(NeutData(batch*neut*(mmnt+1),N+1))
	allocate(GamaData(batch*gama*(mmnt+1),N+1))

	!Initializing file names

	do i = 1,mmnt+1

		write(filenum,'(i1)') i-1
		filename(i) = fid(i)//filenum//'.mmnt'
		filenamegamma(i) = fidg(i)//filenum//'.mmnt'

	enddo

	tf = 20
	t0 = 0

	dt = (tf-t0)/N

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

	! This is the fuckup right here...

	do i = 1, batch

		do j = 1, 1 !neut

				do k = 1, N + 1

			do l = 1+(i-1)*batchint, batchint*i

		

					if (PnData(l,k) .eq. (j - 1) ) then

						Pn(j,k,i) = Pn(j,k,i) + 1

					endif
					
				enddo
				
			enddo
		
		enddo

		!print *, Pn(j,:,i)

		do p = 1, 1 !neut

			write(3,*) Pn(p,:,i)

		enddo	

		!print *, i 

	enddo					

	call sleep(3600)

	!do i = 1, batch

		!do j = 1, neut

			!do k = 1, N+1



	!do i = 1, batch

		!print *, 1+batch*(i-1), batch*i

		!do j = 1, neut

			!do k = 1, N+1

				!do l = 1+100*(i-1), 100*i
				!do l = 1+batch*(i-1), batch*i

					!if ( PnData(l,k) .eq. (j-1) ) then

						!Pn(j,k,i) = Pn(j,k,i) + 1

					!endif

				!enddo

			!enddo

		!enddo

		!do p = 1, neut

			!write(3,*) Pn(p,:,i)/real(batch)

		!enddo 

	!enddo

	close( unit = 3 )

	!call sleep(3600)

	open( unit = 4, file = 'PnMmnt0.txt' )

	!print *, batch*neut

	do i = 1, batch*neut

		read( unit = 4, FMT = * ) PnMmntData(i,:)

	enddo

	close( unit = 4 )

	!print *, PnMmntData(1,:)

	!call sleep(3600)

	!!! This is where the problem is...

	do i = 1, 1 !mmnt + 1

		do j = 1, 1 !neut

			do k = 1+batch*(j-1), batch*j

				!PnMmntMatrix(k+(j-1)*batch,:,i) = PnMmntData(j+neut*(j-1),:)*(j-1)**(i-1)
				print *, j + neut*(k-1)
				PnMmntMatrix(k,:,i) = PnMmntData(j+neut*(k-1),:)*(j-1)**(i-1)

			enddo

		enddo
		
	enddo		

	!do i = 1, mmnt+1

		!do j = 1, neut

			!do k = 1, batch

				!print *, j + neut*(k-1)

				!PnMmntMatrix(k,:,i,j) = PnMmntData(j+neut*(k-1),:)*(j-1)**(i-1)

			!enddo

		!enddo

	!enddo

	!print *, PnMmntMatrix(101,:,1)

	call sleep(3600)

	open( unit = 120, file = 'ntrnmoment.data' )

	do i = 1, mmnt+1
				
		!do j = 1, neut

			do k = 1, batch*neut

				!write(120,*) PnMmntMatrix(k,:,i,j)
				!write(120,*) PnMmntMatrix(k,:,i)

			enddo

		!enddo

	enddo	

	close(unit = 120)

	call sleep(3600)

	open( unit = 5, file = 'PgMmnt0.txt' )

	do i = 1, batch

		do j = 1, gama

			do k = 1, N+1

				!do l = 1+100*(i-1), 100*i
				do l = 1+batch*(i-1), batch*i

					if ( PgData(l,k) .eq. (j-1) ) then

						Pg(j,k,i) = Pg(j,k,i) + 1

					endif

				enddo

			enddo

		enddo

		do p = 1, gama

			write(5,*) Pg(p,:,i)/real(batch)

		enddo 

	enddo

	close( unit = 5 )

	open( unit = 4, file = 'PgMmnt0.txt')

	do i = 1, batch*gama

		read( unit = 4, FMT = * ) PgMmntData(i,:)

	enddo

	close( unit = 4 )

	do i = 1, mmnt+1

		do j = 1, gama

			do k = 1, batch

				PgMmntMatrix(k,:,i,j) = PgMmntData(j+gama*(k-1),:)*(j-1)**(i-1)

			enddo

		enddo

	enddo

	open( unit = 121, file = 'gammamoment.data' )

	do i = 1, mmnt + 1
		
		do j = 1, gama

			do k = 1, batch

				write(121,*) PgMmntMatrix(k,:,i,j)

			enddo

		enddo

	enddo

	close( unit = 121 )

	!call sleep(3600)

	open( unit = 120, file = 'ntrnmoment.data' )

	!print *, neut*batch*(mmnt+1)
	!print *, batch

	do i = 1, neut*batch*(mmnt+1)

		read( unit = 120, FMT = * ) NeutData(i,:)

	enddo

	close( unit = 120 )

	do k = 1, mmnt + 1

		do i = 1, neut

			print *, 1+batch*(i-1)+(k-1)*batch*neut, batch*i+(k-1)*batch*neut

			do j = 1+batch*(i-1)+(k-1)*batch*neut, batch*i+(k-1)*batch*neut

				PnMeans(i,:,k) = PnMeans(i,:,k) + NeutData(j,:)

			enddo

			PnMeans(i,:,k) = PnMeans(i,:,k)/real(batch)

		enddo

	enddo

	print *, PnMeans(2,:,1)

	call sleep(3600)

	open( unit = 17, file = 'ntrnmeanmoment.data')

	do k = 1, mmnt+1

		do i = 1, neut

			write( 17, * ) PnMeans(i,:,k)

		enddo

	enddo

	close( unit = 17 )

	open( unit = 120, file = 'gammamoment.data')

	do i = 1, gama*batch*(mmnt+1)

		read( unit = 120, FMT = * ) GamaData(i,:)

	enddo

	close( unit = 120 )

	do k = 1, mmnt + 1

		do i = 1, gama

			do j = 1 + batch*(i-1)+(k-1)*batch*gama, batch*i+(k-1)*batch*gama

				PgMeans(i,:,k) = PgMeans(i,:,k) + GamaData(j,:)

			enddo

			PgMeans(i,:,k) = PgMeans(i,:,k)/real(batch)

		enddo
		
	enddo

	open( unit = 18, file = 'gammameanmoment.data')

	do k = 1, mmnt+1

		do i = 1, gama

			write(18,*) PgMeans(i,:,k)

		enddo

	enddo

	close( unit = 18 )


	print *, "I'm here!!!"
	!call sleep(3600)

	do k = 1, mmnt+1

		do i = 1, neut

			!print *, 1 + batch*(i-1)+(k-1)*batch*neut, batch*i+(k-1)*batch*neut

			do j = 1 + batch*(i-1)+(k-1)*batch*neut, batch*i+(k-1)*batch*neut

				varneut(i,:,k) = varneut(i,:,k) + (NeutData(j,:) - PnMeans(i,:,k))**2

			enddo

			varneut(i,:,k) = varneut(i,:,k)/(real(batch)-1)
			varneut(i,:,k) = 2*sqrt(varneut(i,:,k))/sqrt(real(batch))
			!varneut(i,:,k) = (1/sqrt(real(batch)))*sqrt(varneut(i,:,k))!/real(batch))

		enddo
	
	enddo

	!call sleep(3600)

	do k = 1, mmnt+1

		do i = 1, gama

			do j = 1 + batch*(i-1)+(k-1)*batch*gama, batch*i+(k-1)*batch*gama

				vargama(i,:,k) = vargama(i,:,k) + (GamaData(j,:) - PgMeans(i,:,k))**2

			enddo
			
			!vargama(i,:,k) = (1/sqrt(real(batch)))*sqrt(vargama(i,:,k))!/real(batch))
			vargama(i,:,k) = vargama(i,:,k)/(real(batch)-1)
			vargama(i,:,k) = 2*sqrt(vargama(i,:,k))/sqrt(real(batch))

		enddo

	enddo
	
	!print *, vargama(5,:,5)

	!open( unit = 23, file = 'final.txt')

	do k = 1, mmnt + 1

		do i = 1, neut

			do j = 1, N + 1

				PnMeanVar(i,1+2*(j-1),k) = PnMeans(i,j,k)
				PnMeanVar(i,2+2*(j-1),k) = varneut(i,j,k)

			enddo

		enddo

	enddo	

	do k = 1, mmnt + 1
	
		do i = 1, gama

			do j = 1, N + 1

				PgMeanVar(i,1+2*(j-1),k) = PgMeans(i,j,k)
				PgMeanVar(i,2+2*(j-1),k) = vargama(i,j,k)

			enddo	

		enddo

	enddo

	do k = 1, mmnt + 1

		open( unit = 50+k, file = filename(k) )
		!print *, filename(k)

		do i = 1, neut

			write( 50+k, * ) PnMeanVar(i,:,k)

		enddo

		close( unit = 50 + k )

	enddo	

	do k = 1, mmnt + 1

		open( unit = 60+k, file = filenamegamma(k) )
		!print *, filenamegamma(k)

		do i = 1, gama

			write( 60+k, * ) PgMeanVar(i,:,k)

		enddo

		close( unit = 60 + k )

	enddo	

	!print *, batch

END PROGRAM Png
