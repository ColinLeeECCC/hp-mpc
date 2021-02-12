  program cluster

    ! version 2.0
    ! Paul Makar, September 21, 2016
    ! Joana Soares, August 14, 2017
    ! Colin Lee,    January 4, 2021

    ! Find full hierarchical clustering information for gridded 
    ! NO2 data based on 1-R metric.

    use MPI
    use NETCDF
    use MHEAP
    
    implicit none

    ! selected data types
    integer, parameter          :: i8 = selected_int_kind(15)
    integer, parameter          :: sp = kind(1.0)
    integer, parameter          :: dp = selected_real_kind(2 * precision(1.0_sp))

    ! file unit for input arglist file
    integer                    :: fu_in

    ! input parameters
    character(256)             :: inDir, outDir
    integer                    :: start_year, start_mon, start_day
    integer                    ::   end_year,   end_mon,   end_day
    integer                    :: useFractionOfRegion
    character(8)               :: gemMachFieldName
    character(8)               :: fieldName
    integer                    :: strLen
    integer                    :: start_date, end_date, start_vec(8), end_vec(8)
    integer                    :: current_date, current_year, current_mon, current_day
    real                       :: start_julian, end_julian

    ! data file variables
    integer                    :: forecastHour
    character(256)             :: dataFileName
    integer                    :: ncId, variableId
    integer, dimension(NF90_MAX_VAR_DIMS) &
                               :: dimIds
    integer                    :: grid_ni, grid_nj, grid_nt
    integer                    :: gridi, gridj, gridi2, gridj2
    integer                    :: numPoints
    real*8, allocatable, dimension(:,:) &
                               :: buffer
    real*8, allocatable, dimension(:) &
                               :: tracer_sum, tracer_sqsum, &
                               tracer_xysum(:,:)
    ! integer, allocatable, dimension(:) &
    !                            :: tracer_n, tracer_xn
    real*8                     :: tracer_tmp, tracer_tmp2
    integer                    :: numDays, numTimesteps
    real*8                     :: SXY, SXX, SYY
    

    ! clustering information
    real*8,allocatable,dimension(:,:,:) &
                               :: NODES
    real*8,allocatable,dimension(:,:) &
                               :: NODESK1, NODESK2
    TYPE(THEAP), allocatable, dimension(:) &
                               :: PQueue
    integer                    :: K1, K2
    logical,allocatable,dimension(:) &
                               :: live
    integer,allocatable,dimension(:) &
                               :: clusterSize
    integer,allocatable,dimension(:,:) &
                               :: heapIdx
    real*8                     :: NODE(2), NODE1(2), NODE2(2)
    real*4                     :: mpiNode(3)
    real*4,allocatable,dimension(:,:) &
                               :: mpiNodes

    integer                    :: I, J, K, N, M
    integer                    :: IDAY, IHR, HRS_SINCE_START
    integer                    :: ierr

    real*8                     :: R, RMIN, RMAX, Nt

    integer*8                  :: countStart, countEnd, countRate, countMax
    integer                    :: uid
    logical                    :: saveDissMatrix, loadDissMatrix

    ! mpi variables
    integer, parameter         :: ROOT = 0
    integer                    :: myRank, numProcs
    integer                    :: status(MPI_STATUS_SIZE)
    real*8                     :: startTimer, endTimer
    integer                    :: nodeType

    ! mpi clustering variables
    integer                    :: II
    integer                    :: myPoints, numClustersThisNode
    integer,allocatable,dimension(:) &
                               :: myClusters,clusterRanks
    real*8                     :: KK1, KK2
    logical                    :: K1ThisRank, K2ThisRank
    integer                    :: K1Rank, K2Rank

    ! debugging
    logical,allocatable,dimension(:) :: visited 
    
    ! begin program
    
    ! initialize MPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, ierr)
    if (myRank == ROOT) WRITE(*,*) 'MPI Initizlied with ', numProcs, ' processors'

    ! Read in argument list from arglist.lst file
    fu_in = 15
    open(unit=fu_in, file='arglist.lst', status='old', iostat=ierr)
    IF (ierr .ne. 0) THEN
       WRITE(6,*) ' Couldn''t open "arglist.list"'
       call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
    ENDIF

    READ(fu_in, '(a256)') inDir ! data where compressed netcdf GEM-MACH data is stored
    READ(fu_in, '(i4,2(1x,i2.2))') start_year, start_mon, start_day
    READ(fu_in, '(i4,2(1x,i2.2))')   end_year,   end_mon,   end_day
    READ(fu_in, *) useFractionOfRegion
    READ(fu_in, '(a8)') gemMachFieldName
    READ(fu_in, '(a8)') fieldName
    READ(fu_in, '(a256)') outDir

    ! make sure inDir and outDir end in a '/'
    strLen = len_trim(inDir)
    IF (inDir(strLen:strLen) .ne. '/') THEN
       inDir = trim(inDir) // '/'
    ENDIF
    strLen = len_trim(outDir)
    IF (outDir(strLen:strLen) .ne. '/') THEN
       outDir = trim(outDir) // '/'
    ENDIF
    

    WRITE(*,*) ' Reading compressed netcdf data from ', trim(inDir)
    WRITE(*,101) start_year, start_mon, start_day, end_year, end_mon, end_day
101 FORMAT(' Doing ', i4, '-', i2, '-', i2, ' to ', i4, '-', i2, '-', i2)
    WRITE(*,102) useFractionOfRegion
102 FORMAT(' Using only 1/', i3, ' of available GEM-MACH domain')
    WRITE(*,*) ' Using field ', fieldName, '(', gemMachFieldName, ')'
    WRITE(*,*) ' Ouputting to ', trim(outDir)

    forecastHour = 18

    ! Read in a single file to get model domain size
103 FORMAT(A, i4.4, i2.2, i2.2, i2.2 '_', i6.6, 'p.netcdf4.compressed')
    WRITE(dataFileName, 103), trim(inDir), start_year, start_mon, start_day, &
         forecastHour, 60
    ierr = nf90_opeN( trim(dataFileName), NF90_NOWRITE, ncId )
    IF (ierr < 0) THEN
       WRITE(*,*) 'Error while opening ', trim(dataFileName)
       call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
    ENDIF
    
    ierr = nf90_inq_varid( ncId, trim(gemMachFieldName), variableId )
    ierr = nf90_inquire_variable( ncId, variableId, dimIds = dimIds )
    ierr = nf90_inquire_dimension( ncID, dimIds(1), len = grid_ni )
    ierr = nf90_inquire_dimension( ncID, dimIds(2), len = grid_nj )
    ierr = nf90_inquire_dimension( ncID, dimIds(3), len = grid_nt )

    ierr = nf90_close( ncId )

    ! grid_ni = grid_ni / useFractionOfRegion
    ! grid_nj = grid_nj / useFractionOfRegion
    grid_ni = useFractionOfRegion
    grid_nj = useFractionOfRegion
    numPoints = grid_ni * grid_nj
    
    IF ( myRank == ROOT) THEN
       write(*,*) 'dims', grid_ni, grid_nj, grid_nt
       write(*,*) 'numPoints = ', numPoints
    ENDIF

    ALLOCATE(myClusters(int(numPoints/numProcs) + 1))
    ALLOCATE(clusterRanks(numPoints))
    numClustersThisNode = int(numPoints / numProcs) + 1
    myPoints = -1
    
    !$OMP PARALLEL DO             &
    !$OMP DEFAULT ( SHARED )          &
    !$OMP PRIVATE ( I, II )           &
    !$OMP REDUCTION ( MAX: myPoints )
    DO II = 1, numClustersThisNode
       I = myRank + (II - 1) * numProcs + 1
       IF (I .LE. numPoints) THEN
          myClusters(II) = I
          myPoints = II
       ENDIF
    ENDDO
    !$OMP END PARALLEL DO
    numClustersThisNode = myPoints
    !$OMP PARALLEL DO             &
    !$OMP DEFAULT ( SHARED )          &
    !$OMP PRIVATE ( I, II )           &
    !$OMP REDUCTION ( MAX: myPoints )
    DO I = 1, numPoints
       clusterRanks(I) = MOD( I - 1, numProcs )
    ENDDO
    !$OMP END PARALLEL DO

    WRITE(*,*) ' Rank ', myRank, ' has ', myPoints, ' clusters to start.'
    WRITE(*,*) '      ', myClusters(1), ', ', myClusters(2), ', ', &
               myClusters(3), '...,', myClusters(numClustersThisNode)
    
    WRITE(dataFileName,107) trim(outDir), start_year, start_mon,     &
         start_day, end_year, end_mon, end_day, forecastHour,  &
         grid_ni, grid_nj
    ! disabling saving and of data for now
    saveDissMatrix = .false.
    uid = 15
    ! WRITE(*,*) ' Checking if ', trim(dataFileName), ' exists...'
    ! inquire(file=dataFileName, exist=loadDissMatrix)
    call SYSTEM_CLOCK( countStart, countRate, countMax )
    
    loadDissMatrix = .false.
    IF ( loadDissMatrix ) THEN
       WRITE(*,*) ' Yes. Loading dissimilarity matrix'
       saveDissMatrix = .false.

       WRITE(*,*) ' Allocating pQueue with len=', numPoints
       ALLOCATE(PQueue(numPoints))
       ALLOCATE(NODES(numPoints, numPoints,2))
       ALLOCATE(LIVE(numPoints))
       ALLOCATE(clusterSize(numPoints))
       ALLOCATE(heapIdx(numPoints,numPoints))
       WRITE(*,*) ' Allocated'

       OPEN(uid, file=trim(dataFileName), form='unformatted', status='old')
       DO N = 1, numPoints
          CALL PQueue(N)%INIT( int(numPoints,4), 2, GREATER1 )
          READ(uid) NODES(N,:,1)
          DO I = 1, numPoints
             IF (N .eq. I) CYCLE
             
             NODES(N,I,2) = real(I,8)

             CALL PQueue(N)%INSERT( NODES(N,I,:) )
             heapIdx(N,I) = PQueue(N)%SIZE()
          ENDDO
          LIVE(N) = .true.
          clusterSize(N) = 1
       ENDDO
       CLOSE(UID)
    ELSE
       WRITE(*,*) ' No. Computing dissimilarity matrix'

    ! read in all the time data from start to end
    start_vec = (/start_year,start_mon,start_day,forecastHour,00,00,00,00/)
    end_vec   = (/end_year,    end_mon,  end_day,forecastHour,00,00,00,00/)
    call d2j(start_vec, start_julian, ierr)
    call d2j(  end_vec,   end_julian, ierr)
    numDays = int(end_julian - start_julian + 1)
    numTimesteps = numDays * 24

108 FORMAT(A, i8.8, i2.2 '_', i6.6, 'p.netcdf4.compressed')
    
    ALLOCATE(BUFFER(GRID_NI, GRID_NJ))

    ALLOCATE(  TRACER_SUM( numPoints ))
    ALLOCATE(TRACER_SQSUM( numPoints ))
    ALLOCATE(TRACER_XYSUM( numClustersThisNode, numPoints ))

    TRACER_SUM = 0d0
    TRACER_SQSUM = 0d0
    TRACER_XYSUM = 0d0
    
    DO IDAY = 1, numDays
       call j2d(start_julian+IDAY*1.0, current_date, ierr)
       current_year = current_date / 10000
       current_mon  = current_date / 100 - current_year * 100
       current_day  = current_date - (current_year * 100 + current_mon) * 100
       WRITE(*,110) current_year, current_mon, current_day
110    FORMAT(' Loading data for day ', i4.4, '/', i2.2, '/', i2.2)

       ! each file contains 1 hourly dataset, so read in all 24 to get the full day
       DO IHR = 1, 24
          ! WRITE(*,'(10x,a,i2.2,a,i2.2)') 'Hr: ', IHR, ':', 0
          WRITE(dataFileName, 108), trim(inDir), current_date, &
               forecastHour, IHR * 60
          ierr = nf90_open( trim(dataFileName), NF90_NOWRITE, ncId, &
               comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
          IF (ierr < 0) THEN
             WRITE(*,*) 'Error while opening ', trim(dataFileName)
             call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
          ENDIF

          ierr = nf90_inq_varid( ncId, trim(gemMachFieldName), variableId )
          ierr = nf90_inquire_variable( ncId, variableId, dimIds = dimIds )
          ierr = nf90_get_var( ncid, variableId, buffer)

          ierr = nf90_close( ncId )

          HRS_SINCE_START = 24 * ( IDAY - 1 ) + IHR
          !$OMP PARALLEL DO                       &
          !$OMP DEFAULT( SHARED )                 &
          !$OMP PRIVATE( II, I, J, GRIDI, GRIDJ ) &
          !$OMP PRIVATE( TRACER_TMP, GRIDI2 )     &
          !$OMP PRIVATE( GRIDJ2, TRACER_TMP2 )
          DO II = 1,numClustersThisNode
             I = myClusters(II)
             ! GRIDI = (I - 1) / GRID_NJ + 1
             ! GRIDJ = (I - 1) - (GRIDI - 1) * GRID_NJ + 1
             GRIDJ = (I - 1) / GRID_NJ + 1
             GRIDI = (I - 1) - (GRIDJ - 1) * GRID_NJ + 1
             TRACER_TMP = BUFFER( GRIDI, GRIDJ )

             DO J = 1,numPoints
                ! GRIDI2 = (J - 1) / GRID_NJ + 1
                ! GRIDJ2 = (J - 1) - ( GRIDI2 - 1 ) * GRID_NJ + 1
                GRIDJ2 = (J - 1) / GRID_NJ + 1
                GRIDI2 = (J - 1) - ( GRIDJ2 - 1 ) * GRID_NJ + 1
                TRACER_TMP2 = BUFFER( GRIDI2, GRIDJ2 )
                
                TRACER_XYSUM( II, J ) = TRACER_XYSUM( II, J ) + &
                     TRACER_TMP * TRACER_TMP2
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          
          !$OMP PARALLEL DO                       &
          !$OMP DEFAULT( SHARED )                 &
          !$OMP PRIVATE( I, GRIDI, GRIDJ )        &
          !$OMP PRIVATE( TRACER_TMP )  
          DO I = 1,numPoints
             ! GRIDI = (I - 1) / GRID_NJ + 1
             ! GRIDJ = (I - 1) - (GRIDI - 1) * GRID_NJ + 1
             GRIDJ = (I - 1) / GRID_NJ + 1
             GRIDI = (I - 1) - (GRIDJ - 1) * GRID_NJ + 1
             TRACER_TMP = BUFFER( GRIDI, GRIDJ )

             TRACER_SUM( I )   = TRACER_SUM( I )   + TRACER_TMP
             TRACER_SQSUM( I ) = TRACER_SQSUM( I ) + TRACER_TMP * TRACER_TMP
          ENDDO
          !$OMP END PARALLEL DO
          
       ENDDO
    ENDDO

    call SYSTEM_CLOCK( countEnd, countRate, countMax )
    WRITE(*,*) 'Loading data and precalculating took ', &
         real(countEnd - countStart) / real(countRate), ' sec'
    countStart = countEnd

    RMIN = 9.999d9
    RMAX = -9.999d9
    Nt = real(numTimesteps, 8)
    !
    ! Okay,this is where things get a little bit mind-bending. We need to allocate
    ! a priority queue for each MPI process that will initially store 1 cluster
    ! for each pair. The cluster will get an ID, a score, and a cardinality. The
    ! priority queue sorts the clusters by score and returns the minimum quickly.
    ! As points are clustered, we need to store how many points are in a cluster
    ! so we can average the scores appropriately.
    WRITE(*,*) ' Allocating pQueue with len=', numPoints
    ALLOCATE(PQueue(numClustersThisNode))
    ALLOCATE(NODES(numClustersThisNode, numPoints,2))
    ALLOCATE(LIVE(numPoints))
    ALLOCATE(heapIdx(numClustersThisNode,numPoints))
    ALLOCATE(clusterSize(numPoints))
    WRITE(*,*) ' Allocated'

    clusterSize = 1
    live = .true.
    ! I have an inkling that doing this in the OMP loop is causing
    ! a runtime error
    DO N = 1,numClustersThisNode
       ! Create min-heap for this node. Will be zero if there are no edges
       CALL PQueue(N)%INIT( int(numPoints,4), 2, GREATER1 )
    ENDDO

    WRITE(*,*) "Initialized all pQueues"
    
    ! !$OMP PARALLEL DO               &
    ! !$OMP DEFAULT( SHARED )         &
    ! !$OMP PRIVATE( N, I, SXX, SYY ) &
    ! !$OMP PRIVATE( SXY, R, II )     &
    ! !$OMP REDUCTION(MAX: RMAX)      &
    ! !$OMP REDUCTION(MIN: RMIN)
    DO II = 1,numClustersThisNode
       ! WRITE(*,*) ' ', II
       N = myClusters(II)
       DO I = 1,numPoints
          IF ( I .eq. N ) CYCLE ! don't compute self-similarity
          ! WRITE(*,*) '   ', I
          ! compute the 1-R for this pair of clusters
          SXY = tracer_xysum(II, I) * Nt - &
               tracer_sum(N) * tracer_sum(I)
          SXX = tracer_sqsum(I) * Nt - &
               ( tracer_sum(I) * tracer_sum(I))
          SYY = tracer_sqsum(N) * Nt - &
               ( tracer_sum(N) * tracer_sum(N) )

          ! WRITE(*,*) '   ', SXX, ', ', SYY, ', ', SXY
          IF ( SXX .eq. 0d0 .or. SYY .eq. 0d0 ) THEN
             R = 0d0
             WRITE(*,*) 'Got bad statistics at ', N, ',', I
             WRITE(*,*) 'Ex = ', tracer_sum(N), 'Ey = ', tracer_sum(I)
             WRITE(*,*) 'Ex2 = ', tracer_sqsum(N), 'Ey2 = ', tracer_sqsum(I)
             WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
          ELSE
             R = SXY / SQRT( SXX * SYY )
             IF ( R < -1d0 .or. R > 1d0 ) THEN
                WRITE(*,*) 'Got weird R at ',II, '(', N, '),', I
                WRITE(*,*) ' R = ', R
                WRITE(*,*) 'Exy = ', tracer_xysum(II, I), ' Nt = ', Nt
                WRITE(*,*) 'Ex = ', tracer_sum(N), 'Ey = ', tracer_sum(I)
                WRITE(*,*) 'Ex2 = ', tracer_sqsum(N), 'Ey2 = ', tracer_sqsum(I)
                WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
             ENDIF
             IF (R < RMIN) RMIN = R
             IF (R > RMAX) RMAX = R
          ENDIF

          ! WRITE(*,*) ' Storing in Node ', II, I
          NODES(II,I,:) = (/ -1d0 * (1d0 - R), real(I,8) /)

112       FORMAT(' Inserting <', F6.3, ', ', I2,'> into queue ', I2)
          ! WRITE(*,112) NODES(II,I,1), INT(NODES(II,I,2)), II
          CALL PQueue(II)%INSERT( NODES(II,I,:), heapIdx(II,I) )
          ! WRITE(*,*) ' Done '
113       FORMAT(' At index(', i2, ',', i2, ') = ', i2)
          ! WRITE(*,113) II, I, heapIdx(II,I)
       ENDDO

    ENDDO
    ! !$OMP END PARALLEL DO
    call SYSTEM_CLOCK( countEnd, countRate, countMax )

    DEALLOCATE( TRACER_SUM )
    DEALLOCATE( TRACER_SQSUM )
    DEALLOCATE( TRACER_XYSUM )

 ENDIF ! .not. loadDissimilarityMatrix

    WRITE(*,*) 'R ranged from ', RMIN, ' to ', RMAX
    WRITE(*,*) 'Calculating initial dissimilarity scores took ', &
         real(countEnd - countStart) / real(countRate), ' sec'

    if (saveDissMatrix) THEN
       WRITE(*,*) 'Saving dissimilarity matrix...'
       WRITE(dataFileName,107) trim(outDir), start_year, start_mon,     &
            start_day, end_year, end_mon, end_day, forecastHour,  &
            grid_ni, grid_nj
107    FORMAT(a,i4.4,i2.2,i2.2, '_', i4.4,i2.2,i2.2,'_', i2.2, '_', &
            i4.4, 'x', i4.4, '.bin')
       open( uid, file=trim(dataFileName),form='unformatted')
       DO N = 1,numPoints
!          DO I = 1,numPoints
             WRITE(uid) NODES(N,:,1)
!          ENDDO
       ENDDO
       close(UID)
       WRITE(*,*) 'Done.'
    ENDIF
    
    ALLOCATE( mpiNodes( 3, numProcs ) )
    call MPI_TYPE_CONTIGUOUS( 3, MPI_REAL, nodeType, ierr )
    call MPI_TYPE_COMMIT( nodeType, ierr )
    ALLOCATE(NODESK1(2,numPoints))
    ALLOCATE(NODESK2(2,numPoints))
    ALLOCATE(VISITED(numPoints))

    call SYSTEM_CLOCK( countStart, countRate, countMax )

    ! Need to iterate this many times to get to a single cluster
    DO K = 1, numPoints-1
       
       ! find the queue with the min value
       RMIN = 9.999d9
       KK1 = -1
       DO II = 1,numClustersThisNode
          IF (.not. LIVE(myClusters(II))) CYCLE
          call PQueue(II).PEEK( 1, NODE )
          IF ( NODE(1) .le. RMIN ) THEN
             NODE1 = NODE
             RMIN = NODE1(1)
             KK1 = II
          ENDIF
       ENDDO
       IF (KK1 .eq. -1) THEN
          ! this mpi node has run out of clusters
          K1 = -1
          K2 = -1
       ELSE

          K1 = myClusters(KK1)
          ! now K1 is the cluster with the smallest dissimilarity
          K2 = NODE1(2)
          IF (K1 .lt. 1 .or. &
               K1 .eq. K2 .or.      &
               .not. LIVE(K1) .or. &
               .not. LIVE(K2) ) THEN
             WRITE(*,*) 'Something has gone wrong a. Print all queues'
             WRITE(*,*) 'K1 = ', K1, ', K2 = ', K2
             WRITE(*,*) 'LIVE(K1) = ',LIVE(K1), ', K2 = ', LIVE(K2)
             DO I = 1,numClustersThisNode
                IF (.not. LIVE(I)) CYCLE
                WRITE(*,*) 'Printing queue ', I
                DO WHILE( PQueue(I)%SIZE() .ne. 0)
                   call PQueue(I).POP( NODE )
                   WRITE(*,*), NODE
                ENDDO
             ENDDO
             call MPI_Abort(MPI_COMM_WORLD, 1, IERR)
          ENDIF
       ENDIF

       mpiNode = (/ real(K1,4), real(RMIN,4), real(K2,4) /)
       ! WRITE(*,*) ' Rank', myRank, ' sending ', mpiNode

       call MPI_ALLGATHER( mpiNode, 3, MPI_REAL, mpiNodes, 1, nodeType, &
            MPI_COMM_WORLD, ierr )

       K1 = -1
       rmin = 9.999d9
       DO I = 1,numProcs
117       FORMAT('R = ', f7.4, ', K1 = ', I1)
          ! WRITE(*,117) mpiNodes(2,I), INT(mpiNodes(1,I))
          IF ( mpiNodes(2,I) .le. rmin ) THEN
                K1   = mpiNodes(1,I)
                RMIN = mpiNodes(2,I)
                K2   = mpiNodes(3,I)
          ENDIF
       ENDDO
       ! To compare with old versions of the code,
       ! we should maintain the order, which always had
       ! K1 > K2. Since everything is symmetric, this should
       ! have no effect apart from that.
       IF (K1 .lt. K2) THEN
          J = K1
          K1 = K2
          K2 = J
       ENDIF
       IF (K1 .lt. 1 .or. &
            K1 .eq. K2 .or.      &
            .not. LIVE(K1) .or. &
            .not. LIVE(K2) ) THEN
          WRITE(*,*) 'Something has gone wrong b. Print all queues'
          WRITE(*,*) 'K1 = ', K1, ', K2 = ', K2
          WRITE(*,*) 'LIVE(K1) = ',LIVE(K1), ', K2 = ', LIVE(K2)
          DO I = 1,numClustersThisNode
             IF (.not. LIVE(I)) CYCLE
             WRITE(*,*) 'Printing queue ', I
             DO WHILE( PQueue(I)%SIZE() .ne. 0)
                call PQueue(I).POP( NODE )
                WRITE(*,*), NODE
             ENDDO
          ENDDO
          call MPI_Abort(MPI_COMM_WORLD, 1, IERR)
       ENDIF
       NODE1(1) = RMIN

       IF (myRank .eq. ROOT) THEN
          IF (numPoints .lt. 1000) THEN
             WRITE(*,106) k1, k2, NODE1(1)
          ELSEIF (numPoints .lt. 1000000) THEN
             WRITE(*,105) k1, k2, NODE1(1)
          ELSE
             WRITE(*,104) k1, k2, NODE1(1)
          ENDIF
       ENDIF
          
104    FORMAT('Clustering ', i9, ' and ', i9, ' (1-R) = ', F9.6)
105    FORMAT('Clustering ', i6, ' and ', i6, ' (1-R) = ', F9.6)
106    FORMAT('Clustering ', i3, ' and ', i3, ' (1-R) = ', F9.6)

       KK1 = (K1 - 1d0 - myRank)/numProcs + 1d0
       KK2 = (K2 - 1d0 - myRank)/numProcs + 1d0

       ! K1ThisRank = MOD( KK1, 1d0 ) .eq. 0d0
       K1Rank = clusterRanks(K1)
       K1ThisRank = K1Rank .eq. myRank

       K2Rank = clusterRanks(K2)
       K2ThisRank = K2Rank .eq. myRank

       ! Every node needs the NODES(K1,:,:) and NODES(:,K1,:)
       ! We can take advantage of the symmertry of the matrix
       ! and just send one of them
       IF ( K1ThisRank ) THEN
          ! WRITE(*,*) ' Rank ', myRank, ' filling NODESK1 to bcast'
          !$OMP PARALLEL DO         &
          !$OMP DEFAULT( SHARED )   &
          !$OMP PRIVATE( I )
          DO I = 1,numPoints
             NODESK1(:,I) = NODES(KK1, I, :)
          ENDDO
          !$OMP END PARALLEL DO
          ! WRITE(*,*) ' Rank ', myRank, ' done filling NODESK1'
       ENDIF
       ! WRITE(*,*), ' Rank', myRank, ' getting ready to bcast k1rank =', K1Rank
       call MPI_BCast( NODESK1, numPoints, MPI_2DOUBLE_PRECISION, &
            K1Rank, MPI_COMM_WORLD, IERR)
       ! WRITE(*,*), ' Rank', myRank, ' back from bcast. IERR = ', &
       !      IERR, ' NK1(1,:) = ', NODESK1(1,:), 'NK2(2,:) = ', NODESK1(2,:)
       
       IF ( K2ThisRank ) THEN
          ! WRITE(*,*) ' Rank ', myRank, ' filling NODESK2 to bcast'
          !$OMP PARALLEL DO         &
          !$OMP DEFAULT( SHARED )   &
          !$OMP PRIVATE( I )
          DO I = 1,numPoints
             NODESK2(:,I) = NODES(KK2, I, :)
          ENDDO
          !$OMP END PARALLEL DO
          ! WRITE(*,*) ' Rank ', myRank, ' done filling NODESK2'
       ENDIF
       ! WRITE(*,*), ' Rank', myRank, ' getting ready to bcast k2rank =', K2Rank
       call MPI_BCast( NODESK2, numPoints, MPI_2DOUBLE_PRECISION, &
            K2Rank, MPI_COMM_WORLD, IERR)
       ! WRITE(*,*), ' Rank', myRank, ' back from bcast. IERR = ', &
       !      IERR, ' NK2(1,:) = ', NODESK2(1,:), 'NK2(2,:) = ', NODESK2(2,:)

       
       ! WRITE(*,*)'rank ',  myRank, ', deleting K2 = ', K2
       LIVE(K2) = .false.

       IF ( K1ThisRank ) THEN
          ! WRITE(*,*) 'K1 is on rank ', myRank, ' at index ', KK1, '. Clearing.'
          call PQueue(KK1)%CLEAR()
       ENDIF

       ! debugging information: check we visisted every cluster
       visited(:) = .false.
       
       !$OMP PARALLEL DO                    &
       !$OMP DEFAULT( SHARED )              &
       !$OMP PRIVATE( I, II, NODE1, NODE2)  &
       !$OMP PRIVATE( M, R, NODE )          
       DO I = 1,numPoints
          
          ! WRITE(*,'(a, i1)') ' I = ', I
          
114       FORMAT('Cluster ', i2, ' LIVE? ', L1, ' .ne. K1? ', L1)
          ! WRITE(*,114) I, LIVE(I), I.ne.K1
          ! IF (.not. LIVE(I)) CYCLE
          ! IF (I .eq. K1) CYCLE
          IF ( LIVE(I) .and. I .ne. K1 .and. &
               clusterRanks(I) .eq. myRank) THEN
             VISITED(I) = .true.
             
             II = (I - 1d0 - myRank)/numProcs + 1d0
             IF (myClusters(II) .ne. I) THEN
                WRITE(*,*) ' CR(I) = ', clusterRanks(I), ' myRank = ', myRank
                WRITE(*,*) ' I =', I, ' II = ', II, ' mC(II) = ', myClusters(II)
                call MPI_Abort(MPI_COMM_WORLD, 3, IERR)
             ENDIF

             NODE1(:) = NODES(II,K1,:)
             M = PQueue(II)%SIZE()
111          FORMAT('Cluster ', i2, ' deleting node ', 2f6.3, ' at index (', i2, ', ', i2, ') = ', i2)
             ! WRITE(*,111) I, NODE1, II, K1, PQueue(II)%INDEXAT(heapIdx(II,K1))
             call PQueue(II)%DELETE( K=heapIdx(II,K1), DNODE=NODE )
             IF (.not. ALL(NODE1 .eq. NODE) ) THEN
120             FORMAT('Step ', i5, ', k1 = ', i5, ', k2 = ', i5)
121             FORMAT('N1 = ', f7.4, ', ', i5)
122             FORMAT('N  = ', f7.4, ', ', i5)
                WRITE(*,*) ' Deleted node not equal to the one we wanted.'
                WRITE(*,120) K, K1, K2
                WRITE(*,121) NODE1(1), INT(NODE1(2))
                WRITE(*,122) NODE(1), INT(NODE(2))
             ENDIF
             IF (M - PQueue(II)%SIZE() .ne. 1) THEN
                WRITE(*,*) ' Failed to delete node1 ', NODE1, ' from queue ', I
                WRITE(*,*) ' M = ', M, ' size = ', PQueue(II)%SIZE()
             ENDIF
             NODE2(:) = NODES(II,K2,:)
             ! WRITE(*,111) I, NODE2, II, K2, PQueue(II)%INDEXAT(heapIdx(II,K2))
             M = PQueue(II)%SIZE()
             call PQueue(II)%DELETE( K=heapIdx(II,K2), DNODE=NODE )
             IF (.not. ALL(NODE2 .eq. NODE) ) THEN
                WRITE(*,*) ' Deleted node not equal to the one we wanted.'
                WRITE(*,120) K, K1, K2
                WRITE(*,121) NODE2(1), INT(NODE2(2))
                WRITE(*,122) NODE(1), INT(NODE(2))
             ENDIF
             IF (M - PQueue(II)%SIZE() .ne. 1)  THEN
                WRITE(*,*) ' Failed to delete node2 ', NODE2, ' from queue ', I
                WRITE(*,*) ' M = ', M, ' size = ', PQueue(II)%SIZE()
             ENDIF

             ! R = (clusterSize(K1) * NODE1(1) + clusterSize(K2) * NODE2(1)) &
             !      / ( clusterSize(K1) + clusterSize(K2) )
             R = MIN( NODE1(1), NODE2(1) )

             NODES(II,K1,1) = R
116          FORMAT('PQueue(',i2,')%INSERT( [', f7.4, ', ', i2, '])')
             ! WRITE(*,116) I, R, K1
             call PQueue(II)%INSERT( NODES(II,K1,:), heapIdx(II,K1) )

             IF ( K1ThisRank ) THEN
                ! WRITE(*,116) K1, R, I
                NODES(KK1, I,1) = R
                !$OMP CRITICAL
                call PQueue(KK1)%INSERT( NODES(KK1, I,:), heapIdx(KK1,I) )
                !$OMP END CRITICAL
             ENDIF
          ELSEIF ( LIVE(I) .and. I .ne. K1 ) THEN

             VISITED(I) = .true.

             NODE1(:) = NODESK1(:,I)
             NODE2(:) = NODESK2(:,I)
             
             R = MIN( NODE1(1), NODE2(1) )

             IF ( K1ThisRank ) THEN
                ! WRITE(*,116) K1, R, I+J
                NODES(KK1, I,1) = R
                !$OMP CRITICAL
                call PQueue(KK1)%INSERT( NODES(KK1, I,:), heapIdx(KK1,I) )
                !$OMP END CRITICAL
             ENDIF
          ENDIF
       ENDDO
       !$OMP END PARALLEL DO

       ! check that we hit each cluster
123    FORMAT('WARNING: Step', i5, ': Did not visit cluster ', i5, ', K1 = ', i5, ', K2 = ', i5)
       DO I = 1,numPoints
          IF (.not. LIVE(I)) CYCLE
          IF (I .eq. K1) CYCLE
          IF (.not. VISITED(I)) THEN
             WRITE(*,123) K, I, K1, K2
          ENDIF
       ENDDO
       
       
       ! Check the length of each queue. THey should all be numPoints - K - 1.
       N = (numPoints - K - 1)
       DO II = 1,numClustersThisNode
          IF (.not. LIVE(myClusters(II))) CYCLE
          M = PQueue(II)%SIZE()
          IF ( M .ne. N ) THEN
118          FORMAT('WARNING: Step', i5,': Cluster ', i5, '(', i4, ') has size ', i5, ' not ', i5)
             WRITE(*,118) k, myClusters(II), II, M, N
          ENDIF
       ENDDO
       
       clusterSize(K1) = clusterSize(K1) + clusterSize(K2)
       clusterSize(K2) = 0
    ENDDO
    call SYSTEM_CLOCK( countEnd, countRate, countMax )
    WRITE(*,*) 'Clustering took ', &
         real(countEnd - countStart) / real(countRate), ' sec'
    countStart = countEnd

    call MPI_FINALIZE(IERR)

  CONTAINS

    subroutine d2j(dat,julian,ierr)
      !---------------------------------------------------------------
      ! * Author:    John S. Urban
      ! * Version:   1.0 2015-12-21
      ! * Reference: From Wikipedia, the free encyclopedia 2015-12-19
      ! * There is no year zero
      ! * Julian Day must be non-negative
      ! * Julian Day starts at noon; while Civil Calendar date starts at midnight
      !---------------------------------------------------------------
      !      character(len=:),parameter :: ident=,  "Converts proleptic Gregorian date array to Julian Day
      integer,intent(in)         :: dat(8)   ! array like returned by DATE_AND_TIME(3f)

      real(kind=4),intent(out)  :: julian   ! Julian Day (non-negative, but may be non-integer)
      integer,intent(out)        :: ierr    ! Error return, 0 for successful execution,-1=invalid year,-2=invalid month,-3=invalid day,
      ! -4=invalid date (29th Feb, non leap-year)
      integer                 :: year, month, day, utc, hour, minute
      real(kind=4)           :: second
      integer                 :: A, Y, M, JDN

      year   = dat(1)                        ! Year
      month  = dat(2)                        ! Month
      day    = dat(3)                        ! Day
      utc    = dat(4)*60                     ! Delta from UTC, convert from minutes to seconds
      hour   = dat(5)                        ! Hour
      minute = dat(6)                        ! Minute
      second = dat(7)-utc+dat(8)/1000.0d0    ! Second   ! correction for time zone and milliseconds
      !
      julian = -HUGE(99999)                  ! this is the date if an error occurs and IERR is < 0
      !
      if(year==0 .or. year .lt. -4713) then
         ierr=-1
         return
      endif
      !
      !  You must compute first the number of years (Y) and months (M) since March 1st -4800 (March 1, 4801 BC)
      A=(14-month)/12 ! A will be 1 for January or Febuary, and 0 for other months, with integer truncation
      Y=year+4800-A
      M=month+12*A-3  ! M will be 0 for March and 11 for Febuary
      !  All years in the BC era must be converted to astronomical years, so that 1BC is year 0, 2 BC is year "-1", etc.
      !  Convert to a negative number, then increment towards zero
      !  Staring from a Gregorian calendar date
      JDN=day + (153*M+2)/5 + 365*Y + Y/4 - Y/100 + Y/400 - 32045  !  with integer truncation
      !  Finding the Julian date given the JDN (Julian day number) and time of day
      julian=JDN + dble(hour-12)/24.0d0 + dble(minute)/1440.0d0 + second/86400.0d0
      !
      if(julian.lt.0.d0) then                  ! Julian Day must be non-negative
         ierr=1
      else
         ierr=0
      endif

      return

    end subroutine d2j
    !
    !************************************************************************************
    !
    subroutine j2d(julian,dat,ierr)
      !---------------------------------------------------------------
      ! * Author:    John S. Urban
      ! * Version:   1.0 2015-12-21
      ! * Reference: From Wikipedia, the free encyclopedia 2015-12-19
      !---------------------------------------------------------------

      real(kind=4),intent(in)  :: julian   ! Julian Day (non-negative, but may be non-integer)
      integer,intent(out)      :: dat      ! array like returned by DATE_AND_TIME(3f)
      integer,intent(out)      :: ierr     ! Error return, 0 for successful execution,-1=invalid year,-2=invalid month,-3=invalid day,
      ! -4=invalid date (29th Feb, non leap-year)
      integer             :: year, month, day, hour, minute
      real(kind=4)        :: second
      integer             :: timezone(8), tz, jalpha,ja,jb,jc,jd,je,ijul
      integer,dimension(8) :: datAr
      real, parameter :: secday = 86400.0

      if(julian.lt.0.d0) then                  ! Julian Day must be non-negative
         ierr=1
      else
         ierr=0
      endif

      ijul=int(julian)                           ! Integral Julian Day
      !---------------------------------------------
      jalpha=idint((dble(ijul-1867216)-0.25d0)/36524.25d0) ! Correction for Gregorian Calendar
      ja=ijul+1+jalpha-idint(0.25d0*dble(jalpha))
      !---------------------------------------------

      jb=ja+1524
      jc=idint(6680.d0+(dble(jb-2439870)-122.1d0)/365.25d0)
      jd=365*jc+idint(0.25d0*dble(jc))
      je=idint(dble(jb-jd)/30.6001d0)
      day=jb-jd-idint(30.6001d0*dble(je))
      month=je-1

      if(month.gt.12)then
         month=month-12
      endif

      year=jc-4715
      if(month.gt.2)then
         year=year-1
      endif

      if(year.le.0)then
         year=year-1
      endif

      datAr(1)=year
      datAr(2)=month
      datAr(3)=day
      dat = datAr(1)*10000+datAr(2)*100+datAr(3)
      ierr=0

      return
    end subroutine j2d


    LOGICAL FUNCTION GREATER1( NODE1, NODE2 )
      DOUBLE  PRECISION, INTENT(IN) :: NODE1(:), NODE2(:)
      GREATER1 = NODE1(1) < NODE2(1)
    END FUNCTION GREATER1

    LOGICAL FUNCTION GREATER2( NODE1, NODE2 )
      DOUBLE PRECISION, INTENT(IN) :: NODE1(:), NODE2(:)
      GREATER2 = NODE1(2) < NODE2(2)
    END FUNCTION GREATER2

  END PROGRAM CLUSTER
