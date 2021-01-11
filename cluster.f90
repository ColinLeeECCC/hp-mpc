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
    real*8, allocatable, dimension(:,:) &
                               :: buffer, tracer_time_ser
    real*8, allocatable, dimension(:) &
                               :: tracer_sum, tracer_sqsum, &
                               tracer_xysum
    ! integer, allocatable, dimension(:) &
    !                            :: tracer_n, tracer_xn
    real*8                     :: tracer_tmp, tracer_tmp2
    integer                    :: numDays, numTimesteps
    real*8                     :: SXY, SXX, SYY
    

    ! clustering information
    integer(kind=8)            :: numPoints, numPairs
    real*8                     :: minNode(2), globalMinNode(2)
    integer                    :: minRank

    ! related to individual processes
    integer(kind=8)            :: myPoints,     myPointMin, myPointMax
    integer(kind=8)            :: myPairMin,    myPairMax
    integer(kind=8)            :: myPairsStart, myPairsEnd, myPairs
    integer,allocatable,dimension(:) &
                               :: myPairsI, myPairsJ
    integer,allocatable,dimension(:) &
                               :: numEdges ! store the number of Js for each I since
                                           ! many will be zero we can avoid allocating those
    
    ! mpi information
    integer, parameter         :: ROOT = 0
    integer                    :: myRank, numProcs
    integer                    :: ierr
    integer                    :: status(MPI_STATUS_SIZE)
    real*8                     :: startTimer, endTimer

    ! priority queue variables
    real*8                     :: node(3)
    TYPE(THEAP), allocatable, dimension(:) &
                               :: PQueue
    
    ! miscellanceous variables
    integer(kind=8)            :: I, J, II, JJ
    integer(kind=8)            :: L, M, N
    integer(kind=8)            :: IND, IND1
    integer                    :: IDAY, IHR, HRS_SINCE_START
    real*8                     :: RMIN, RMAX, R, Nt

    
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
       STOP
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
    

    IF (myRank == ROOT) THEN
       WRITE(*,*) ' Reading compressed netcdf data from ', trim(inDir)
       WRITE(*,101) start_year, start_mon, start_day, end_year, end_mon, end_day
101    FORMAT(' Doing ', i4, '-', i2, '-', i2, ' to ', i4, '-', i2, '-', i2)
       WRITE(*,102) useFractionOfRegion
102    FORMAT(' Using only 1/', i3, ' of available GEM-MACH domain')
       WRITE(*,*) ' Using field ', fieldName, '(', gemMachFieldName, ')'
       WRITE(*,*) ' Ouputting to ', trim(outDir)
    ENDIF

    forecastHour = 18

    ! Read in a single file to get model domain size
103 FORMAT(A, i4.4, i2.2, i2.2, i2.2 '_', i6.6, 'p.netcdf4.compressed')
    WRITE(dataFileName, 103), trim(inDir), start_year, start_mon, start_day, &
         forecastHour, 60
    ierr = nf90_opeN( trim(dataFileName), NF90_NOWRITE, ncId, &
         comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
    IF (ierr < 0) THEN
       WRITE(*,*) 'Error while opening ', trim(dataFileName)
       STOP
    ENDIF
    
    ierr = nf90_inq_varid( ncId, trim(gemMachFieldName), variableId )
    ierr = nf90_inquire_variable( ncId, variableId, dimIds = dimIds )
    ierr = nf90_inquire_dimension( ncID, dimIds(1), len = grid_ni )
    ierr = nf90_inquire_dimension( ncID, dimIds(2), len = grid_nj )
    ierr = nf90_inquire_dimension( ncID, dimIds(3), len = grid_nt )

    ierr = nf90_close( ncId )

    grid_ni = grid_ni / useFractionOfRegion
    grid_nj = grid_nj / useFractionOfRegion

    IF ( myRank == ROOT) THEN
       write(*,*) 'dims', grid_ni, grid_nj, grid_nt
    ENDIF

    startTimer = MPI_Wtime()
    ! Initiallie we compute the pairwise dissimilarity between each station
    numPoints = grid_ni * grid_nj
    numPairs  = numPoints * (numPoints - 1) / 2

    IF ( myRank == ROOT) THEN
       write(*,*) 'dims', grid_ni, grid_nj, grid_nt
       write(*,*) 'numPoints = ', numPoints, ' numPairs = ', numPairs
    ENDIF
    
    myPairsStart = myRank * numPairs / (numProcs*1_i8) + 1
    myPairsEnd   = (myRank + 1) * numPairs / (numProcs*1_i8)
    myPairsEnd   = min(myPairsEnd, numPairs)
    myPairs      = myPairsEnd - myPairsStart + 1

    ! pair indices in point space (not grid space)
    ALLOCATE(myPairsI(myPairs))
    ALLOCATE(myPairsJ(myPairs))

    !
    ! save the subscripts of the upper-triangular dissimilarity matrix
    !
    !$OMP PARALLEL DO                &
    !$OMP DEFAULT( SHARED )          &
    !$OMP PRIVATE( I, J, IND, IND1 )
    DO I = 1, numPoints-1
       IND1 = (I - 1_i8) * numPoints - (I * (I + 1_i8)) / 2_i8
       DO J = I+1,numPoints
          IND = J + IND1
          IF (IND >= myPairsStart .and. &
               IND <= myPairsEnd ) THEN
             myPairsI(IND - myPairsStart + 1) = I
             myPairsJ(IND - myPairsStart + 1) = J
          ELSEIF (IND > myPairsEnd) THEN
             EXIT
          ENDIF
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    myPointMin = MIN( MINVAL( myPairsI ), MINVAL( myPairsJ ) )
    myPointMax = MAX( MAXVAL( myPairsI ), MAXVAL( myPairsJ ) )
    myPoints   = myPointMax - myPointMin + 1
!    if (myRank == ROOT) &
         WRITE(*,*) ' I have ', myPoints, ' points to care about. (', myPointMin,&
                    ' to ', myPointMax

    endTimer = MPI_Wtime()
    
    IF (numProcs > 1) THEN
       IF ( myRank == ROOT ) THEN
          call MPI_SEND( I, 1, MPI_INTEGER, myRank + 1, 0, MPI_COMM_WORLD, ierr)
       ELSE
          call MPI_PROBE( MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          call MPI_RECV( I, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, ierr)
          IF ( myRank < numProcs - 1) THEN
             call MPI_SEND( I, 1, MPI_INTEGER, myRank + 1, 0, MPI_COMM_WORLD, ierr)
          ENDIF
       ENDIF
    ENDIF

    ! in order to not have a million blank spaces but also allow these to print
    ! information over the wide range of possible pairs, use formats conditionally
    IF (numPairs < 1000000) THEN
       WRITE(*,104) myRank, myPairsStart, myPairsEnd, myPairsI(1), myPairsJ(1), &
            myPairsI(size(myPairsI)), myPairsJ(size(myPairsJ))
    ELSEIF (numPairs < 1000000000) THEN
       WRITE(*,106) myRank, myPairsStart, myPairsEnd, myPairsI(1), myPairsJ(1), &
            myPairsI(size(myPairsI)), myPairsJ(size(myPairsJ))
    ELSE
       WRITE(*,107) myRank, myPairsStart, myPairsEnd, myPairsI(1), myPairsJ(1), &
            myPairsI(size(myPairsI)), myPairsJ(size(myPairsJ))
    ENDIF

    WRITE(*,105) myRank, endTimer - startTimer
    startTimer = endTimer
    
104 FORMAT(' Rank ', i3, ' has pairs ', i6, ' to ', i6, ' going from I,J = ', i3, &
         ', ', i3, ' to I,J = ', i3, ', ', i3)
105 FORMAT(' Rank ', i3, ' took ', f12.4, ' seconds.')
106 FORMAT(' Rank ', i3, ' has pairs ', i9, ' to ', i9, ' going from I,J = ', i6, &
         ', ', i6, ' to I,J = ', i6, ', ', i6)
107 FORMAT(' Rank ', i3, ' has pairs ', i14, ' to ', i14, ' going from I,J = ', i9, &
         ', ', i9, ' to I,J = ', i9, ', ', i9)

    ! read in all the time data from start to end
    start_vec = (/start_year,start_mon,start_day,forecastHour,00,00,00,00/)
    end_vec   = (/end_year,    end_mon,  end_day,forecastHour,00,00,00,00/)
    call d2j(start_vec, start_julian, ierr)
    call d2j(  end_vec,   end_julian, ierr)
    numDays = int(end_julian - start_julian + 1)
    numTimesteps = numDays * 24

108 FORMAT(A, i8.8, i2.2 '_', i6.6, 'p.netcdf4.compressed')
! 109 FORMAT(' Allocating ', i5, ', ', i5, ' for ', A)
!     IF (myRank == ROOT) THEN
!        WRITE(*,109) GRID_NI, GRID_NJ, 'BUFFER'
!        WRITE(*,109) myPairs, numDays*24, 'TRACER_TIME_SER'
!     ENDIF
       
    
    ALLOCATE(BUFFER(GRID_NI, GRID_NJ))
!     IF (myRank == ROOT) &
!          WRITE(*,*) ' Allocated BUFFER'

    ALLOCATE(TRACER_TIME_SER( myPoints, numDays * 24 ))
!     IF (myRank == ROOT) &
!          WRITE(*,*) ' Allocated TRACER_TIME_SER'
    ALLOCATE(  TRACER_SUM( myPoints ))
    ! ALLOCATE(    TRACER_N( myPoints ))
    ALLOCATE(TRACER_SQSUM( myPoints ))
    ALLOCATE(    numEdges( myPoints ))
    ALLOCATE(TRACER_XYSUM( myPairs ))
    ! ALLOCATE(   TRACER_XN( myPairs ))

    TRACER_SUM = 0d0
    TRACER_SQSUM = 0d0
    TRACER_XYSUM = 0d0
    numEdges = 0
    ! TRACER_N = 0
    ! TRACER_XN = 0
    
    DO IDAY = 1, numDays
       call j2d(start_julian+IDAY*1.0, current_date, ierr)
       IF (myRank == ROOT) THEN
          current_year = current_date / 10000
          current_mon  = current_date / 100 - current_year * 100
          current_day  = current_date - (current_year * 100 + current_mon) * 100
          WRITE(*,110) current_year, current_mon, current_day
110       FORMAT(' Loading data for day ', i4.4, '/', i2.2, '/', i2.2)
       ENDIF
       ! each file contains 1 hourly dataset, so read in all 24 to get the full day
       DO IHR = 1, 24
           IF (myRank == ROOT) &
                WRITE(*,'(10x,a,i2.2,a,i2.2)') 'Hr: ', IHR, ':', 0
          WRITE(dataFileName, 108), trim(inDir), current_date, &
               forecastHour, IHR * 60
!           IF (myRank == ROOT) &
!                WRITE(*,*) '   Opening "', trim(dataFileName), '"'
          ierr = nf90_opeN( trim(dataFileName), NF90_NOWRITE, ncId, &
               comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
!           IF (myRank == ROOT) &
!                WRITE(*,*) '   Opened!'
          IF (ierr < 0) THEN
             WRITE(*,*) 'Error while opening ', trim(dataFileName)
             STOP
          ENDIF

!           IF (myRank == ROOT) &
!                WRITE(*,*) '   Inquring variable "', trim(gemMachFieldName), '" id'
          ierr = nf90_inq_varid( ncId, trim(gemMachFieldName), variableId )
!           IF (myRank == ROOT) &
!                WRITE(*,*) '   Inquring variable ', variableId
          ierr = nf90_inquire_variable( ncId, variableId, dimIds = dimIds )
!           IF (myRank == ROOT) &
!                WRITE(*,*) '   Getting data.'
          ierr = nf90_get_var( ncid, variableId, buffer)

          ierr = nf90_close( ncId )

          HRS_SINCE_START = 24 * ( IDAY - 1 ) + IHR
!           IF (myRank == ROOT) &
!               WRITE(*,*) '   Hrs_since_start = ', HRS_SINCE_START
          !$OMP PARALLEL DO                      &
          !$OMP DEFAULT( SHARED )                &
          !$OMP PRIVATE( I, J, GRIDI, GRIDJ, L ) &
          !$OMP PRIVATE( TRACER_TMP, GRIDI2 )    &
          !$OMP PRIVATE( GRIDJ2, TRACER_TMP2 )
          DO II = 1,myPoints
             I = II + myPointMin - 1
             GRIDI = (I - 1) / GRID_NJ + 1
             GRIDJ = (I - 1) - (GRIDI - 1) * GRID_NJ + 1
             TRACER_TMP = BUFFER( GRIDI, GRIDJ )
!             IF (myRank == ROOT .and. IHR == 1) &
!                  WRITE(*,*) 'TR(', GRIDI, ',', GRIDJ,') = ', TRACER_TMP
             TRACER_TIME_SER( II, HRS_SINCE_START ) = TRACER_TMP
             TRACER_SUM( II )   = TRACER_SUM( II )   + TRACER_TMP
             TRACER_SQSUM( II ) = TRACER_SQSUM( II ) + TRACER_TMP * TRACER_TMP
             ! TRACER_N(II) = TRACER_N(II) + 1
             DO L = 1,myPairs
                IF ( myPairsI(L) .EQ. I ) THEN
                   numEdges(II) = numEdges(II) + 1
                   J = myPairsJ(L)
                   GRIDI2 = (J - 1) / GRID_NJ + 1
                   GRIDJ2 = (J - 1) - ( GRIDI2 - 1 ) * GRID_NJ + 1
                   TRACER_TMP2 = BUFFER( GRIDI2, GRIDJ2 )
                
                   TRACER_XYSUM( L ) = TRACER_XYSUM( L ) + &
                        TRACER_TMP * TRACER_TMP2
                   ! TRACER_XN(L) = TRACER_XN(L) + 1
                ENDIF
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          
       ENDDO
    ENDDO

    endTimer = MPI_Wtime()
    WRITE(*,105) myRank, endTimer - startTimer
    startTimer = endTimer

    ! J = 0
    ! DO I = 1,myPoints
    !    IF (TRACER_N(I) .NE. numTimesteps .and. &
    !         TRACER_N(I) .NE. 0) &
    !         WRITE(*,*) 'N(', I, ') = ', TRACER_N(I)
    !    IF (TRACER_N(I) .EQ. numTimesteps) &
    !         J = J + 1
    ! ENDDO
    ! WRITE(*,*) J, ' points had valid sums'
    ! J = 0
    ! DO L = 1,myPairs
    !    IF (TRACER_XN(L) .NE. numTimesteps .and. &
    !         TRACER_XN(L) .NE. 0) &
    !         WRITE(*,*) 'XN(', L, ') = ', TRACER_XN(L)
    !    IF (TRACER_XN(L) .EQ. numTimesteps) &
    !         J = J + 1
    ! ENDDO
    ! WRITE(*,*) J, ' points had valid sums'
    
    IF (myRank == ROOT) THEN
         WRITE(*,*) 'TRACER_SUM(1:10) = ',   TRACER_SUM(1:10)
         WRITE(*,*) 'TRACER_SQSUM(1:10) = ', TRACER_SQSUM(1:10)
      ENDIF

    DEALLOCATE( BUFFER )

    RMIN = 999d99
    RMAX = -999d99
    Nt = real(numTimesteps, 8)
    !
    ! Okay,this is where things get a little bit mind-bending. We need to allocate
    ! a priority queue for each MPI process that will initially store 1 cluster
    ! for each pair. The cluster will get an ID, a score, and a cardinality. The
    ! priority queue sorts the clusters by score and returns the minimum quickly.
    ! As points are clustered, we need to store how many points are in a cluster
    ! so we can average the scores appropriately.
    if (myRank == ROOT) &
         WRITE(*,*) ' Allocating pQueue with len=', myPoints
    ALLOCATE(PQueue(myPoints))
    if (myRank == ROOT) &
         WRITE(*,*) ' Allocated'
    ! This loop may not be OMP-able as it stands
    DO II = 1,myPoints
       I = II + myPointMin - 1
       ! Create min-heap for this node. Will be zero if there are no edges
       CALL PQueue(II)%INIT( numEdges(II), 3, GREATER1 )
       
       IF (numEdges(II) .GT. 0) THEN
          DO L = 1,myPairs
             IF (I .EQ. myPairsI(L)) THEN
                J = myPairsJ(L)
                JJ = J - myPointMin + 1
                ! Store the ID of the cluster we're pairing with
                NODE(2) = real(J, 8)
                ! Cardinality is 1 since there's only one node per cluster
                NODE(3) = 1d0
                ! compute the 1-R for this pair of clusters
                ! I'm having trouble understanding the math in the original code.
                ! It doesn't seem to actually compute pearson R. Optimised
                ! to avoid divisions
                SXY = tracer_xysum(L) * Nt - &
                     tracer_sum(II) * tracer_sum(JJ)
                SXX = tracer_sqsum(II) * Nt - &
                     ( tracer_sum(II) * tracer_sum(II))
                SYY = tracer_sqsum(JJ) * Nt - &
                     ( tracer_sum(JJ) * tracer_sum(JJ) )

                ! IF (myRank == ROOT) THEN
                !    WRITE(*,*) 'numTimesteps = ', numTimesteps, ' Exy = ',tracer_xysum(L)
                !    WRITE(*,*) 'Ex = ', tracer_sum(I), 'Ey = ', tracer_sum(J)
                !    WRITE(*,*) 'Ex2 = ', tracer_sqsum(I), 'Ey2 = ', tracer_sqsum(J)
                !    WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
                !    WRITE(*,*) 'R = ', SXY / SQRT( SXY * SYY )
                ! ENDIF
                IF ( SXX .eq. 0d0 .or. SYY .eq. 0d0 ) THEN
                   NODE(1) = 1d0
                   WRITE(*,*) 'Rank ', myRank, ' got bad statistics at ', I, ',', J
                   WRITE(*,*) 'Ex = ', tracer_sum(I), 'Ey = ', tracer_sum(J)
                   WRITE(*,*) 'Ex2 = ', tracer_sqsum(I), 'Ey2 = ', tracer_sqsum(J)
                   WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
                ELSE
                   R = SXY / SQRT( SXX * SYY )
                   IF ( R < -1d0 .or. R > 1d0 ) THEN
                      WRITE(*,*) 'Rank ', myRank, ' got weird R at ', I, ',', J
                      WRITE(*,*) ' R = ', R
                      WRITE(*,*) 'Exy = ', tracer_xysum(L), ' Nt = ', Nt
                      WRITE(*,*) 'Ex = ', tracer_sum(I), 'Ey = ', tracer_sum(J)
                      WRITE(*,*) 'Ex2 = ', tracer_sqsum(I), 'Ey2 = ', tracer_sqsum(J)
                      WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
                      WRITE(*,*) 'Tx(', MINVAL(TRACER_TIME_SER( N, : ) ), ', ', &
                           MAXVAL( TRACER_TIME_SER( N, : ) ), ')'
                      WRITE(*,*) 'Ty(', MINVAL(TRACER_TIME_SER( J, : ) ), ', ', &
                           MAXVAL( TRACER_TIME_SER( J, : ) ), ')'
                   ENDIF
                   IF (R < RMIN) RMIN = R
                   IF (R > RMAX) RMAX = R
                   NODE(1) = (1d0 - R )
                ENDIF

                CALL PQueue(II)%INSERT( NODE )
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    
    DEALLOCATE( TRACER_TIME_SER )
    endTimer = MPI_Wtime()
    WRITE(*,105) myRank, endTimer - startTimer
    WRITE(*,*) ' Rank ', myRank, ' saw R range from ', RMIN, ' to ', RMAX

    II = 0
    IF (myRank == ROOT) THEN
       DO I = 1,numPoints
          IF ( I .GE. myPointMin ) THEN
             N = I - myPointMin + 1
             L = PQueue(N)%SIZE()
             J = MAX( L / 2, 1 )
             IF ( L .GT. 0_i8 ) THEN
                WRITE(*,*) 'Pqueue(', I, ') = ['
                call PQueue(N)%PEEK(1, NODE)
                WRITE(*,*), NODE, '...'
                call PQueue(N)%PEEK(int(J,4), NODE)
                WRITE(*,*), NODE, '...'
                call PQueue(N)%PEEK(INT(L,4), NODE)
                WRITE(*,*), NODE, '](', int(L, 4), ')'
                II = II + 1
             ENDIF
          ENDIF
          IF ( II .GE. 10 ) EXIT
       ENDDO
    ENDIF

    !
    !
    !
    ! At this point we have a priority queue for each point with all the
    ! R values to other points. We need to start combining the minimum
    ! R value points by finding out which MPI process has the smallest
    ! R value and then popping that and combining that cluster with
    ! the I and recomputing all the cluster R values
    !
    ! We need to iteratre as many times ars there are points because
    ! each iteration decreases the number of clusters by 1. We start
    ! with [numPoints] clusters and want to end with 1. 
    DO L = 1, numPoints
       RMIN = 999d99
       ! Iterate over all my priority Queues
       DO II = 1, numPoints
          IF ( II .LT. myPointMin ) CYCLE
          
          I = II - myPointMin + 1
          
          IF ( PQueue(I)%SIZE() .EQ. 0 ) CYCLE
          
          ! Check the min R value
          call PQueue(I)%PEEK(1, NODE)

          IF ( NODE(1) < RMIN ) THEN
             RMIN = NODE(1)
             J = II
          ENDIF
       ENDDO

       call PQueue(J - myPointMin + 1)%PEEK(1, NODE)
       minNode = (/ NODE(1), real(myRank + J, 8) /)

       call MPI_AllReduce( minNode, globalMinNode, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, &
            MPI_COMM_WORLD, ierr)

       ! There is a distinct possibility that real*8 doesn't have enough precision to
       ! do this correctly. 
       minRank = INT( globalMinNode(2) / numPoints )       
       
       ! Now the globalMinNode contains the R value and Rank of the global lowest R value
       IF ( minRank .eq. myRank) THEN
          ! for now we'll just pop the lowest one off and print it out
          call PQueue(J)%POP(NODE)
          call MPI_SEND( J, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr )
       ENDIF
       IF ( myRank .eq. ROOT ) THEN
          call MPI_RECV( J, 1, MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, ierr )
       ENDIF
    ENDDO
    
        
    call MPI_FINALIZE(ierr)
    
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
  
