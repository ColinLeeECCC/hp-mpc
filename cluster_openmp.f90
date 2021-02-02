  program cluster

    ! version 2.0
    ! Paul Makar, September 21, 2016
    ! Joana Soares, August 14, 2017
    ! Colin Lee,    January 4, 2021

    ! Find full hierarchical clustering information for gridded 
    ! NO2 data based on 1-R metric.

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

    integer                    :: I, J, K, N, M
    integer                    :: IDAY, IHR, HRS_SINCE_START
    integer                    :: ierr

    real*8                     :: R, RMIN, RMAX, Nt

    integer*8                  :: countStart, countEnd, countRate, countMax
    integer                    :: uid
    logical                    :: saveDissMatrix, loadDissMatrix
    
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
       STOP
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
    
    write(*,*) 'dims', grid_ni, grid_nj, grid_nt
    write(*,*) 'numPoints = ', numPoints
    
    WRITE(dataFileName,107) trim(outDir), start_year, start_mon,     &
         start_day, end_year, end_mon, end_day, forecastHour,  &
         grid_ni, grid_nj
    saveDissMatrix = .true.
    uid = 15
    WRITE(*,*) ' Checking if ', trim(dataFileName), ' exists...'
    inquire(file=dataFileName, exist=loadDissMatrix)
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

    call SYSTEM_CLOCK( countStart, countRate, countMax )
    
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
    ALLOCATE(TRACER_XYSUM( numPoints, numPoints ))

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
          ierr = nf90_opeN( trim(dataFileName), NF90_NOWRITE, ncId )
          IF (ierr < 0) THEN
             WRITE(*,*) 'Error while opening ', trim(dataFileName)
             STOP
          ENDIF

          ierr = nf90_inq_varid( ncId, trim(gemMachFieldName), variableId )
          ierr = nf90_inquire_variable( ncId, variableId, dimIds = dimIds )
          ierr = nf90_get_var( ncid, variableId, buffer)

          ierr = nf90_close( ncId )

          HRS_SINCE_START = 24 * ( IDAY - 1 ) + IHR
          !$OMP PARALLEL DO                      &
          !$OMP DEFAULT( SHARED )                &
          !$OMP PRIVATE( I, J, GRIDI, GRIDJ )    &
          !$OMP PRIVATE( TRACER_TMP, GRIDI2 )    &
          !$OMP PRIVATE( GRIDJ2, TRACER_TMP2 )
          DO I = 1,numPoints
             ! GRIDI = (I - 1) / GRID_NJ + 1
             ! GRIDJ = (I - 1) - (GRIDI - 1) * GRID_NJ + 1
             GRIDJ = (I - 1) / GRID_NJ + 1
             GRIDI = (I - 1) - (GRIDJ - 1) * GRID_NJ + 1
             TRACER_TMP = BUFFER( GRIDI, GRIDJ )
             TRACER_SUM( I )   = TRACER_SUM( I )   + TRACER_TMP
             TRACER_SQSUM( I ) = TRACER_SQSUM( I ) + TRACER_TMP * TRACER_TMP
             DO J = 1,numPoints
                ! GRIDI2 = (J - 1) / GRID_NJ + 1
                ! GRIDJ2 = (J - 1) - ( GRIDI2 - 1 ) * GRID_NJ + 1
                GRIDJ2 = (J - 1) / GRID_NJ + 1
                GRIDI2 = (J - 1) - ( GRIDJ2 - 1 ) * GRID_NJ + 1
                TRACER_TMP2 = BUFFER( GRIDI2, GRIDJ2 )

                TRACER_XYSUM( I, J ) = TRACER_XYSUM( I, J ) + &
                     TRACER_TMP * TRACER_TMP2
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          
       ENDDO
    ENDDO

    call SYSTEM_CLOCK( countEnd, countRate, countMax )
    WRITE(*,*) 'Loading data and precalculating took ', &
         real(countEnd - countStart) / real(countRate), ' sec'
    countStart = countEnd

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
    WRITE(*,*) ' Allocating pQueue with len=', numPoints
    ALLOCATE(PQueue(numPoints))
    ALLOCATE(NODES(numPoints, numPoints,2))
    ALLOCATE(LIVE(numPoints))
    ALLOCATE(clusterSize(numPoints))
    WRITE(*,*) ' Allocated'

    ! I have an inkling that doing this in the OMP loop is causing
    ! a runtime error
    DO N = 1,numPoints
       ! Create min-heap for this node. Will be zero if there are no edges
       CALL PQueue(N)%INIT( int(numPoints,4), 2, GREATER1 )
    ENDDO
    
    !$OMP PARALLEL DO               &
    !$OMP DEFAULT( SHARED )         &
    !$OMP PRIVATE( N, I, SXX, SYY ) &
    !$OMP PRIVATE( SXY, R )         &
    !$OMP REDUCTION(MAX: RMAX)      &
    !$OMP REDUCTION(MIN: RMIN)
    DO N = 1,numPoints
       DO I = 1,numPoints
          IF ( I .eq. N ) CYCLE ! don't compute self-similarity
          ! compute the 1-R for this pair of clusters
          SXY = tracer_xysum(N, I) * Nt - &
               tracer_sum(N) * tracer_sum(I)
          SXX = tracer_sqsum(I) * Nt - &
               ( tracer_sum(I) * tracer_sum(I))
          SYY = tracer_sqsum(N) * Nt - &
               ( tracer_sum(N) * tracer_sum(N) )

          IF ( SXX .eq. 0d0 .or. SYY .eq. 0d0 ) THEN
             R = 0d0
             WRITE(*,*) 'Got bad statistics at ', N, ',', I
             WRITE(*,*) 'Ex = ', tracer_sum(N), 'Ey = ', tracer_sum(I)
             WRITE(*,*) 'Ex2 = ', tracer_sqsum(N), 'Ey2 = ', tracer_sqsum(I)
             WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
          ELSE
             R = SXY / SQRT( SXX * SYY )
             IF ( R < -1d0 .or. R > 1d0 ) THEN
                WRITE(*,*) 'Got weird R at ', N, ',', I
                WRITE(*,*) ' R = ', R
                WRITE(*,*) 'Exy = ', tracer_xysum(N, I), ' Nt = ', Nt
                WRITE(*,*) 'Ex = ', tracer_sum(N), 'Ey = ', tracer_sum(I)
                WRITE(*,*) 'Ex2 = ', tracer_sqsum(N), 'Ey2 = ', tracer_sqsum(I)
                WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
             ENDIF
             IF (R < RMIN) RMIN = R
             IF (R > RMAX) RMAX = R
          ENDIF

          NODES(N,I,:) = (/ -1d0 * (1d0 - R), real(I,8) /)

          CALL PQueue(N)%INSERT( NODES(N,I,:) )
          heapIdx(N,I) = PQueue(N)%SIZE()
       ENDDO
       LIVE(N) = .true.
       clusterSize(N) = 1

    ENDDO
    !$OMP END PARALLEL DO
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
    
    call SYSTEM_CLOCK( countStart, countRate, countMax )

    ! Need to iterate this many times to get to a single cluster
    DO K = 1, numPoints-1
       
       ! find the queue with the min value
       RMIN = 999d99
       K1 = -1
       DO I = 1,numPoints
          IF (.not. LIVE(I)) CYCLE
          call PQueue(I).PEEK( 1, NODE )
          IF ( NODE(1) .le. RMIN ) THEN
             NODE1 = NODE
             RMIN = NODE1(1)
             K1 = I
          ENDIF
       ENDDO
       ! now K1 is the cluster with the smallest dissimilarity
       K2 = NODE1(2)
       IF (K1 .eq. K2 .or.      &
            .not. LIVE(K1) .or. &
            .not. LIVE(K2) ) THEN
          WRITE(*,*) 'Something has gone wrong. Print all queues'
          DO I = 1,numPoints
             IF (.not. LIVE(I)) CYCLE
             WRITE(*,*) 'Printing queue ', I
             DO WHILE( PQueue(I)%SIZE() .ne. 0)
                call PQueue(I).POP( NODE )
                WRITE(*,*), NODE
             ENDDO
          ENDDO
          STOP 1
       ENDIF

       IF (numPoints .lt. 1000) THEN
          WRITE(*,106) k1, k2, NODE1(1)
       ELSEIF (numPoints .lt. 1000000) THEN
          WRITE(*,105) k1, k2, NODE1(1)
       ELSE
          WRITE(*,104) k1, k2, NODE1(1)
       ENDIF
          
104    FORMAT('Clustering ', i9, ' and ', i9, ' (1-R) = ', F9.6)
105    FORMAT('Clustering ', i6, ' and ', i6, ' (1-R) = ', F9.6)
106    FORMAT('Clustering ', i3, ' and ', i3, ' (1-R) = ', F9.6)
       
       LIVE(K2) = .false.
       call PQueue(K1)%CLEAR()

       !$OMP PARALLEL DO                    &
       !$OMP DEFAULT( SHARED )              &
       !$OMP PRIVATE( I, NODE1, NODE2)      &
       !$OMP PRIVATE( M, R )
       DO I = 1,numPoints
          IF (.not. LIVE(I)) CYCLE
          IF (I .eq. K1) CYCLE
          
          NODE1(:) = NODES(I,K1,:)
          M = PQueue(I)%SIZE()
          ! call PQueue(I)%DELETE( NODE1 )
          call PQueue(I)%DELETE( K=heapIdx(I,K1) )
          IF (M - PQueue(I)%SIZE() .ne. 1) THEN
             WRITE(*,*) ' Failed to delete node1 ', NODE1, ' from queue ', I
             WRITE(*,*) ' M = ', M, ' size = ', PQueue(I)%SIZE()
          ENDIF
          NODE2(:) = NODES(I,K2,:)
          M = PQueue(I)%SIZE()
          ! call PQueue(I)%DELETE( NODE2 )
          call PQueue(I)%DELETE( K=heapIdx(K1,I) )
          IF (M - PQueue(I)%SIZE() .ne. 1)  THEN
             WRITE(*,*) ' Failed to delete node2 ', NODE2, ' from queue ', I
             WRITE(*,*) ' M = ', M, ' size = ', PQueue(I)%SIZE()
          ENDIF
          
          
          ! R = (clusterSize(K1) * NODE1(1) + clusterSize(K2) * NODE2(1)) &
          !      / ( clusterSize(K1) + clusterSize(K2) )
          R = MIN( NODE1(1), NODE2(1) )

          NODES( I,K1,1) = R
          NODES(K1, I,1) = R
          call PQueue( I)%INSERT( NODES( I,K1,:) )
          heapIdx( I,K1) = PQueue(I)%SIZE()
          !$OMP CRITICAL
          call PQueue(K1)%INSERT( NODES(K1, I,:) )
          heapIdx(K1,I) = PQueue(K1)%SIZE()
          !$OMP END CRITICAL
       ENDDO
       !$OMP END PARALLEL DO
       clusterSize(K1) = clusterSize(K1) + clusterSize(K2)
       clusterSize(K2) = 0
    ENDDO
    call SYSTEM_CLOCK( countEnd, countRate, countMax )
    WRITE(*,*) 'Clustering took ', &
         real(countEnd - countStart) / real(countRate), ' sec'
    countStart = countEnd


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
