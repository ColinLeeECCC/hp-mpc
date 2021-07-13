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
  character(256)             :: arglist, arg
  integer                    :: fu_in

  ! input parameters
  logical                    :: is_aircraft_data
  logical                    :: only_save_dissMat,tmpLogical
  character(256)             :: inDir, outDir
  integer                    :: start_year, start_mon, start_day
  integer                    ::   end_year,   end_mon,   end_day
  integer                    :: useFractionOfRegion
  character(8)               :: gemMachFieldName
  character(8)               :: fieldName
  character(256)             :: ncFieldName
  integer                    :: strLen
  integer                    :: start_date, end_date, start_vec(8), end_vec(8)
  integer                    :: current_date, current_year, current_mon, current_day
  real                       :: start_julian, end_julian

  ! data file variables
  integer                    :: forecastHour
  character(256)             :: dataFileName, acFileName
  integer                    :: ncId, variableId
  integer, dimension(NF90_MAX_VAR_DIMS) &
       :: dimIds
  integer                    :: grid_ni, grid_nj, grid_nt
  integer                    :: gridi, gridj, gridi2, gridj2
  integer                    :: numPoints, numSpectral
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

  integer,allocatable,dimension(:,:) &
                             :: clusterPairs
  real*4,allocatable,dimension(:) &
                             :: clusterDissimilarities
  
  integer                    :: I, J, K, N, M
  integer                    :: IDAY, IHR, HRS_SINCE_START
  integer                    :: ierr

  real*8                     :: R, RMIN, RMAX, RTMP, Nt

  integer*8                  :: countStart, countEnd, countRate, countMax
  integer                    :: uid
  logical                    :: saveDissMatrix, loadDissMatrix

  ! mpi variables
  integer, parameter         :: ROOT = 0
  integer                    :: myRank, numProcs
  integer                    :: status(MPI_STATUS_SIZE), readCount
  integer                    :: nodeType, dissMatRowType, dissMatBlockType
  integer(KIND=MPI_OFFSET_KIND):: disp
  integer                    :: mpi_uid

  ! mpi clustering variables
  integer                    :: II,JJ
  integer                    :: myPoints, numClustersThisNode
  integer                    :: clusterPrintPoint
  integer,allocatable,dimension(:) &
                             :: myClusters,clusterRanks
  real*8                     :: KK1, KK2
  logical                    :: K1ThisRank, K2ThisRank
  integer                    :: K1Rank, K2Rank

  ! Linkage variables
  integer                    :: linkage
  real*8                     :: D_IJ, alpha1, alpha2
  real*8                     :: beta, gamma
  

  ! MPI timing variables
  real*8                     :: startTimer, endTimer

  ! variables for tiling calculating dissimilarity matrix
  logical                    :: tileDissMat
  integer                    :: dissN, dissM
  integer                    :: dissSizeI, dissSizeJ
  integer                    :: dissStartI, dissStartJ
  real*8,allocatable,dimension(:,:) &
                             :: dissMat

  ! Variables for dissilimarity matrix file striping
  integer                    :: numStripes
  
  ! debugging
  logical,allocatable,dimension(:) :: visited 

  ! begin program

  ! initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numProcs, ierr)
  if (myRank == ROOT) WRITE(*,*) 'MPI Initizlied with ', numProcs, ' processors'

  ! Set default parameters
  arglist = 'arglist.lst'
  only_save_dissMat = .false.

  I = 1
  do while ( I .le. command_argument_count())
     call get_command_argument(i, arg)
     select case (arg)
     case ('-s', '--only-save')
        only_save_dissMat = .true.
     case ('f', '--argfile')
        I = I + 1
        call get_command_argument(i, arglist)
     case ('a', '--aircraft')
        is_aircraft_data = .true.
     case default
        IF ( I .eq. command_argument_count() ) THEN
           arglist = trim( arg )
        ELSE
           print '(2a, /)', 'unrecognised command-line option: ', arg
           stop
        ENDIF
     end select
     I = I + 1
  ENDDO
  

  WRITE(*,*) '  Reading arguments from file ''', trim(arglist), ''''
  ! Read in argument list from arglist.lst file
  fu_in = 15
  open(unit=fu_in, file=trim(arglist), status='old', iostat=ierr)
  IF (ierr .ne. 0) THEN
     WRITE(6,*) ' Couldn''t open "', trim(arglist), '"'
     call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
  ENDIF

  if ( .not. is_aircraft_data ) then
     READ(fu_in, '(a256)') inDir ! data where compressed netcdf GEM-MACH data is stored
     READ(fu_in, '(i4,2(1x,i2.2))') start_year, start_mon, start_day
     READ(fu_in, '(i4,2(1x,i2.2))')   end_year,   end_mon,   end_day
     READ(fu_in, *) linkage
     READ(fu_in, *) useFractionOfRegion
     READ(fu_in, '(a8)') gemMachFieldName
     READ(fu_in, '(a8)') fieldName
  else
     READ(fu_in, '(a256)') acFileName ! data where compressed netcdf GEM-MACH data is stored
     READ(fu_in, *) linkage
     READ(fu_in, *) useFractionOfRegion
     READ(fu_in, '(a256)') ncFieldName
  end if
  READ(fu_in, '(a256)', iostat = ierr) outDir

  if (ierr .eq. 0) THEN
     READ(fu_in, '(l1)', iostat = ierr) tmpLogical
     if ( ierr .eq. 0 ) THEN
        only_save_dissMat = tmpLogical
     ENDIF
  endif

  if (linkage .lt. 0 .or. linkage .gt. 6) THEN
     WRITE(*,*) 'Linkage must be between 0 and 6, inclusive, not ', linkage
     call MPI_Abort(MPI_COMM_WORLD, 7, IERR)
  endif

  ! make sure outDir ends in a '/'
  strLen = len_trim(outDir)
  IF (outDir(strLen:strLen) .ne. '/') THEN
     outDir = trim(outDir) // '/'
  ENDIF

  if (.not. is_aircraft_data) then
     strLen = len_trim(inDir)
     IF (inDir(strLen:strLen) .ne. '/') THEN
        inDir = trim(inDir) // '/'
     ENDIF


     WRITE(*,*) ' Reading compressed netcdf data from ', trim(inDir)
     WRITE(*,101) start_year, start_mon, start_day, end_year, end_mon, end_day
101  FORMAT(' Doing ', i4, '-', i2, '-', i2, ' to ', i4, '-', i2, '-', i2)
     WRITE(*,102) useFractionOfRegion
102  FORMAT(' Using only 1/', i3, ' of available GEM-MACH domain')
     WRITE(*,*) ' Using field ', fieldName, '(', gemMachFieldName, ')'
     WRITE(*,*) ' Ouputting to ', trim(outDir)

     forecastHour = 18

     ! Read in a single file to get model domain size
103  FORMAT(A, i4.4, i2.2, i2.2, i2.2 '_', i6.6, 'p.netcdf4.compressed')
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
     if (useFractionOfRegion .ne. -1) then
        grid_ni = useFractionOfRegion
        grid_nj = useFractionOfRegion
     end if
     numPoints = grid_ni * grid_nj
  else

     dataFileName = acFileName
     WRITE(*,*) ' Reading data from ', trim(dataFileName)
     WRITE(*,*) ' Using field ', trim(ncFieldName)
     WRITE(*,*) ' Ouputting to ', trim(outDir)
     ! Open the aircraft data file to get dimensions

     call check( nf90_open( trim(dataFileName), NF90_NOWRITE, ncId ) )
     
     WRITE(*,*) ' Using field ',  trim(ncFieldName)
     WRITE(*,*) ' Ouputting to ', trim(outDir)

     forecastHour = 18


     call check( nf90_inq_varid( ncId, trim(ncFieldName), variableId ) )
     call check( nf90_inquire_variable( ncId, variableId, dimIds = dimIds ) )
     call check( nf90_inquire_dimension( ncID, dimIds(1), len = grid_ni ) )
     call check( nf90_inquire_dimension( ncID, dimIds(2), len = grid_nj ) )

     if (useFractionOfRegion .ne. -1) then
        grid_ni = useFractionOfRegion
     end if
     numPoints = grid_ni
     numSpectral = grid_nj
     grid_nt = -1
     
     call check( nf90_close( ncId ) )
  end if

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

  ! for very large numPoints, our files get too big. The file system seems
  ! more than happy to go as big as we would like, but the problem is
  ! when we try to MPI_SET_VIEW we can overflow MPI_OFFSET_KIND. The limit seems
  ! to be
  call get_dissMat_filename(dataFileName, outDir, is_aircraft_data,     &
                            linkage, grid_ni, grid_nj,                  &
                            start_year, start_mon, start_day,           &
                            end_year, end_mon, end_day, forecastHour   )
  
  ! in most cases, we will load the dissMatrix from file, either
  ! because it already exists, or because we computed it tilewise and
  ! saved it out to be read in by the appropriate nodes
  loadDissMatrix = .true.
  uid = 15
  WRITE(*,*) ' Checking if ', trim(dataFileName), ' exists...'
  inquire(file=dataFileName, exist=saveDissMatrix)
  startTimer = MPI_Wtime()
  saveDissMatrix = .not. saveDissMatrix

  IF ( saveDissMatrix ) THEN
     WRITE(*,*) ' No. Computing dissimilarity matrix'

     ! In order to save memory, we can split the dissimilarity matrix up
     ! into roughly square tiles and each MPI node only computes that
     ! small tile

     tileDissMat = .false.
     ! This only works if our number of MPI nodes is a multiple of 2 or
     ! a square number (I think)
     dissN = INT(SQRT(REAL( numProcs )))
     if ( numProcs .ge. 4 .and. ( &
          MOD(numProcs, 2) .eq. 0 .or. &
          dissN**2 .eq. numProcs ) ) THEN
        tileDissMat = .true.
        dissM = INT(CEILING( REAL(numProcs) / REAL( dissN ) ))
        DO WHILE ( dissN * dissM .ne. numProcs )
           dissN = dissN - 1
           dissM = INT(CEILING( REAL(numProcs) / REAL( dissN ) ))
           IF ( dissN .le. 1 ) THEN
              WRITE(*,*) 'Please choose a number of MPI_THREADS that is a round or square number'
              call MPI_Abort( MPI_COMM_WORLD, 9, IERR )
           ENDIF
        ENDDO
        dissSizeI = CEILING( numPoints / REAL( dissN ) )
        dissSizeJ = CEILING( numPoints / REAL( dissM ) )
        dissStartI =      MOD(myRank,  dissN)  * dissSizeI + 1
        dissStartJ = INT(REAL(myRank / dissN)) * dissSizeJ + 1
        IF ( dissStartI + dissSizeI .gt. numPoints ) THEN
           dissSizeI = numPoints - dissStartI + 1
        END IF
        IF ( dissStartJ + dissSizeJ .gt. numPoints ) THEN 
           dissSizeJ = numPoints - dissStartJ + 1
        END IF
124     FORMAT('Node ', i2, ' I = ', i4, ' to ', i4, ', J = ', i4, ' to ', i4)
127     FORMAT('Node ', i2, ' I = ', i7, ' to ', i7, ', J = ', i7, ' to ', i7)
128     FORMAT('Node ', i2, ' I = ', i10, ' to ', i10, ', J = ', i10, ' to ', i10)
        IF ( numPoints .lt. 1000 ) THEN
           WRITE(*,124) myRank, dissStartI, dissStartI + dissSizeI - 1, &
                dissStartJ, dissStartJ + dissSizeJ - 1
        ELSEIF ( numPoints .lt. 1000000 ) THEN
           WRITE(*,127) myRank, dissStartI, dissStartI + dissSizeI - 1, &
                dissStartJ, dissStartJ + dissSizeJ - 1
        ELSE
           WRITE(*,128) myRank, dissStartI, dissStartI + dissSizeI - 1, &
                dissStartJ, dissStartJ + dissSizeJ - 1
        ENDIF
        
     ENDIF


     ALLOCATE(BUFFER(GRID_NI, GRID_NJ))

     IF ( tileDissMat ) THEN
        ALLOCATE(  TRACER_SUM( dissSizeI + dissSizeJ ))
        ALLOCATE(TRACER_SQSUM( dissSizeI + dissSizeJ ))
        ALLOCATE(TRACER_XYSUM( dissSizeI,  dissSizeJ ))
     ELSE
        ALLOCATE(  TRACER_SUM( numPoints ))
        ALLOCATE(TRACER_SQSUM( numPoints ))
        ALLOCATE(TRACER_XYSUM( numClustersThisNode, numPoints ))
     ENDIF

     TRACER_SUM = 0d0
     TRACER_SQSUM = 0d0
     TRACER_XYSUM = 0d0

     if (is_aircraft_data) then
        dataFileName = acFileName
        write(*,*) 'Opening ', trim(dataFileName)
        ierr = nf90_open( trim(dataFileName), NF90_NOWRITE, ncId)!, & 
!             comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)
        IF (ierr < 0) THEN
           WRITE(*,*) 'Error while opening ', trim(dataFileName), ' ierr = ', ierr
           call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
        ENDIF
        ierr = nf90_inq_varid( ncId, trim(ncFieldName), variableId )
        IF (ierr < 0) THEN
           WRITE(*,*) 'Error getting varid ', trim(ncFieldName)
           call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
        ENDIF
        ierr = nf90_inquire_variable( ncId, variableId, dimIds = dimIds )
        IF (ierr < 0) THEN
           WRITE(*,*) 'Error getting dimids '
           call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
        ENDIF
        write(*,*) 'Reading ', trim(ncFieldName)
        ierr = nf90_get_var( ncid, variableId, buffer)
        IF (ierr < 0) THEN
           WRITE(*,*) 'Error getting data'
           call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
        ENDIF
        write(*,*) 'Done reading ', trim(ncFieldName)
        ierr = nf90_close( ncId )
        IF (ierr < 0) THEN
           WRITE(*,*) 'Error closing file.'
           call MPI_Abort(MPI_COMM_WORLD, 2, IERR)
        ENDIF
        write(*,*) 'Closed ', trim(dataFileName)
        
        IF ( useFractionOfRegion .ne. -1 ) then
           DO I = 1,numPoints
              WRITE(*,1111) BUFFER(I,1:10)
           END DO
1111       FORMAT(G11.4, 9(', ',G11.4))
        endif
        numTimesteps = numSpectral

        DO IHR = 1, numSpectral
           !write(*,*) 'IHR = ', IHR

           IF ( tileDissMat ) THEN
              !$OMP PARALLEL DO                       &
              !$OMP DEFAULT( SHARED )                 &
              !$OMP PRIVATE( II, I, J, GRIDI, GRIDJ ) &
              !$OMP PRIVATE( TRACER_TMP, GRIDI2 )     &
              !$OMP PRIVATE( GRIDJ2, TRACER_TMP2 )
              DO II = 1,dissSizeI + dissSizeJ
                 IF ( II .le. dissSizeI ) THEN
                    I = II - 1 + dissStartI
                 ELSE
                    I = II - 1 - dissSizeI + dissStartJ
                 ENDIF
                 TRACER_TMP = BUFFER( I, IHR )

                 TRACER_SUM(II)   = TRACER_SUM(II)   + TRACER_TMP
                 TRACER_SQSUM(II) = TRACER_SQSUM(II) + TRACER_TMP * TRACER_TMP

                 IF ( II .le. dissSizeI ) THEN
                    DO JJ = 1,dissSizeJ
                       J = JJ - 1 + dissStartJ
                       TRACER_TMP2 = BUFFER( J, IHR )

                       TRACER_XYSUM( II, JJ ) = TRACER_XYSUM( II, JJ ) + &
                            TRACER_TMP * TRACER_TMP2
                    ENDDO
                 ENDIF
              ENDDO
              !$OMP END PARALLEL DO

           ELSE
              !$OMP PARALLEL DO                       &
              !$OMP DEFAULT( SHARED )                 &
              !$OMP PRIVATE( II, I, J, GRIDI, GRIDJ ) &
              !$OMP PRIVATE( TRACER_TMP, GRIDI2 )     &
              !$OMP PRIVATE( GRIDJ2, TRACER_TMP2 )
              DO II = 1,numClustersThisNode
                 I = myClusters(II)
                 ! GRIDI = (I - 1) / GRID_NJ + 1
                 ! GRIDJ = (I - 1) - (GRIDI - 1) * GRID_NJ + 1
                 TRACER_TMP = BUFFER( I, IHR )

                 DO J = 1,numPoints
                    ! GRIDI2 = (J - 1) / GRID_NJ + 1
                    ! GRIDJ2 = (J - 1) - ( GRIDI2 - 1 ) * GRID_NJ + 1
                    TRACER_TMP2 = BUFFER( J, IHR )

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
                 TRACER_TMP = BUFFER( I, IHR )

                 TRACER_SUM( I )   = TRACER_SUM( I )   + TRACER_TMP
                 TRACER_SQSUM( I ) = TRACER_SQSUM( I ) + TRACER_TMP * TRACER_TMP
              ENDDO
              !$OMP END PARALLEL DO
           ENDIF ! ( tileDissMat )
        ENDDO
     else


        ! read in all the time data from start to end
        start_vec = (/start_year,start_mon,start_day,forecastHour,00,00,00,00/)
        end_vec   = (/end_year,    end_mon,  end_day,forecastHour,00,00,00,00/)
        call d2j(start_vec, start_julian, ierr)
        call d2j(  end_vec,   end_julian, ierr)
        numDays = int(end_julian - start_julian + 1)
        numTimesteps = numDays * 24

108  FORMAT(A, i8.8, i2.2 '_', i6.6, 'p.netcdf4.compressed')

     DO IDAY = 1, numDays
        call j2d(start_julian+IDAY*1.0, current_date, ierr)
        current_year = current_date / 10000
        current_mon  = current_date / 100 - current_year * 100
        current_day  = current_date - (current_year * 100 + current_mon) * 100
        WRITE(*,110) current_year, current_mon, current_day
110     FORMAT(' Loading data for day ', i4.4, '/', i2.2, '/', i2.2)

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
           IF ( tileDissMat ) THEN
              !$OMP PARALLEL DO                       &
              !$OMP DEFAULT( SHARED )                 &
              !$OMP PRIVATE( II, I, J, GRIDI, GRIDJ ) &
              !$OMP PRIVATE( TRACER_TMP, GRIDI2 )     &
              !$OMP PRIVATE( GRIDJ2, TRACER_TMP2 )
              DO II = 1,dissSizeI + dissSizeJ
                 IF ( II .le. dissSizeI ) THEN
                    I = II - 1 + dissStartI
                 ELSE
                    I = II - 1 - dissSizeI + dissStartJ
                 ENDIF
                 GRIDJ = (I - 1) / GRID_NJ + 1
                 GRIDI = (I - 1) - (GRIDJ - 1) * GRID_NJ + 1
                 TRACER_TMP = BUFFER( GRIDI, GRIDJ )

                 TRACER_SUM(II)   = TRACER_SUM(II)   + TRACER_TMP
                 TRACER_SQSUM(II) = TRACER_SQSUM(II) + TRACER_TMP * TRACER_TMP

                 IF ( II .le. dissSizeI ) THEN
                    DO JJ = 1,dissSizeJ
                       J = JJ - 1 + dissStartJ
                       GRIDJ2 = (J - 1) / GRID_NJ + 1
                       GRIDI2 = (J - 1) - ( GRIDJ2 - 1 ) * GRID_NJ + 1
                       TRACER_TMP2 = BUFFER( GRIDI2, GRIDJ2 )

                       TRACER_XYSUM( II, JJ ) = TRACER_XYSUM( II, JJ ) + &
                            TRACER_TMP * TRACER_TMP2
                    ENDDO
                 ENDIF
              ENDDO
              !$OMP END PARALLEL DO

           ELSE
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
           ENDIF ! ( tileDissMat )
        ENDDO
     ENDDO
     end if ! if (is_aircraft_data)

     IF ( tileDissMat ) THEN
        write(*,*) ' Allocating dissMat '
        ALLOCATE( dissMat( dissSizeI, dissSizeJ ) )
        write(*,*) ' Allocated.'

        RMIN = 9.999d9
        RMAX = -9.999d9
        Nt = real(numTimesteps, 8)

        !$OMP PARALLEL DO               &
        !$OMP DEFAULT( SHARED )         &
        !$OMP PRIVATE( N, I, SXX, SYY ) &
        !$OMP PRIVATE( SXY, R, II )     &
        !$OMP REDUCTION(MAX: RMAX)      &
        !$OMP REDUCTION(MIN: RMIN)
        DO II = 1,dissSizeI
           N = II - 1 + dissStartI

           DO JJ = 1, dissSizeJ
              J = JJ - 1 + dissStartJ
              IF ( J .eq. N ) CYCLE ! don't compute self-similarity
              ! WRITE(*,*) '   ', I
              ! compute the 1-R for this pair of clusters
              SXY = tracer_xysum(II, JJ) * Nt - &
                   tracer_sum(II) * tracer_sum(JJ + dissSizeI)
              SXX = tracer_sqsum(JJ + dissSizeI) * Nt - &
                   ( tracer_sum(JJ + dissSizeI) * tracer_sum(JJ + dissSizeI))
              SYY = tracer_sqsum(II) * Nt - &
                   ( tracer_sum(II) * tracer_sum(II) )
              ! WRITE(*,*) '   ', SXX, ', ', SYY, ', ', SXY
              IF ( SXX .eq. 0d0 .or. SYY .eq. 0d0 ) THEN
                 R = -1d0
                 if (.not. is_aircraft_data) then
                    WRITE(*,*) 'Got bad statistics at ', N, '(', II, '), ', JJ
                    WRITE(*,*) 'Ex = ', tracer_sum(II), 'Ey = ', tracer_sum(JJ + dissSizeI)
                    WRITE(*,*) 'Ex2 = ', tracer_sqsum(II), 'Ey2 = ', tracer_sqsum(JJ + dissSizeI)
                    WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
                 end if
              ELSE
                 R = SXY / SQRT( SXX * SYY )
                 IF ( R < -1d0 .or. R > 1d0 ) THEN
                    WRITE(*,*) 'Got weird R at ',II, '(', N, '),', JJ
                    WRITE(*,*) ' R = ', R
                    WRITE(*,*) 'Exy = ', tracer_xysum(II, JJ), ' Nt = ', Nt
                    WRITE(*,*) 'Ex = ', tracer_sum(II), 'Ey = ', tracer_sum(JJ + dissSizeI)
                    WRITE(*,*) 'Ex2 = ', tracer_sqsum(II), 'Ey2 = ', tracer_sqsum(JJ + dissSizeI)
                    WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
                 ENDIF
                 IF (R < RMIN) RMIN = R
                 IF (R > RMAX) RMAX = R
              ENDIF
              dissMat(II,JJ) = ( 1d0 - R )

           ENDDO

        ENDDO
        !$OMP END PARALLEL DO

        ! Now we save out the tile of the dissimilarity matrix we computed.
        call get_dissMat_filename(dataFileName, outDir, is_aircraft_data,  &
             linkage, grid_ni, grid_nj, start_year, start_mon, start_day,  &
             end_year, end_mon, end_day, forecastHour)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(dataFileName), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, mpi_uid, ierr)
        IF (IERR .ne. MPI_SUCCESS) THEN
           WRITE(*,*) ' ERROR: ', IERR, ' Opening file ', trim(dataFileName)
           call MPI_Abort(MPI_COMM_WORLD, 4, IERR)
        ENDIF
        WRITE(*,*) 'Creating block type'
        call MPI_Type_vector(dissSizeI, dissSizeJ, numPoints, &
             MPI_DOUBLE_PRECISION, dissMatBlockType, ierr)
        call MPI_Type_commit(dissMatBlockType, ierr)
        ! I think multiplying by 8 for real*8 should be right...
        disp = ((INT(dissStartI, kind=MPI_OFFSET_KIND) - 1) * &
             INT(numPoints, kind=MPI_OFFSET_KIND) + &
             INT(dissStartJ, kind=MPI_OFFSET_KIND) - 1_MPI_OFFSET_KIND) &
             * 8_MPI_OFFSET_KIND
        
        WRITE(*,*) 'Setting view to start at ', disp
        call MPI_FILE_SET_VIEW(mpi_uid, disp, MPI_DOUBLE_PRECISION, &
             dissMatBlockType, 'native', MPI_INFO_NULL, ierr)

        WRITE(*,*) 'Creating row type'
        ! Now we need a data type for the individual row
        call MPI_Type_contiguous( dissSizeJ, MPI_DOUBLE_PRECISION, dissMatRowType, ierr )
        call MPI_Type_commit( dissMatRowType, ierr )            

        DO II = 1,dissSizeI
           call MPI_File_write( mpi_uid, dissMat(II,:), 1, dissMatRowType, status, ierr )
           IF (ierr .ne. MPI_SUCCESS) THEN
              WRITE(*,*) 'MPI_File_write retured ierr = ', ierr
              call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
           ENDIF
           call MPI_Get_count(status, MPI_DOUBLE_PRECISION, readCount, ierr)
           IF ( readCount .ne. dissSizeJ ) THEN
              WRITE(*,*) 'MPI_File_write wrote ', readCount, ' not ', dissSizeJ, &
                   ' at step ', II, ' of ', dissSizeI
              call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
           ENDIF
           WRITE(*,*) ' Wrote ', readCount, ' doubles to file.'
        ENDDO
        call MPI_File_Close( mpi_uid, ierr )
        WRITE(*,*) 'Finished writing dissimilarity matrix. File closed.'
        saveDissMatrix = .false.

        call MPI_Barrier( MPI_COMM_WORLD, ierr )
     ELSE
        loadDissMatrix = .false.

        RMIN = 9.999d9
        RMAX = -9.999d9
        Nt = real(numTimesteps, 8)

        ! I have an inkling that doing this in the OMP loop is causing
        ! a runtime error
        DO N = 1,numClustersThisNode
           ! Create min-heap for this node. Will be zero if there are no edges
           CALL PQueue(N)%INIT( int(numPoints,4), 2, GREATER1 )
        ENDDO

        WRITE(*,*) "Initialized all pQueues"
        !$OMP PARALLEL DO               &
        !$OMP DEFAULT( SHARED )         &
        !$OMP PRIVATE( N, I, SXX, SYY ) &
        !$OMP PRIVATE( SXY, R, II )     &
        !$OMP REDUCTION(MAX: RMAX)      &
        !$OMP REDUCTION(MIN: RMIN)
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
                 R = -1d0
                 if (.not. is_aircraft_data) then
                    WRITE(*,*) 'Got bad statistics at ', N, ',', I
                    WRITE(*,*) 'Ex = ', tracer_sum(N), 'Ey = ', tracer_sum(I)
                    WRITE(*,*) 'Ex2 = ', tracer_sqsum(N), 'Ey2 = ', tracer_sqsum(I)
                    WRITE(*,*) 'SXY = ', SXY, 'SXX = ', SXX, 'SYY = ', SYY
                 endif
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
              NODES(II,I,:) = (/ (1d0 - R), real(I,8) /)

112           FORMAT(' Inserting <', F6.3, ', ', I2,'> into queue ', I2)
              ! WRITE(*,112) NODES(II,I,1), INT(NODES(II,I,2)), II
              CALL PQueue(II)%INSERT( NODES(II,I,:), heapIdx(II,I) )
              ! WRITE(*,*) ' Done '
113           FORMAT(' At index(', i2, ',', i2, ') = ', i2)
              ! WRITE(*,113) II, I, heapIdx(II,I)
           ENDDO

        ENDDO
        !$OMP END PARALLEL DO
     ENDIF

     DEALLOCATE( TRACER_SUM )
     DEALLOCATE( TRACER_SQSUM )
     DEALLOCATE( TRACER_XYSUM )
     endTimer = MPI_Wtime()
     WRITE(*,*) ' Calclating dissimilarity matrix took ', endTimer - startTimer, ' seconds'
  END IF ! ( saveDissMatrix )

  if (only_save_dissMat) THEN
     if (myRank .eq. ROOT) THEN
        WRITE(*,*) ' Only saving dissimilarity matrix. Exiting'
     ENDIF

     call MPI_Finalize( IERR )
     STOP
     
  ENDIF

  IF ( loadDissMatrix ) THEN
     startTimer = MPI_Wtime()

     WRITE(*,*) 'Opening'

     call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(dataFileName), MPI_MODE_RDONLY, MPI_INFO_NULL, mpi_uid, ierr)
     IF (IERR .ne. MPI_SUCCESS) THEN
        WRITE(*,*) ' ERROR: Opening file ', trim(dataFileName)
        call MPI_Abort(MPI_COMM_WORLD, 4, IERR)
     ENDIF
     ! This is a little bit experimental at the moment. I think we should be able
     ! to set the view with this vector with COUNT=numClustersThisNode, BLOCKLENGTH=numPoints,
     ! STRIDE=numPoints*numProcs
     WRITE(*,*) 'Creating block type'
     call MPI_Type_vector(numClustersThisNode, numPoints, numPoints*numProcs, &
          MPI_DOUBLE_PRECISION, dissMatBlockType, ierr)
     call MPI_Type_commit(dissMatBlockType, ierr)
     ! I think multiplying by 8 for real*8 should be right...
     disp = myRank * numPoints * 8_MPI_OFFSET_KIND
     WRITE(*,*) 'Setting view'
     call MPI_FILE_SET_VIEW(mpi_uid, disp, MPI_DOUBLE_PRECISION, &
          dissMatBlockType, 'native', MPI_INFO_NULL, ierr)

     WRITE(*,*) 'Creating row type'
     ! Now we need a data type for the individual row
     call MPI_Type_contiguous( numPoints, MPI_DOUBLE_PRECISION, dissMatRowType, ierr )
     call MPI_Type_commit( dissMatRowType, ierr )

     RMIN = 9.999d9
     RMAX = -9.999d9

     clusterPrintPoint = int(numClustersThisNode / 10d0)
    
     DO II = 1, numClustersThisNode
        I = myClusters(II)
        IF ( myRank .eq. ROOT ) THEN
           IF ( MOD( II, clusterPrintPoint ) .eq. 1 ) THEN
              WRITE(*,*) '  Initializing pQ ', II, ' out of ', numClustersThisNode
           ENDIF
        ENDIF
        CALL PQueue(II)%INIT( int(numPoints,4), 2, GREATER1 )

        IF ( myRank .eq. ROOT .and. II .eq. numClustersThisNode ) THEN
           WRITE(*,*) ' pQ ', II, ' initialized.'
        ENDIF

        call MPI_File_Read( mpi_uid, NODES(II,:,1), 1, dissMatRowType, status, ierr )
        IF (ierr .ne. MPI_SUCCESS) THEN
           WRITE(*,*) 'MPI_File_Read retured ierr = ', ierr
           call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
        ENDIF
        call MPI_Get_count(status, MPI_DOUBLE_PRECISION, readCount, ierr)
        IF ( readCount .ne. numPoints ) THEN
           WRITE(*,*) 'MPI_File_Read read ', readCount, ' not ', numPoints, &
                ' at step ', II, '(', I, ')'
           call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
        ENDIF

        IF ( myRank .eq. ROOT .and. II .eq. numClustersThisNode ) THEN
           WRITE(*,*) ' Row ', II, ' read from file.'
        ENDIF
        
        DO J = 1, numPoints
           IF (I .eq. J) CYCLE

           NODES(II,J,2) = real(J,8)
           CALL PQueue(II)%INSERT( NODES(II,J,:), heapIdx(II,J) )
        ENDDO

        IF ( myRank .eq. ROOT .and. II .eq. numClustersThisNode ) THEN
           WRITE(*,*) ' All values inserted into pQ ', II
        ENDIF
        
        IF (I .eq. 1) THEN
           RTMP = MINVAL( NODES(II,2:numPoints,1) )
        ELSEIF (I .eq. numPoints) THEN
           RTMP = MINVAL( NODES(II,1:numPoints-1,1) )
        ELSE
           RTMP = MIN( MINVAL( NODES(II,1:(I-1), 1) ), &
                MINVAL( NODES(II, (I+1):numPoints,1) ) )
        ENDIF
        IF (RTMP < RMIN) &
             RMIN = RTMP
        IF (I .eq. 1) THEN
           RTMP = MAXVAL( NODES(II,2:numPoints,1) )
        ELSEIF (I .eq. numPoints) THEN
           RTMP = MAXVAL( NODES(II,1:numPoints-1,1) )
        ELSE
           RTMP = MAX( MAXVAL( NODES(II,1:(I-1), 1) ), &
                MAXVAL( NODES(II, (I+1):numPoints,1 ) ) )
        ENDIF
        IF (RTMP > RMAX) &
             RMAX = RTMP
     ENDDO
     call MPI_File_Close( mpi_uid, ierr )
     WRITE(*,*) 'Finished reading dissimilarity matrix. File closed.'

     endTimer = MPI_Wtime()
     WRITE(*,*) 'Loading dissimilarity matrix took ', endTimer - startTimer, ' seconds'
     
     RMIN = 1d0 - RMIN
     RMAX = 1d0 - RMAX
  ENDIF

  WRITE(*,*) 'R ranged from ', RMIN, ' to ', RMAX

  if (saveDissMatrix) THEN
     startTimer = MPI_Wtime()
     WRITE(*,*) 'Saving dissimilarity matrix...'
     call get_dissMat_filename(dataFileName,outDir, is_aircraft_data,   &
          linkage, grid_ni, grid_nj, start_year, start_mon, start_day,  &
          end_year, end_mon, end_day, forecastHour)

     WRITE(*,*) 'Opening'
     call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(dataFileName), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, mpi_uid, ierr)
     IF (IERR .ne. MPI_SUCCESS) THEN
        WRITE(*,*) ' ERROR: ', IERR, ' Opening file ', trim(dataFileName)
        call MPI_Abort(MPI_COMM_WORLD, 4, IERR)
     ENDIF
     WRITE(*,*) 'Creating block type'
     call MPI_Type_vector(numClustersThisNode, numPoints, numPoints*numProcs, &
          MPI_DOUBLE_PRECISION, dissMatBlockType, ierr)
     call MPI_Type_commit(dissMatBlockType, ierr)
     ! I think multiplying by 8 for real*8 should be right...
     disp = myRank * numPoints * 8_MPI_OFFSET_KIND
     WRITE(*,*) 'Setting view'
     call MPI_FILE_SET_VIEW(mpi_uid, disp, MPI_DOUBLE_PRECISION, &
          dissMatBlockType, 'native', MPI_INFO_NULL, ierr)

     WRITE(*,*) 'Creating row type'
     ! Now we need a data type for the individual row
     call MPI_Type_contiguous( numPoints, MPI_DOUBLE_PRECISION, dissMatRowType, ierr )
     call MPI_Type_commit( dissMatRowType, ierr )            

     DO II = 1,numClustersThisNode
        call MPI_File_write( mpi_uid, NODES(II,:,1), 1, dissMatRowType, status, ierr )
        IF (ierr .ne. MPI_SUCCESS) THEN
           WRITE(*,*) 'MPI_File_write retured ierr = ', ierr
           call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
        ENDIF
        call MPI_Get_count(status, MPI_DOUBLE_PRECISION, readCount, ierr)
        IF ( readCount .ne. numPoints ) THEN
           WRITE(*,*) 'MPI_File_write wrote ', readCount, ' not ', numPoints
           call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
        ENDIF
        WRITE(*,*) ' Wrote ', readCount, ' doubles to file.'
     ENDDO
     call MPI_File_Close( mpi_uid, ierr )
     WRITE(*,*) 'Finished writing dissimilarity matrix. File closed.'
     WRITE(*,*) 'Done.'
     endTimer = MPI_Wtime()
     WRITE(*,*) 'Saving dissimilarity matrix took ', endTimer - startTimer, ' seconds'
  ENDIF

  ALLOCATE( mpiNodes( 3, numProcs ) )
  call MPI_TYPE_CONTIGUOUS( 3, MPI_REAL, nodeType, ierr )
  call MPI_TYPE_COMMIT( nodeType, ierr )
  ALLOCATE(NODESK1(2,numPoints))
  ALLOCATE(NODESK2(2,numPoints))
  ALLOCATE(VISITED(numPoints))
  IF ( myRank .eq. ROOT ) THEN
     ALLOCATE(clusterPairs(numPoints - 1, 2))
     ALLOCATE(clusterDissimilarities(numPoints - 1))
  ENDIF

  startTimer = MPI_Wtime()

  ! to save time on some metrics, save the Lance-Williams coefficients
  ! in advance
  call init_linkage(linkage)

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
117     FORMAT('R = ', f7.4, ', K1 = ', I1)
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
     ! save the distance of the clustered pair
     D_IJ     = RMIN

     IF ( myRank .eq. ROOT ) THEN
        clusterPairs(K,1) = k1
        clusterPairs(K,2) = k2
        clusterDissimilarities(K) = NODE1(1)
     ENDIF
     IF (myRank .eq. ROOT) THEN
        IF (numPoints .lt. 1000) THEN
           WRITE(*,106) k1, k2, NODE1(1)
        ELSEIF (numPoints .lt. 1000000) THEN
           WRITE(*,105) k1, k2, NODE1(1)
        ELSE
           WRITE(*,104) k1, k2, NODE1(1)
        ENDIF
     ENDIF

104  FORMAT('Clustering ', i9, ' and ', i9, ' (1-R) = ', F9.6)
105  FORMAT('Clustering ', i6, ' and ', i6, ' (1-R) = ', F9.6)
106  FORMAT('Clustering ', i3, ' and ', i3, ' (1-R) = ', F9.6)

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

114     FORMAT('Cluster ', i2, ' LIVE? ', L1, ' .ne. K1? ', L1)
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
111        FORMAT('Cluster ', i2, ' deleting node ', 2f6.3, ' at index (', i2, ', ', i2, ') = ', i2)
           ! WRITE(*,111) I, NODE1, II, K1, PQueue(II)%INDEXAT(heapIdx(II,K1))
           call PQueue(II)%DELETE( K=heapIdx(II,K1), DNODE=NODE )
           IF (.not. ALL(NODE1 .eq. NODE) ) THEN
120           FORMAT('Step ', i5, ', k1 = ', i5, ', k2 = ', i5)
121           FORMAT('N1 = ', f7.4, ', ', i5)
122           FORMAT('N  = ', f7.4, ', ', i5)
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

           ! calculate the new distance metric between the combined clusters
           ! K1 and K2 and the old cluster II
           call calc_new_distance( linkage, R, NODE1(1), NODE2(1), D_IJ,     &
                                   clusterSize(K1), clusterSize(K2),         &
                                   clusterSize(I)                           )
           
           NODES(II,K1,1) = R
116        FORMAT('PQueue(',i2,')%INSERT( [', f7.4, ', ', i2, '])')
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

           ! calculate the new distance metric between the combined clusters
           ! K1 and K2 and the old cluster II
           call calc_new_distance( linkage, R, NODE1(1), NODE2(1), D_IJ,     &
                                   clusterSize(K1), clusterSize(K2),         &
                                   clusterSize(I)                           )

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
123  FORMAT('WARNING: Step', i5, ': Did not visit cluster ', i5, ', K1 = ', i5, ', K2 = ', i5)
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
118        FORMAT('WARNING: Step', i5,': Cluster ', i5, '(', i4, ') has size ', i5, ' not ', i5)
           WRITE(*,118) k, myClusters(II), II, M, N
        ENDIF
     ENDDO

     clusterSize(K1) = clusterSize(K1) + clusterSize(K2)
     clusterSize(K2) = 0
  ENDDO
  endTimer = MPI_Wtime()
  WRITE(*,*) 'Clustering took ', endTimer - startTimer,' sec'
  countStart = countEnd

  IF ( myRank .eq. ROOT ) THEN
     if (is_aircraft_data) then
        WRITE(dataFileName,129) trim(outDir), grid_ni, linkage
     else
        WRITE(dataFileName,126) trim(outDir), start_year, start_mon,     &
             start_day, end_year, end_mon, end_day, forecastHour,  &
             linkage, grid_ni, grid_nj
     end if
126  FORMAT(a,i4.4,i2.2,i2.2, '_', i4.4,i2.2,i2.2,'_', 2(i2.2, '_'), &
          i4.4, 'x', i4.4, '_clusters.dat')
129  FORMAT(a,i7.7,'_',i2.2, '_clusters.dat')
     open( uid, file=trim(dataFileName), form='formatted', action='write', &
          status='replace', iostat=ierr)
     write(uid,*) 'xx\tlki\tlkj'
     DO K = 1, numPoints - 1
        write(uid,125) clusterDissimilarities(K), clusterPairs(K,1), clusterPairs(K,2)
     ENDDO
125  FORMAT(f8.4, 10x, i5, 5x, i5)
     close(uid)
  ENDIF

  call MPI_FINALIZE(IERR)

CONTAINS
  subroutine get_dissMat_filename(fileName, outDir, is_aircraft_data, &
       linkage, grid_ni, grid_nj, start_year, start_mon, start_day,   &
       end_year, end_mon, end_day, forecastHour)

    character(256), intent(out) :: fileName
    character(256), intent(in)  :: outDir
    logical, intent(in)         :: is_aircraft_data
    integer, intent(in)         :: linkage
    integer, intent(in)         :: grid_ni, grid_nj
    integer, intent(in)         :: start_year, start_mon, start_day
    integer, intent(in)         :: end_year, end_mon, end_day
    integer, intent(in)         :: forecastHour

107 FORMAT(a,i4.4,i2.2,i2.2, '_', i4.4,i2.2,i2.2,'_', 2(i2.2, '_'), &
         i4.4, 'x', i4.4, '_mpio.bin')
108 FORMAT(a, i7.7,'_',i2.2, '_mpio.bin')

    if (is_aircraft_data) then
       WRITE(fileName,108) trim(outDir), grid_ni, linkage
    else
       WRITE(dataFileName,107) trim(outDir), start_year, start_mon,     &
            start_day, end_year, end_mon, end_day, forecastHour,  &
            linkage, grid_ni, grid_nj
    end if
  end subroutine get_dissMat_filename

  subroutine init_linkage(linkage)
    ! Set up the variables for calculating the combined cluster
    ! distance metric based on one of 7 linkage methods

    ! Integer from 0-6 indicating which linkage
    integer, intent(in) :: linkage

    ! Single linkage, complete linkage and average linkage
    ! all have non-Lance-Williams update methods which
    ! should be faster, so we just skip those. Centroid
    ! linkage and Ward's method have dynamic coefficients
    ! so we also won't initialize those here.
    ! 
    ! All the methods we're using the Lance-Williams
    ! equation for have gamma=0, so we're just gonna leave
    ! that term out.
    
       ! Weighted linkage
    if ( linkage .eq. 3 ) then
       alpha1 = 0.5
       alpha2 = 0.5
       beta   = 0.0
       ! Median linkage
    else if (linkage .eq. 5 ) then
       alpha1 = 0.5
       alpha2 = 0.5
       beta   = -0.25
    else if (linkage .lt. 0 .or. linkage .gt. 6) then
       WRITE(*,*) 'ERROR in INIT_LINKAGE. Linkage = ', linkage, ' not supported.'
       call MPI_Abort(MPI_COMM_WORLD, 7, IERR)
    endif
  end subroutine init_linkage

  subroutine calc_new_distance(linkage, D_IJK, D_IK, D_JK, D_IJ, N_I, N_J, N_K)
    ! update the distance between new cluster IJ and old cluster K from
    ! the old distances betwen I and K and J and K

    ! linkage method 0-6
    integer, intent(in)       :: linkage
    ! return value, new distance between IJ and K
    real*8, intent(out)        :: D_IJK
    ! old distances between I and K, J and K, and I and J
    real*8, intent(in)        :: D_IK, D_JK, D_IJ
    ! cluster sizes of I, J and K
    integer, intent(in)       :: N_I, N_J, N_K

    select case (LINKAGE)
       ! Single linkage
    case (0)
       D_IJK = MIN( D_IK, D_JK )
       ! Complete linkage
    case (1)
       D_IJK = MAX( D_IK, D_JK )
       ! Average linkage
    case (2)
       D_IJK =  ( N_I * D_IK + N_J * D_JK ) &
            / ( N_I + N_J )
       ! Weighted linkage
    case (3)
       D_IJK = alpha1 * D_IK + alpha2 * D_JK + &
               beta * D_IJ
       ! Centroid linkage
    case (4)
       alpha1 = N_I / ( N_I + N_J )
       alpha2 = N_J / ( N_I + N_J )
       beta   = - N_I * N_J / &
            ( N_I + N_J ) ** 2
       D_IJK = alpha1 * D_IK + alpha2 * D_JK + &
               beta * D_IJ
       ! Median linkage
    case (5)
       D_IJK = alpha1 * D_IK + alpha2 * D_JK + &
               beta * D_IJ
       ! Ward's linkage
    case (6)
       alpha1 = ( N_I + N_K ) / &
            ( N_I + N_J + N_K )
       alpha2 = ( N_J + N_K ) / &
            ( N_I + N_J + N_K )
       beta   = - N_K / &
            ( N_I + N_J + N_K )
       D_IJK = alpha1 * D_IK + alpha2 * D_JK + &
               beta * D_IJ
       ! Default to single linkage
    case default
       WRITE(*,*) 'ERROR in calc_new_distance. Linkage = ', linkage, ' not supported.'
       call MPI_Abort(MPI_COMM_WORLD, 7, IERR)
    end select
  end subroutine calc_new_distance
    
    
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
  
  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check

END PROGRAM CLUSTER
