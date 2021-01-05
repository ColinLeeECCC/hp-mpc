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

    ! data file variables
    integer                    :: forecastHour
    character(256)         :: dataFileName
    integer                    :: ncId, variableId
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIds
    integer                    :: grid_ni, grid_nj, grid_nt

    ! clustering information
    integer(kind=8)            :: numPoints, numPairs

    ! mpi information
    integer, parameter         :: ROOT = 0
    integer                    :: myRank, numProcs
    integer                    :: ierr
    integer                    :: myPointsStart, myPoints
    integer(kind=8)            :: myPairsStart, myPairsEnd
    integer,allocatable,dimension(:) &
         :: myPairsI, myPairsJ
    integer                    :: status(MPI_STATUS_SIZE)
    real*8                     :: startTimer, endTimer
    ! miscellanceous variables
    integer                    :: I, J, L, M, N
    integer(kind=8)            :: IND, IND1
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
103 FORMAT(A, '/', i4.4, i2.2, i2.2, i2.2 '_', i6.6, 'p.netcdf4.compressed')
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

    ! pair indices in point space (not grid space)
    ALLOCATE(myPairsI(myPairsEnd - myPairsStart + 1))
    ALLOCATE(myPairsJ(myPairsEnd - myPairsStart + 1))
    ! save the subscripts of the upper-triangular dissimilarity matrix
    
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
    
104 FORMAT(' Rank ', i3, ' has pairs ', i6, ' to ', i6, ' going from I,J = ', i3, &
         ', ', i3, ' to I,J = ', i3, ', ', i3)
105 FORMAT(' Rank ', i3, ' took ', f12.4, ' seconds.')
106 FORMAT(' Rank ', i3, ' has pairs ', i9, ' to ', i9, ' going from I,J = ', i6, &
         ', ', i6, ' to I,J = ', i6, ', ', i6)
107 FORMAT(' Rank ', i3, ' has pairs ', i14, ' to ', i14, ' going from I,J = ', i9, &
         ', ', i9, ' to I,J = ', i9, ', ', i9)

    call MPI_FINALIZE(ierr)
    
CONTAINS

   LOGICAL FUNCTION GREATER1( NODE1, NODE2 )
      DOUBLE  PRECISION, INTENT(IN) :: NODE1(:), NODE2(:)
      GREATER1 = NODE1(1) < NODE2(1)
   END FUNCTION GREATER1

   LOGICAL FUNCTION GREATER2( NODE1, NODE2 )
      DOUBLE PRECISION, INTENT(IN) :: NODE1(:), NODE2(:)
      GREATER2 = NODE1(2) < NODE2(2)
   END FUNCTION GREATER2

END PROGRAM CLUSTER
