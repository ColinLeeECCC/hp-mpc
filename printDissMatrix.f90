PROGRAM	PRINTDISSMATRIX

  implicit none

  real(kind=4),allocatable,dimension(:,:)  &
                     :: E
  real(kind=8),allocatable,dimension(:,:)  &
                     :: R
  integer(kind=1),allocatable,dimension(:) &
                     :: T
  integer(kind=1)    :: tmpi(8)
  real(kind=8)       :: tmpr
  equivalence (tmpi, tmpr)
  integer*8          :: i, j,l

  character(len=256) :: argStr, fileName, fmt
  integer            :: ftype, argLen, stat
  integer            :: uid, argCount
  logical            :: fileExist
  integer*8          :: nPoints

  argCount = 1
  uid = 15
  ! first get the filename we want to dump
  call get_command_argument( argCount, argStr, argLen, stat)
  fileName = argStr
  inquire( file=trim(fileName), exist=fileExist )
  if (stat .ne. 0 .or. .not. fileExist) &
       call PRINT_USAGE()

  argCount = argCount + 1
  ! first get the size of the side of the tile we used
  ! from the GEM-MACH background fields
  call get_command_argument( argCount, argStr, argLen, stat)
  READ( argStr(1:argLen), '(i4)' ) npoints ! nside
  if (npoints .le. 0 .or. stat .ne. 0) &
       call PRINT_USAGE()

  argCount = argCount + 1
  ! first get the real type (ie did the writing file use
  ! real*4 or real*8
  call get_command_argument( argCount, argStr, argLen, stat)
  READ( argStr(1:argLen), '(i4)' ) ftype

  if ( stat .ne. 0 .or. .not. (ftype .eq. 4 .or. &
       ftype .eq. 8 .or. &
       ftype .eq. 9 ) ) &
       call PRINT_USAGE()

  ! NPOINTS = NSIDE * NSIDE
  if ( ftype .eq. 9 ) then
     open(uid, file=trim(fileName), status='old', form='unformatted', &
          access='direct', recl=NPOINTS*8, convert='BIG_ENDIAN')
  else
     open(uid, file=trim(fileName), status='old', form='unformatted')
  endif
  IF (FTYPE .eq. 4) THEN
     ALLOCATE(E(NPOINTS, NPOINTS))
     DO I = 1,NPOINTS
        DO J = 1, NPOINTS
           READ(uid) E(I,J)
        ENDDO
     ENDDO
     
  ELSE IF ( FTYPE .eq. 8 ) THEN
     ALLOCATE(R(NPOINTS, NPOINTS))
     DO I = 1,NPOINTS
        READ(uid) R(I,:)
        ! DO J = 1, NPOINTS
        !    READ(uid) R(I,J)
        ! ENDDO
     ENDDO
  ELSE IF ( FTYPE .eq. 9 ) THEN
     ALLOCATE(R(NPOINTS, NPOINTS))
     ALLOCATE(T(NPOINTS*8))
     DO I = 1,NPOINTS
        ! there has to be a better way to do this but I can't think of it right now
        ! i've tried a few obvious things and the don't work.
        READ(uid,rec=I, iostat=stat) (T(J),J=1,NPOINTS*8)
        if (stat < 0) then
           write(*,*) ' Reached EOF at I = ', I
           EXIT
        else
           
        ENDIF
        DO J = 1,NPOINTS
           TMPI = T(((J-1)*8+1):(J*8))
           R(I,J) = tmpr
        ENDDO
     ENDDO
  ENDIF

  close(uid)
  
101 FORMAT( 4( G11.4, 2x ), ' . . . ', 4( G11.4, 2x ) )
103 FORMAT( a )
    if (nPoints .GE. 8) then
       DO I = 1,4
          IF ( ftype .eq. 4) THEN
             WRITE(*,101) E(I, 1:4), E(I, nPoints-3:nPoints)
          ELSE
             WRITE(*,101) R(I, 1:4), R(I, nPoints-3:nPoints)
          ENDIF
       ENDDO
       WRITE(*,'(a, 20x, a)') ' . . . ', '. . .'
       DO I = nPoints-3,nPoints
          IF ( ftype .eq. 4) THEN
             WRITE(*,101) E(I, 1:4), E(I, nPoints-3:nPoints)
          ELSE
             WRITE(*,101) R(I, 1:4), R(I, nPoints-3:nPoints)
          ENDIF
       ENDDO
    else
       write(fmt, "('(', i1, '( G11.4, 2x ) )')") nPoints
       DO I = 1,nPoints
          IF ( ftype .eq. 4) THEN
             WRITE(*,fmt) E(I, :)
          ELSE
             WRITE(*,fmt) R(I, :)
          ENDIF
       ENDDO
    endif


CONTAINS
  SUBROUTINE PRINT_USAGE()
    WRITE(*,*) ' printDissMatrix <matFileName> <nPoints> <realType> '
    WRITE(*,*) '        matFileName: the name of the file with the matrix data'
    WRITE(*,*) '        nPoints:     the number of points being clustered'
    WRITE(*,*) '        realType:    the size of real (real*8 or real*4) used (use 9 for MPI files)'
    STOP 1
  END SUBROUTINE PRINT_USAGE
END PROGRAM PRINTDISSMATRIX
    
