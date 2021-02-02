PROGRAM	PRINTDISSMATRIX

  implicit none

  real(kind=4),allocatable,dimension(:,:) &
                     :: E
  real(kind=8),allocatable,dimension(:,:) &
                     :: R

  integer            :: i,j

  character(len=256) :: argStr, fileName
  integer            :: nside, ftype, argLen, stat
  integer            :: nPoints, uid, argCount
  logical            :: fileExist


  argCount = 1
  ! first get the filename we want to dump
  call get_command_argument( argCount, argStr, argLen, stat)
  fileName = argStr
  inquire( file=trim(fileName), exist=fileExist )
  if (stat .ne. 0 .or. .not. fileExist) &
       call PRINT_USAGE()
  open(uid, file=fileName, status='old', form='unformatted')

  argCount = argCount + 1
  ! first get the size of the side of the tile we used
  ! from the GEM-MACH background fields
  call get_command_argument( argCount, argStr, argLen, stat)
  READ( argStr(1:argLen), '(i4)' ) nside
  if (nside .le. 0 .or. stat .ne. 0) &
       call PRINT_USAGE()

  argCount = argCount + 1
  ! first get the real type (ie did the writing file use
  ! real*4 or real*8
  call get_command_argument( argCount, argStr, argLen, stat)
  READ( argStr(1:argLen), '(i4)' ) ftype

  if ( stat .ne. 0 .or. .not. (ftype .eq. 4 .or. &
       ftype .eq. 8 ) ) &
       call PRINT_USAGE()

  NPOINTS = NSIDE * NSIDE
  IF (FTYPE .eq. 4) THEN
     ALLOCATE(E(NPOINTS, NPOINTS))
     DO I = 1,NPOINTS
        DO J = 1, NPOINTS
           READ(uid) E(I,J)
        ENDDO
     ENDDO
     
  ELSE
     ALLOCATE(R(NPOINTS, NPOINTS))
     DO I = 1,NPOINTS
        DO J = 1, NPOINTS
           READ(uid) R(I,J)
        ENDDO
     ENDDO
  ENDIF

  close(uid)
  
101 FORMAT( 4( G11.4, 2x ), ' . . . ', 4( G11.4, 2x ) )
103 FORMAT( a )
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

CONTAINS
  SUBROUTINE PRINT_USAGE()
    WRITE(*,*) ' printDissMatrix <matFileName> <nSide> <realType> '
    WRITE(*,*) '        matFileName: the name of the file with the matrix data'
    WRITE(*,*) '        nSide:       the size of the tile from the GEM-MACH background fields'
    WRITE(*,*) '        realType:    the size of real (real*8 or real*4) used'
    STOP 1
  END SUBROUTINE PRINT_USAGE
END PROGRAM PRINTDISSMATRIX
    
