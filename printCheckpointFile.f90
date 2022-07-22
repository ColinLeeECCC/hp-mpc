PROGRAM	printcheckpoint

  implicit none

  real(kind=4),allocatable,dimension(:,:)  &
                     :: E
  real(kind=8),allocatable,dimension(:,:)  &
                     :: R
  integer(kind=1),allocatable,dimension(:) &
                     :: T
  logical, allocatable, dimension(:) &
                     :: live
  integer, allocatable, dimension(:,:) &
                     :: clusterPairs
  real(kind=4), allocatable, dimension(:) &
                     :: clusterDissimilarities
  
  integer(kind=1)    :: tmpi(8), tmpj(4)
  real(kind=8)       :: tmpr
  equivalence (tmpi, tmpr)
  integer            :: tmpint
  equivalence (tmpj, tmpint)
  integer            :: i,j,l
  real               :: tmpf
  equivalence (tmpj, tmpf)

  character(len=256) :: argStr, fileName, fmt
  integer            :: ftype, argLen, stat
  integer            :: nPoints, uid, argCount
  logical            :: fileExist
  integer            :: step, nPointsLive
  integer            :: dmStart
  


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

  if ( stat .ne. 0  ) &
       call PRINT_USAGE()

  open(uid, file=trim(fileName), status='old', form='unformatted', &
       access='stream') !, convert='BIG_ENDIAN')

  ! Read the step this was output at
  READ(uid, iostat=stat) tmpj(1:4)
  ! becuase the byte order is reversed between MPI writes and fortran reads
  ! we need to use this 'equivalence' step to get the real integer out
  step = tmpint
  write(*,*) ' Step = ', step

  ALLOCATE(LIVE(NPOINTS))
  READ(uid, iostat=stat) live
  write(*,*) ' Live = ', live

  ALLOCATE(clusterPairs(nPoints - 1, 2))
  DO J = 1, 2
  DO I = 1, nPoints - 1
     READ(uid, iostat=stat) tmpj(1:4)
     clusterPairs(I,J) = tmpint
  ENDDO
  ENDDO
  
  WRITE(*,*) ' Pairings: '
  DO I = 1, nPoints - 1
     WRITE(*,'(i0, ", ", i0)') clusterPairs(I,:)
  ENDDO
  

  ALLOCATE(clusterDissimilarities(nPoints - 1))
  DO I = 1, nPoints - 1
     READ(uid, iostat=stat) tmpj(1:4)
     clusterDissimilarities(I) = tmpf
  ENDDO
  
  WRITE(*,*) ' Dissimilarities: ', clusterDissimilarities

  nPointsLive = 0
  DO I = 1,nPoints
     if (live(I)) nPointsLive = nPointsLive + 1
  ENDDO

  WRITE(*,*) ' number of Live points: ', nPointsLive
  

  ALLOCATE(R(NPOINTSLIVE, NPOINTSLIVE))
  ALLOCATE(T(NPOINTSLIVE*8))
  L = 1
  DO I = 1,NPOINTS
     ! there has to be a better way to do this but I can't think of it right now
     ! i've tried a few obvious things and the don't work.
     READ(uid, iostat=stat) (T(J),J=1,NPOINTSLIVE*8)
     if (stat < 0) then
        write(*,*) ' Reached EOF at I = ', I
        EXIT
     else

     ENDIF
     if (live(I)) then
     DO J = 1,NPOINTSLIVE
        TMPI = T(((J-1)*8+1):(J*8))
        R(L,J) = tmpr
     ENDDO
     L = L + 1
     endif
  ENDDO

  close(uid)
  
101 FORMAT( 4( G11.4, 2x ), ' . . . ', 4( G11.4, 2x ) )
103 FORMAT( a )
    if (nPointsLive .GE. 8) then
       DO I = 1,4
             WRITE(*,101) R(I, 1:4), R(I, nPointsLive-3:nPointsLive)
       ENDDO
       WRITE(*,'(a, 20x, a)') ' . . . ', '. . .'
       DO I = nPointsLive-3,nPointsLive
             WRITE(*,101) R(I, 1:4), R(I, nPointsLive-3:nPointsLive)
       ENDDO
    else
       write(fmt, "('(', i1, '( G11.4, 2x ) )')") nPointsLive
       DO I = 1,nPointsLive
             WRITE(*,fmt) R(I, :)
       ENDDO
    endif


CONTAINS
  SUBROUTINE PRINT_USAGE()
    WRITE(*,*) ' printDissMatrix <matFileName> <nPoints> <realType> '
    WRITE(*,*) '        matFileName: the name of the file with the matrix data'
    WRITE(*,*) '        nPoints:     the number of points being clustered'
    STOP 1
  END SUBROUTINE PRINT_USAGE
END PROGRAM printcheckpoint
    
