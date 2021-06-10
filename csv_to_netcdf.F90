program csv_to_netcdf

  ! version 1.0
  ! Colin Lee,    May 21, 2021

  ! Take aircraft campaign CIMS and PTR-MS data from 
  ! processed Matrix.csv files to NetCDF for processing by
  ! hierarchical clustering algorithm

  use NETCDF
  use ISO_FORTRAN_ENV

  implicit none

  integer       :: iun, ncId, variableId

  integer       :: nPts, nSpec
  real, dimension(:,:), allocatable &
                :: massSpec
  real, dimension(:), allocatable &
                :: minMS, maxMS, meanMS
  real          :: tmp
  
  character(:), allocatable :: line, name, a(:)
  character(20) :: fmt
  double precision, allocatable :: v(:)
  integer :: n, nrow, ncol, i
  integer :: status

  ! NetCDF storage variables
  integer, parameter :: NDIMS = 2
  character (len = *), parameter :: pts_NAME = "sampling points"
  character (len = *), parameter :: spec_NAME = "mass spec signal"
  integer :: pts_dimid, spec_dimid

  ! NetCDF coordinate variables
  integer, dimension(:), allocatable &
                :: ptsIdx, specIdx
  integer :: pts_varid, spec_varid
  real, parameter :: START_LAT = 25.0, START_LON = -125.0

  ! Main netCDF variable
  integer :: ms_varid
  integer :: dimids(NDIMS)

  ! It's good practice for each variable to carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: MS_UNITS = "signal"
  character (len = *), parameter :: pts_UNITS = "step"
  character (len = *), parameter :: spec_UNITS = "m/z"
  
  
  call get_command_argument(1, length=n)
  allocate(character(n) :: name)
  call get_command_argument(1, name)
  open(unit=10, file=name, action="read", form="formatted", access="stream")
  deallocate(name)

  call get_command_argument(2, length=n)
  allocate(character(n) :: name)
  call get_command_argument(2, name)
  ! Create the NetCDF file. 
  call check( nf90_create(name, nf90_clobber, ncid) )
  deallocate(name)

  ! Zip through the file and find out how many rows and columns
  ! it has
  nrow = 0
  ncol = 0
  OUTER: do while (readline(10, line))
     nrow = nrow + 1

     call split(line, a)

     if (nrow == 1) then
        ncol = size(a)
        n = len_trim(a(1))
        allocate(character(n) :: name)
        name = a(1)
     elseif (nrow == 2) then
        ncol = size(a)
        write(fmt, "('(',G0,'(G0,:,''',A,'''))')") ncol, ","
        allocate(minMS(ncol))
        allocate(maxMS(ncol))
        allocate(meanMS(ncol))
        minMS(:) = 99999.9
        maxMS(:) = -99999.9
        meanMS(:) = 0.0        
     else
        do i = 1,ncol
           read(a(i), *, iostat=status) tmp
           if (status .lt. 0) then
              print "(A, ' ', G0)", "Encountered EOF at position ", i, " in row ", nrow
           end if
           if (tmp .eq. -999999.0) then
              write(*,121) a(i), tmp, nrow, i
121           FORMAT('Found NaN (', a, ',', g11.4, ') at row ', i5, ' column ', i3)
              cycle OUTER
           end if
           if (tmp .lt. minMS(i)) minMS(i) = tmp
           if (tmp .gt. maxMS(i)) maxMS(i) = tmp
           meanMS(i) = meanMS(i) + tmp              
        end do
     end if
  end do OUTER
  do i = 1, ncol
     meanMS(i) = meanMS(i) / (nrow - 2)
  end do
  
  ! Allocate an array to save the entire dataset
  nSpec = ncol
  nPts  = nrow - 1
  ! NetCDF indexing m
  allocate(massSpec(nPts,nSpec))
  
  ! Now go back and save all the data in an array
  rewind(unit=10)
  nrow = 0
  ncol = nSpec
  do while (readline(10, line))
     if (nrow == 0) then
        ! do nothing with the first row
     else
        call split(line, a)
        do i = 1,ncol
           read(a(i), *, iostat=status) massSpec(nrow,i)
           if (status .lt. 0) then
              print "(A, ' ', G0)", "Encountered EOF at position ", i, " in row ", nrow
           end if
           if (massSpec(nrow, i) .ne. -999999.0) then
              massSpec(nrow,i) = (massSpec(nrow,i) - minMS(i)) / (maxMS(i) - minMS(i))
           end if
        end do
        if ( nrow < 11 ) then
           print "(A, ' ', G0)", "Row ", nrow, " has sum: ", sum(massSpec(nrow,:))
        end if
     end if
     ! don't increment counter until the end
     ! so we're not stuck with a bunch of
     ! nrow-1 indexing into massSpec
     nrow = nrow+1
  end do
  close(10)

  ! At this point we should have the entire dataset from the CSV file in
  ! the array massSpec. Now we just need to spit that out into a netCDF file

  ! Create the dimension data. In future, we may want to use Time as the first
  ! dimension but for now they're both just ordinal lists
  allocate(ptsIdx(nPts))
  do i = 1,nPts
     ptsIdx(i) = i
  end do
  
  allocate(specIdx(nSpec))
  do i = 1,nSpec
     specIdx(i) = i
  end do

  ! Define the dimensions.
  call check( nf90_def_dim(ncid, pts_NAME, nPts, pts_dimid) )
  call check( nf90_def_dim(ncid, spec_NAME, nSpec, spec_dimid) )

  ! Define the coordinate variables. They will hold the coordinate
  ! information, that is, the latitudes and longitudes. A varid is
  ! returned for each.
  call check( nf90_def_var(ncid, pts_NAME, NF90_INT, pts_dimid, pts_varid) )
  call check( nf90_def_var(ncid, spec_NAME, NF90_INT, spec_dimid, spec_varid) )

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  call check( nf90_put_att(ncid, pts_varid, UNITS, pts_UNITS) )
  call check( nf90_put_att(ncid, spec_varid, UNITS, spec_UNITS) )
  
  ! Define the netCDF variables. The dimids array is used to pass the
  ! dimids of the dimensions of the netCDF variables.
  dimids = (/ pts_dimid, spec_dimid /)
  call check( nf90_def_var(ncid, NAME, NF90_REAL, dimids, ms_varid) )
  
  ! Assign units attributes to the netCDF variable.
  call check( nf90_put_att(ncid, ms_varid, UNITS, MS_UNITS) )
  
  ! End define mode.
  call check( nf90_enddef(ncid) )
  
  ! Write the coordinate variable data.
  call check( nf90_put_var(ncid, pts_varid, ptsIdx) )
  call check( nf90_put_var(ncid, spec_varid, specIdx) )

  ! Write the data.
  call check( nf90_put_var(ncid, ms_varid, massSpec) )

  ! Close the file.
  call check( nf90_close(ncid) )

  deallocate(massSpec)
contains
  function readline(unit, line)
    use iso_fortran_env
    logical :: readline
    integer :: unit, ios, n
    character(:), allocatable :: line
    character(10) :: buffer

    line = ""
    readline = .false.
    do
       read(unit, "(A)", advance="no", size=n, iostat=ios) buffer
       if (ios == iostat_end) return
       readline = .true.
       line = line // buffer(1:n)
       if (ios == iostat_eor) return
    end do
  end function readline

  subroutine split(line, array, separator)
    character(*) line
    character(:), allocatable :: array(:)
    character, optional :: separator
    character :: sep
    integer :: n, m, p, i, k

    if (present(separator)) then
       sep = separator
    else
       sep = ","
    end if

    n = len(line)
    m = 0
    p = 1
    k = 1
    do i = 1, n
       if (line(i:i) == sep) then
          p = p + 1
          m = max(m, i - k)
          k = i + 1
       end if
    end do
    m = max(m, n - k + 1)

    if (allocated(array)) deallocate(array)
    allocate(character(m) :: array(p))

    p = 1
    k = 1
    do i = 1, n
       if (line(i:i) == sep) then
          array(p) = line(k:i-1)
          p = p + 1
          k = i + 1
       end if
    end do
    array(p) = line(k:n)
  end subroutine split

  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check

end program csv_to_netcdf
