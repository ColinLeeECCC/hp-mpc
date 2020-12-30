      program kzfilter_hc
      

!  version 1.0 
!  Paul Makar, September 21, 2016
! Joana Soares August 14,2017
!  Compile with s.compile -src kzfilter-GEM_MACH.ftn90 -o kzfilter.exe -librmn rmn_015.2 -O 2 on Joule
!  Compile with s.compile -src kz_hc_GEM-MACH_v18.f90 -o ../phase4/kz_hc_GEM-MACH_v18.exe -librmn rmn_016.2 -O 2 -optf='-C -traceback -init=snan -init=arrays -openmp -O2'  on Science
!  This code uses a Kolmogrov-Zurbenko moving average filter to 
!  remove high frequencies from an input (assumed hourly) data set
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! NOTE: run *.exe with 4 arguments 
!!!       nst(no.stations),mtype(metric),istd(standardization on or off),
!!!       itmp(clustering method)
!'Select the type of proximity measure:'
! 1 = Euclidean distance
! 2 = Pearson correlations (R)
! 3 = 1-R
! 4 = 1-R*EuD
!      mtype = 1
!Select standardization by range:
! 0 = no
! 1 = yes
!      istd = 0
!select proximity according to dis or dimilarity
! 1 => dissimilarity
! -1 ==> similarity
!       is = -1 ! , if mtype=2,4 the is=1
! Select clustering method
! 1-Single Linkage Method
! 2-Complete Linkage Method
! 3-Group Average Method
! 4-Weighted Average Method
! 5-Centroid Method
! 6-Median Method
! 7-Ward''s Minimum Variance Method'
! 8-Lance & William''s Beta-Flexible Method'
! 9-Belbin, Faith, & Milligan''s Beta-Flexible Method'
! itmp = 4
!!!!!!!!!!!!

!      use ARMN
!      use mpi
      use NETCDF

      implicit none

       include 'mpif.h'
       include 'netcdf.inc'
 
      integer,parameter                       :: i8 = selected_int_kind(15), &
                                                 sp = kind(1.0),             &
                                                 dp = selected_real_kind(2*precision(1.0_sp)) 
      integer                                 :: mtype,itmp,is,istd, &
                                                  itest, &
                                                  iclus1,iclId,iclIdpi, &
                                                  iclTmp,ncol, &
                                                  ir,inf,jr,l,k, &
                                                  nts,w0, &
                                                  i1,kz,mkz,pkz, & 
                                                  nsum, fu_in,fu_in1,fu_in2, fu_bg, &
                                                  key,ier, &
                                                  iIni,iEnd, &
                                                  jIni,jEnd,fcst_hr, &
                                                  iday, ithr, fractC
       integer                                 :: start_date,end_date,&
                                                  strTmp,endTmp, &
                                                  start_year, start_mon, &
                                                  start_day, end_year, &
                                                  end_mon, end_day

       integer :: pi,pi_i,pi_j
       real(dp)                                :: freq,period,avg, &
                                                  datmiss
       character(len=140)                      :: header
       character(len=40)                       :: headername
       character(len=14)                       :: cdum
       character(len=35)                       :: filein,obsheader
       integer,dimension(:),allocatable        :: id !,lk,lki,lkj 
       !real(kind=4),dimension(:),allocatable   :: lat,lon,xx,dissMatrix
       !integer(kind=4),dimension(:),allocatable:: stId1,stId2
       !integer(kind=4),dimension(:,:),allocatable :: cluster_table
       !real(kind=4),allocatable,dimension(:)   :: i_array,j_array
       real(kind=4),allocatable,dimension(:,:) :: data_st,data_stT,varbg
       !areal(kind=4),dimension(:,:),allocatable :: kzout,data_in
       character(len=70) :: obs_file, outfile
!
       real(dp)                                :: i_stack,j_stack, &
                                                  latPS,lonPS, &
                                                  sxy,sxx,syy
       
       character(len=2)                         :: fhr 
       character(len=8)                         :: dateIniStr
       character(len=8)                         :: dateEndStr
       character(len=256)                       :: bg_filename,inDir, &
                                                   outDir,mdl_file, &
                                                   outFNm,outFNmData, &
                                                   outFNmDataIn, &
                                                   outfilename
       integer :: cur_date,idate_next, next_date, dummy_hour
       character(len=8) :: model_date,cur_dateS
       integer :: model_tstep, process_step, fcst_offset
       character(len=6) :: tstep
!
! 4 FST FILES & NETCDF
       integer                     :: ncId, pVarId
       integer, dimension( NF90_MAX_VAR_DIMS) :: dimIds
       integer                     :: bggrid_npts
       integer                     :: dateo_bidon, deet_bidon, npas_bidon
       integer                     :: bggrid_ni, bggrid_nj,  bggrid_nt, &
                                      ip1_bidon, ip2_bidon, ip3_bidon, &
                                      ip1_grid, ip2_grid, ip3_grid, &
                                      ig1_grid, ig2_grid, ig3_grid, ig4_grid
       character(len=2)            :: typvar_bidon, grtyp_grid
       character(len=4)            :: nomvar_bidon,bggrid_type
       character(len=12)           :: etiket_bidon
       character(len=4)            :: bggrid_def_field, nomvar_grid
       integer                     :: bggrid_ig1, bggrid_ig2, bggrid_ig3, &
                                      bggrid_ig4, &
                                      bggrid_id, nbits_bidon, datyp_bidon
       integer                     :: ni_bidon, nj_bidon, nk_bidon, &
                                      e1_bidon, e2_bidon, e3_bidon, &
                                      swa_bidon, lng_bidon, dltf_bidon, &
                                      ubc_bidon
!
!
!  Species:
!
      integer, parameter :: nbg = 1 !#5
      integer, parameter :: nkz = 1 ! 4
      character(len=4), dimension(nbg) :: varbg_name =  &
                   (/ 'N2  '/) !# (/ 'AC  ', 'AF  ', 'O3  ', 'NO  ', 'N2  ', 'S2  '/)
      character(len=4), dimension(nbg) :: varbg_name2 =  &
                   (/ 'NO2'/) !!!'SO2 ', 'NO  ', 'NO2 ', 'O3  ', 'AF  '/)
      character(len=10), dimension(nkz) :: kzname =  &
                   (/'kz_001_001'/) !,'kz_017_003','kz_095_005','kz_523_003'/)
      !real, dimension(1) :: diss = (/0.65/)
      real, parameter :: diss = 0.65 !1-R

! In the following, a value of 1.00 is a placeholder for those species
! which do not need molecular mass to convert units
!  ppb1 = needs conversion from ug/kg to ppbv.
!  ppb2 = original output units are ppb.
!  ugm3 = needs conversion from ug/kg to ug/m3.
      character(len=4), dimension(nbg) :: unitbg  =  &
                   (/'ppb2'/) !#, 'ppb2', 'ppb2', 'ppb2', 'ugm3'/)
      real, dimension(nbg) :: mmspec_nbg = (/1.0/) !#, 1.0, 1.0, 1.0, 1.0/)
!
      real, parameter :: mmair = 28.97, rgasd = 287.00248878
!
      real(dp)             :: start,finish,start_taskS,end_taskS, &
                              start_taskU,end_taskU 
      character(len=150) :: timeFNm

      logical, parameter :: ascii_out = .false.
      logical, parameter :: centroids = .false.
      character(len=1)   :: dataAvail,distMavail,clustersAvai,cntroidsAvail
      logical            :: distMatrix_not_available, &
                            clusters_not_available,centroids_not_available, &
                            writeOutR, writeOutC
       
      integer :: narg,stopt
      character(len=10) :: nstS,stOptS,lonS,latS
      character(len=1) :: metricS,istdS,itmpS
     
      integer(kind=4),dimension(8) :: date2julianS,date2julianE
      real                         :: jDayS,jDayE,jDaySTmp,jDayETmp
      integer                      :: err
      character(len=8)             :: start_dateS,start_dateE,metric

       
      !mpi variables
       integer                             :: ierr,myid,myprocs,master, &
                                              llocId,piId,st1Id,st2Id
       integer :: ni,nj
       integer                             :: i,j,il,jl,ndays,&
                                              ict,icl,icl2, &
                                              minp,maxp,minpst1,maxpst1, &
                                              minT,maxT,prev_i,prev_j, &
                                              st1Tmp, st2Tmp, indTmp, idTmp
       integer                               :: mpiTypei8   ! for MPI kind def
       integer,allocatable,dimension(:)      :: lki,lkj,lk,req
       integer(i8)                           :: nst, nst2, npairs, iclus2, &
                                                igc, igc2, &
                                                inc, ltmp, &
                                                strSi, strSe, &    ! check size of pairs 
                                                strSid, endSid, npairsId,&
                                                jpl1, jpl2,lglobId,ist,&
                                                check, checkL, checkU
       integer,allocatable,dimension(:)      :: station1,station2, & !ij
                                                rbuf_1,rbuf_2, id4cl
       integer(i8),allocatable,dimension(:)  :: lloc, strS, endS, sti, stj
       !integer,allocatable,dimension(:,:):: stations
       real                                  :: isr,ntr
       real(i8)                              :: varTmpi8, dissM_i, dissM_j
       real(dp)                              :: varTmp,sum_st1st2, &
                                                data_st1,data_st2, &
                                                valTmp, minDId, &
                                                dissMetric_tmp ! perR
       real(dp),allocatable,dimension(:)     :: sum_st,sum_stSQR, &
                                                dissMetric,minDloc, &
                                                xx, rbuf_i, rbuf_j
       logical                               :: chk_prog
       logical,allocatable,dimension(:)      :: w,wj
       integer status(MPI_STATUS_SIZE)
       integer, parameter :: root = 0 
      
       !integer statusReq(MPI_STATUS_SIZE,2)

       !type minD
       !  integer       :: ind, rank, st1, st2
       !  real(kind=4)  :: val
       !end type minD
       !type(minD) :: minDin, minDout
       !real(kind=4) :: minDin, minDout
       !real(kind=4),dimension(2) :: minDin, minDout ! val, id, st1, st2, rank
       real(dp),dimension(2) :: minDin, minDout ! val, id, st1, st2, rank
!  Model surface ip1 values.
!
      !integer, parameter :: ip1g_surface = 93423264
      !integer, parameter :: ip1t_surface = 76696048
      !integer, parameter :: ip1m_surface = 75597472
      integer :: ip1g_surface,ip1t_surface,ip1m_surface
!
! parallel
      character(2)            :: ithrS 
      integer                 :: idThr,numThr
      integer,allocatable, &
              dimension(:)    :: strT,endT 
      real,allocatable, &
              dimension(:)    :: strJ,endJ 
      integer, external       :: omp_get_max_threads, &
                                 omp_get_thread_num, &
                                 omp_get_num_threads
!
! FST functions
      integer,external ::  fstinf, fstprm, fstlir, ezqkdef, gdrls, newdate, &
                           ezsetopt, fnom, fstouv, fstfrm, fclos, gdxyfll, &
                           gdllsval, gdllfxy, fstopc

! 
! MPI
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, myprocs, ierr)
      if (myid == root) write(*,*) 'number of processors', myprocs
!    
!  general settings
! 
      ip1g_surface = 93423264
      ip1t_surface = 76696048
      ip1m_surface = 75597472
!
! Read in the passed argument list(now in a file):
       fu_in = 15
      open(unit=fu_in,file='arglist.lst',status="old",iostat=ierr)
      if (ierr /= 0) then
       write(6,*)'Abnormal...error in opening file'
       stop
      endif
      read(fu_in,'(a256)') inDir                                           ! where GEM-MACH data is located
      read(fu_in,'(i4,2(1x,i2.2))') start_year, start_mon, start_day
      read(fu_in,'(i4,2(1x,i2.2))') end_year, end_mon, end_day
      read(fu_in,*) fractC                                          ! fraction of x-y no cells
      nst=fractC
      read(fu_in,'(a256)') outDir
      read(fu_in,*) distMatrix_not_available
                !amtype,itmp, fracti, &
                !dataAvail,distMavail,clustersAvai,cntroidsAvail
      read(fu_in,*) writeOutR
      read(fu_in,*) clusters_not_available
      read(fu_in,*) writeOutC
      close(fu_in)
!
      if (myid == root) then
      write(*,*) 'input'
      write(*,*) inDir
      write(*,*) start_year, start_mon, start_day
      write(*,*) end_year, end_mon, end_day
      write(*,*) fractC
      write(*,*) outDir
      write(*,*) distMatrix_not_available
      end if

      ! hard coded for now
      mtype = 3       
      istd = 0
      itmp = 1
      is = -1
!
!  begining and end of the run, and time steps requested   
      fcst_hr = 18
      date2julianS = (/start_year,start_mon,start_day,fcst_hr,00,00,00,00/)
      date2julianE = (/end_year,end_mon,end_day,fcst_hr,00,00,00,00/)
      call d2j(date2julianS,jDayS,err)
      call d2j(date2julianE,jDayE,err)
      nts = int((jDayE-jDayS+1)*24)
      ndays = int((jDayE-jDayS+1))
      start_date = date2julianS(1)*10000+date2julianS(2)*100+date2julianS(3)
      end_date = date2julianE(1)*10000+date2julianE(2)*100+date2julianE(3)
      write (dateIniStr,'(i8)') start_date ! converting integer to string using a 'internal file'
      write (dateEndStr,'(i8)') end_date ! converting integer to string using a 'internal file'
      write(fhr,'(i2.2)') fcst_hr
!
! ################  starting ############################ 
!
      start_taskS = MPI_Wtime()
      if (myid == root) then
      write(6,*) '!-------------------------------------!'
      write(6,*) 'Hi world! Network analysis in progress.'
      write(6,*) '!-------------------------------------!'
      end if
      if (mtype == 1) then
        metric = 'EuD'
      elseif (mtype == 2) then
        metric = 'EuDN'
      elseif (mtype == 3) then
        metric = 'R'
      elseif (mtype == 4) then
        metric = '1-R'
      elseif (mtype == 5) then
        metric = '1-RxEuD'
      elseif (mtype == 6) then
        metric = 'variance'
      endif
!
      do l = 1,nbg
        headername = '_modelPhase4_v24_'
        outFNmData = trim(outDir)//trim(varbg_name2(l))//&
                 trim(headername)
        outFNm = trim(outDir)//trim(varbg_name2(l))//&
                 trim(headername)//trim(metric)//'_'
        if (myid == root) then
          write(6,*) 'start_date: ',dateIniStr,', end_date: ',dateEndStr
          !write(6,'(a25,1x,a5,1x,a6)') ' central point (lat,lon):',latS,lonS
          write(6,*) 'species: ',trim(varbg_name(l))
          write(6,*) 'time series length:',nts, 'no. stations:', nst
          write(6,*) 'metric(1-5),standardization(0-no,1-yes):',mtype,istd
          write(6,*) 'clustering method (1-9):,',itmp
          !write(6,*) 'no. clusters for optimization:',stOpt
          write(6,*) 'output filename: ',outFNm
        end if
        bggrid_def_field = trim(varbg_name(l))
        write(*,*) 'field requested ',bggrid_def_field        
!
! Get information about the model data
        readNETCDF: if (.true.) then
! 
          fu_bg = 15
          write(model_date, '(i8.8)') start_date
          write(tstep, '(i6.6)') 60
          bg_filename = model_date//trim(fhr)//'_'//trim(tstep)//'p'
          bg_filename = trim(inDir)//model_date//trim(fhr)//&
                        '_'//trim(tstep)//'p.netcdf4.compressed'
!   Reading grid parameters
          write(*,*) 'bg_filename ', bg_filename
          ierr = nf90_open( trim(bg_filename), NF90_NOWRITE, ncId) !, chunksize)
          if (ierr < 0) then
            print *,'Ending,error while opening ',trim(bg_filename)
            stop
          end if
          !write( *, '(a,x,a)' ) 'Reading variable from', trim(bg_filename)
          ierr = nf90_inq_varid( ncId, bggrid_def_field, pVarId )
          ierr = nf90_inquire_variable( ncId, pVarId, dimIds = dimIds )
          ierr = nf90_inquire_dimension( ncID, dimIds(1), len = bggrid_ni )
          ierr = nf90_inquire_dimension( ncID, dimIds(2), len = bggrid_nj )
          ierr = nf90_inquire_dimension( ncID, dimIds(3), len = bggrid_nt )
          write(*,*) 'dims', bggrid_ni, bggrid_nj, bggrid_nt
          !allocate(varbg(bggrid_ni, bggrid_nj))
          !ierr = nf90_get_var( ncId, pVarId, varbg ) !

          ierr = nf90_close( ncId )
!
        end if  readNETCDF

        readFST: if (.false.) then 
!
! Get information about the model data
!       ### turn off FST information messages
        ier = fstopc('MSGLVL', 'SYSTEM', 0)
        fu_bg = 15 
        write(model_date, '(i8.8)') start_date
        write(tstep, '(i6.6)') 60
        bg_filename = model_date//trim(fhr)//'_'//trim(tstep)//'p'
        bg_filename = trim(inDir)//model_date//trim(fhr)//&
                      '/'//trim(bg_filename)
        ier = fnom(fu_bg, trim(bg_filename), 'RND+OLD+R/O', 0)
        if (ier < 0) then
          print *,'Ending,error while opening ',trim(bg_filename)
          stop
        end if
        ier = fstouv(fu_bg, 'RND')
!   Reading grid parameters
        key = fstinf(fu_bg,bggrid_ni,bggrid_nj, nk_bidon, -1, ' ',&
                     -1, -1, -1, ' ', bggrid_def_field)
        if (key < 0) then
          write(6,*) 'Key value negative after opening file',&
                     ' with grid definitions'
          stop
        end if

        ier = fstprm(key, dateo_bidon, deet_bidon, npas_bidon,         &
                     bggrid_ni, bggrid_nj, nk_bidon, nbits_bidon,      &
                     datyp_bidon, ip1_bidon, ip2_bidon, ip3_bidon,     &
                     typvar_bidon,nomvar_bidon,etiket_bidon,bggrid_type, &
                     bggrid_ig1, bggrid_ig2, bggrid_ig3, bggrid_ig4,   &
                     swa_bidon, lng_bidon, dltf_bidon, ubc_bidon,      &
                     e1_bidon, e2_bidon, e3_bidon)
        if (ier < 0) then
          write(6,*) 'ier value negative after attempting to ', &
                     'read grid parameters'
          stop
        end if
!    Initialize the bg grid for ezscint functions
        ier = ezsetopt('VERBOSE', 'NO')
        bggrid_id = ezqkdef(bggrid_ni,bggrid_nj,bggrid_type,bggrid_ig1,&
                            bggrid_ig2, bggrid_ig3, bggrid_ig4, fu_bg)
        ier = gdxyfll(bggrid_id,i_stack,j_stack,latPS,lonPS + 360.0,1) !!! REMOVE
        ier = fstfrm(fu_bg)
        ier = fclos(fu_bg)
!   
       end if readFST
!
!  Get grid and number of clusters
       !nst = bggrid_ni * bggrid_nj
       ni = bggrid_ni !int(bggrid_ni/10)
       nj = bggrid_nj !int(bggrid_nj/10)
       i =  fractC
       ni = int(bggrid_ni/i)
       nj = int(bggrid_nj/i)
!
!
! FOR TESTING with seq
       ni =  fractC
       nj =  fractC
!
!
       nst = ni*nj
!  lower and upper bound of the clusering
       iclus1 = 2
       iclus2 = nst-1
       npairs = map(nst,nst-1,nst)
       if (myid == root) then
         write(*,*) 'ni nj nst npairs', ni, nj, nst, npairs 
       end if
       !if (myid ==0 ) then
       !  write(*,*) 'range: pairs and strS',range(npairs),range(strS) 
       !  !stop
       !end if 
       !if (npairs > 10**range(npairs)) then
       !  write(6,*) 'Overflow: either reduce number of grid cells',&
       !             ' or increase kind'
       !  !stop
       !end if 

       call MPI_Type_create_f90_integer(15, mpiTypei8, ierr)
       allocate(strS(myprocs),endS(myprocs))
       allocate(id4cl(myprocs))
!
!  Set the no of station, pairs and clusters per PI in the root, and
!  and brodcast this information to the other PIs  
!
       if (myid == root) then 
!
!  number of pairs (size) allocated to each pi
         inc = (npairs/(myprocs)*1_i8) 
         strS = 1
         endS = npairs
         do idThr = 1,myprocs-1
           endS(idThr) = strS(idThr) + inc - 1_i8
           strS(idThr+1) = strS(idThr) + inc
         end do
         write(*,*) 'strS',strS
         write(*,*) 'endS',endS
!
! calculate the number of cluster information attributed to each PI, 
         inc = floor(real((nst-1)/(myprocs)*1_i8))
         do idThr = 1,myprocs-1
           id4cl(idThr) = inc * idThr
         end do
         id4cl(myprocs) = nst-1

       end if ! myid = 0
!
!    Brodcast information to the diferent processors
       call MPI_Bcast(strS, myprocs, mpiTypei8, root, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(endS, myprocs, mpiTypei8, root, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(id4cl, myprocs, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
        
!
!  Reading data if the GEM-MACH data is not available in a binary format 
!
       IF_DATA_NAVAIL: if (distMatrix_not_available) then 
!        
         npairsId = (endS(myid+1)-strS(myid+1)) + 1_i8
         write(*,'(a2,1x,i2.2,a36,i5.5)')&
                  'PI',myid,'computing distance metric for npairs ', &
                   npairsId
         write(*,*) 'PI',myid, &
                    'start & end pair', strS(myid+1),endS(myid+1) 
         allocate(station1(npairsId),station2(npairsId))  
         ict = 0
         !!$OMP PARALLEL DO PRIVATE(i,j,igc)
         do i = 1,nst-1
           do j = i+1,nst
             igc = j + (i-1_i8)*nst - (i*(i+1_i8))/2_i8 
             if (igc >= strS(myid+1) .and. igc <=endS(myid+1)) then
               ict = ict + 1
               station1(ict) = i
               station2(ict) = j
             end if
           end do
         end do
         !!$OMP END PARALLEL DO
!         write(*,*) 'station', station1
!         write(*,*) 'station2', station2
         minp = minval(station1)
         maxp = maxval(station2)
         nst2 = maxp-minp+1
         
!
!      store grid indexes for each station
         allocate(sti(nst2),stj(nst2))
         ist = 0
         igc = 0 !minp-1
         if  (writeOutR) then
             fu_in = 30+myid
             write(ithrS,'(i2.2)') myid
             write(nstS,'(i5.5)') nst
             outFNmDataIn = trim(outFNmData)//'check_sti_stj_'//&
                            trim(ithrS)//'_'//trim(nstS)//'.dat'
             open(unit=fu_in,file=trim(outFNmDataIn),action='write', &
                  status='replace', form='formatted', iostat=ierr)
             write(fu_in,*) 'igc,ist,sti,stj'
         end if

         do i = 1,ni
           do j = 1,nj
             igc = igc + 1
             if (igc >= minp .and. igc <= maxp) then
               ist = ist + 1
               sti(ist) = i
               stj(ist) = j
               if (writeOutR) then
                  write(fu_in,*) igc,ist,sti(ist),stj(ist)
               endif
             end if
           end do
         end do

         if  (writeOutR) close(fu_in)

         allocate(data_st(nst2,ndays*24),sum_st(nst2),sum_stSQR(nst2))
         sum_st = 0.
         sum_stSQR = 0.
         cur_date = start_date
         !allocate(varbg(bggrid_ni,bggrid_nj))
         allocate(varbg(ni,nj))
         !write(*,*) &
         ! 'ict igc2 j i varbg(i,j) sum_st(igc2) sum_stSQR(igc2)'
!
         if  (writeOutR) then
             fu_in = 30+myid
             write(ithrS,'(i2.2)') myid
             write(nstS,'(i5.5)') nst
             outFNmDataIn = trim(outFNmData)//'checkSums_data_sum_'//&
                            trim(ithrS)//'_'//trim(nstS)//'.dat'
             open(unit=fu_in,file=trim(outFNmDataIn),action='write', &
                  status='replace', form='formatted', iostat=ierr)
             write(fu_in,*) 'ict,ist,i,j,val,sum_st,sum_stSQR'
         end if
 
        !# reading files for the time period requested
         write(*,*) 'reading files'
         do iday = 1,ndays
           call j2d(jDayS+(iday)*1.0,cur_date,err)
           write(model_date, '(i8.8)') cur_date
          !!$OMP PARALLEL DO SHARED(iday,inDir,model_date,fhr,varbg_name,)REDICTION sum_st,sum_stSQR
           do process_step = 1, 24
             model_tstep = process_step
             ict = 24*(iday-1)+process_step
             write(tstep, '(i6.6)') model_tstep * 60

             bg_filename = model_date//trim(fhr)//'_'//trim(tstep)//'p'
             bg_filename = trim(inDir)//model_date//trim(fhr)//&
                           '_'//trim(tstep)//'p.netcdf4.compressed'
             fu_bg = 10+process_step !+ithrd

!   Open input files and read in the requested field for this model timestep
             readNETCDF2: if (.true.) then
!
               ierr = nf90_open( trim(bg_filename), NF90_NOWRITE, ncId) !, chunksize)
               if (ierr < 0) then
                 print *,'Ending,error while opening ',trim(bg_filename)
                 stop
               end if
               !write( *, '(a,x,a)' ) 'Reading variable from', trim(bg_filename)
               ierr = nf90_inq_varid( ncId, bggrid_def_field, pVarId )
               ierr = nf90_inquire_variable( ncId, pVarId, dimIds = dimIds )
               !ierr = nf90_inquire_dimension( ncID, dimIds(1), len = bggrid_ni )
               !ierr = nf90_inquire_dimension( ncID, dimIds(2), len = bggrid_nj )
               !ierr = nf90_inquire_dimension( ncID, dimIds(3), len = bggrid_nt )
               ierr = nf90_get_var( ncId, pVarId, varbg )
!               write(6,*) 'Input values of varbg parallel code: '
!               write(*,2011) varbg
!2011           format(8(1x,1pe14.7,1x))
               !write(*,*) 'print varbg(1:20,1:20) ', varbg(1:20,1:20)
               ierr = nf90_close( ncId )
!
              end if readNetcdf2 
!
             readFST2: if (.false.) then
!
               ier = fnom(fu_bg, trim(bg_filename), 'RND+OLD+R/O', 0)
               if (ier < 0) then
                 print *,'Ending,error while opening ',trim(bg_filename)
                 stop
               end if
               ier = fstouv(fu_bg, 'RND')
               ier = fstlir(varbg(:,:), fu_bg, ni_bidon, nj_bidon, nk_bidon, &
                             -1, '', ip1t_surface, -1, -1, '', varbg_name)
               ier = fstfrm(fu_bg)
               ier = fclos(fu_bg)
! 
             end if readFST2
!
!!   Start to read data for grid cells and compute the distance matrix 
!
             do ist = 1,nst2
               i = sti(ist)
               j = stj(ist)
               varTmp = varbg(i,j)
               data_st(ist,ict) = varTmp
               sum_st(ist) = sum_st(ist)+varTmp
               sum_stSQR(ist) = sum_stSQR(ist)+varTmp*varTmp
               if  (writeOutR) write(fu_in,*) ict,ist,i,j,varTmp, &
                                              sum_st(ist),sum_stSQR(ist)
             end do
!
           end do !# process_step, read hourly files
          !!$OMP END PARALLEL DO
         end do !# cycle over the period requested (1:ndays)
!
         if  (writeOutR) close(fu_in)
!
         deallocate(varbg)
         deallocate(sti,stj)
!!
         write(*,*) 'writting ended'
         minpst1 = minp
         maxpst1 = maxval(station1)
         allocate(dissMetric(npairsId))
         dissMetric = 0.
! Compute intermediate sums for computing metric below
         ntr = 1/real(ndays*24)
         isr = real(is)
         !!$OMP PARALLEL DO SHARED(npairsId,station1,station2,minpst1,nts,data_st,sum_st,sum_stSQR,ntr,isr,dissMetric)
         write(*,*) 'writeOutR = ', writeOutR
         if  (writeOutR) then
            fu_in = 30+myid
            write(ithrS,'(i2.2)') myid
            write(nstS,'(i5.5)') nst
            outFNmDataIn = trim(outFNmData)//'checkSumsFinal_'//&
                           trim(ithrS)//'_'//trim(nstS)//'.dat'
            open(unit=fu_in,file=trim(outFNmDataIn),action='write', &
                 status='replace', form='formatted', iostat=ierr)
         write(fu_in,*) 'sum_st: ', sum_st
         write(fu_in,*) 'sum_stSQR: ', sum_stSQR
         write(fu_in,*) 'i,station1(i),station2(i),data_st(st1),', &
                        'data_st(st2),st1Tmp,st2Tmp,' , &
                        'data_st(st1Tmp,ict), data_st(st2Tmp,ict)'
         write(fu_in,*) 'i,st1,st2,st1,st2,sum_st1st2', &
                           ' SXY, SXX, SYY,', &
                           'dissMetric(i)'
         end if
         write(*,*) 'compute metric'
         write(*,*) 'fu_in = ', fu_in
         do i = 1,npairsId
            st1Tmp = station1(i)-station1(1)+1
            st2Tmp = station2(i)-minpst1+1
            sum_st1st2 = 0.
            do ict = 1,nts
              sum_st1st2 = sum_st1st2 + &
                           data_st(st1Tmp,ict) * data_st(st2Tmp,ict)
            end do
!            write(*,*) 'st1Tmp,st2Tmp,nst',st1Tmp,st2Tmp,nst
            if (writeOutR) then
               write(fu_in,*) i,station1(i), station2(i), &
                           st1Tmp,st2Tmp, data_st(st1Tmp,nts), &
                           data_st(st2Tmp,nts)
            endif
            data_st1 = sum_st(st1Tmp)
            data_st2 = sum_st(st2Tmp)
            SXY = sum_st1st2 - (data_st1*data_st2*ntr)
            SXX = sum_stSQR(st1Tmp)-(data_st1*data_st1*ntr)
            SYY = sum_stSQR(st2Tmp)-(data_st2*data_st2*ntr)
            dissMetric(i) = (1.0-SXY/SQRT(SXX*SYY))*isr
            if  (writeOutR) then
              write(fu_in,*) i, st1Tmp, st2Tmp, &
                             sum_st1st2, SXY, SXX, SYY, &
                             dissMetric(i)
            end if
         end do
          write(*,*) 'end compute metric'
          if  (writeOutR) close(fu_in)
         !!$OMP END PARALLEL DO
          deallocate(data_st, sum_st, sum_stSQR)
!           
          if  (writeOutR) then
!     ### save data
           fu_in = 20+myid
           write(ithrS,'(i2.2)') myid
           write(nstS,'(i5.5)') nst
           outFNmDataIn = trim(outFNmData)//'OSdomain_pairs_metric_'//&
                          trim(ithrS)//'_'//trim(nstS)//'.dat'

           open(unit=fu_in,file=trim(outFNmDataIn),action='write', &
                status='replace', form='formatted', iostat=ierr)
           if (ierr /= 0) then
             write(*,*)'Abnormal...error in opening output file'
             stop
           endif

           do i = 1,npairsId
             write(fu_in,*) station1(i),station2(i),dissMetric(i)
           end do
           close(fu_in)
!
          end if
!
        end if  IF_DATA_NAVAIL ! data is not avaialable 
!
        end_taskS = MPI_Wtime()
        write(*,'(a2,1x,i2.2,1x,a40,f10.2)') 'PI',myid, &
                          'compute distance matrix running time (s)' , &
                          end_taskS-start_taskS
!
!
! Start clustering
!
        start_taskS = MPI_Wtime()
        IF_CLUSTERS_NAVAIL: if (clusters_not_available) then !! needs input from sub compute_metric
!
         write(*,'(a2,1x,i2.2,a30)') 'PI',myid,'clustering stations'
         !!$OMP PARALLEL SHARED(numThr)
         ! numThr = omp_get_num_threads()
         !!$OMP END PARALLEL
          numThr = 1
!         
          !if (myid /= 0 .and. distMatrix_not_available == .false.) then
          !if (myid /= 0 .and. distMatrix_not_available .eqv. .false.) then
          ! write(ithrS,'(i2.2)') myid
          ! outFNmDataIn = trim(outFNmData)//'OSdomain_pairs_metric_'//&
          !                ithrS//'.bin'
          ! open(unit=31,file=trim(outFNmDataIn),form='unformatted', &
          !                                      action='write')
          ! do i = 1,npairsId
          !   write(31) station1(i),station2(i),dissMetric(i)
          ! end do
          ! close(31)
          !end if
          !if (distMatrix_not_available == .false.) &
          if (distMatrix_not_available .eqv. .false.) &
                      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!

            fu_in = 20+myid
            write(ithrS,'(i2.2)') myid
            write(nstS,'(i10.10)') nst
            outFNmDataIn = trim(outFNmData)//&
                           'check_parallelcheckClusters_'//&
                           trim(ithrS)//'_'//trim(nstS)//'.dat'
            open(unit=fu_in,file=trim(outFNmDataIn),action='write', &
                 status='replace', form='formatted', iostat=ierr)
            if (ierr < 0) then
              print *,'Ending,error while creating ',trim(outFNmDataIn)
              stop
            end if
            write(fu_in,*) ' xx lki lkj st2Id wj(st1Id)'

          if (writeOutC) then
            fu_in1 = 40+myid
            write(ithrS,'(i2.2)') myid
            write(nstS,'(i5.5)') nst
            outFNmDataIn = trim(outFNmData)//'check_clustering_Bcast_'//&
                           trim(ithrS)//'_'//trim(nstS)//'.dat'
            open(unit=fu_in1,file=trim(outFNmDataIn),action='write', &
                 status='replace', form='formatted', iostat=ierr)
            if (ierr < 0) then
              print *,'Ending,error while creating ',trim(outFNmDataIn)
              stop
            end if         
            write(fu_in1,*) 'clustering'

            fu_in2 = 50+myid
            write(ithrS,'(i2.2)') myid
            write(nstS,'(i5.5)') nst
            outFNmDataIn = trim(outFNmData)//'checkClusters_'//&
                           trim(ithrS)//'_'//trim(nstS)//'.dat'
            open(unit=fu_in2,file=trim(outFNmDataIn),action='write', &
                 status='replace', form='formatted', iostat=ierr)
            if (ierr < 0) then
              print *,'Ending,error while creating ',trim(outFNmDataIn)
              stop
            end if
            write(fu_in2,*) '        icl,         i        st1(i)',&
               '    st2 w wj(st1),wj(st2) minDloc(idThr),dissMetric(i)'

          end if
!
          !call percentiler(dissMetric,npairs,0.9,perR)
!
! calculate the number of stations attributed to each PI, 
! note that we need nclusters-1, not nst.
          inc = floor(real((nst-1)/(myprocs)*1_i8))
          if (myid == myprocs-1 .and. mod(nst-1,(myprocs)*1_i8) /= 0.) then
            inc = nst-1 - id4cl(myprocs-1)
          end if
! Track pairs clustered (lki,lkj) and at what level (xx)
! arrays are distributed evenly between PIs, inc=nst/npis
          allocate(lki(inc),lkj(inc),xx(inc))
          lki = 0 
          lkj = 0
          xx = 0 
!
! Track the pairs being clustered (w), same length as no of pairs distributed per pi
          allocate(w(npairsId)) !,wj1(npairsId),wj2(npairsId))
          w = .false.          
! Track individual stations being clustered (wj)
          allocate(wj(nst))
          wj = .false.
!a
! Local 
          allocate(minDloc(numThr))
          allocate(lloc(numThr))
          allocate(rbuf_i(myprocs), rbuf_j(myprocs))
          allocate(rbuf_1(myprocs), rbuf_2(myprocs))
!
! Let clustering begin 
          iclId = 0 !counter to distribute the metric and stations clustered will be located (distributed across the PIs)
          iclIdpi = 0 !counter for individual processors
          do icl = 1,nst-1
            write(*,*) 'cluster step: ', icl, '########################'
            lloc = 0
            minDloc = 999999999.
            lloc = -999999999
           !!$OMP PARALLEL PRIVATE(i,idThr)
           ! idThr = omp_get_thread_num()+1
            idThr = 1
            !write(fu_in2,*) 'B: dissMetric',dissMetric
            !write(fu_in2,*) 'w',w
            !write(fu_in2,*) 'wj',wj
            do i = 1,npairsId
              if (writeOutC) write(fu_in2,*) icl,i,station1(i),station2(i), &
                             w(i),wj(station1(i)),wj(station2(i)),&
                             minDloc(idThr),dissMetric(i)
              if (w(i)) cycle                                 ! excludes pairs that have been clustered
              !if (dissMetric(i) >= perR) cycle
              if (wj(station1(i)) .or. wj(station2(i))) cycle ! excludes "stations" that have been clustered
              if (dissMetric(i) > minDloc(idThr)) cycle
              minDloc(idThr) = dissMetric(i)
              lloc(idThr) = i

              write(fu_in2,*) icl,i,station1(i),station2(i), &
                              w(i),wj(station1(i)),wj(station2(i)),&
                              minDloc(idThr),dissMetric(i)
            end do
            !!$OMP END PARALLEL
!
            minDId = minval(minDloc) ! val, id, st1, st2, rank
            if (writeOutC) write(fu_in2,*) &
                            'minDloc(idThr),lloc(idThr),minDId', &
                             minDloc(idThr),lloc(idThr),minDId
            if (minDId == 999999999.) then
              minDin(1) = minDId
              minDin(2) = myid
            else
              do idThr = 1,numThr
                if(minDloc(idThr) == minDId) then 
                  llocId = lloc(idThr)
                  st1Id = station1(lloc(idThr))
                  st2Id = station2(lloc(idThr))
                end if
              end do
            end if
!
!   Get the global minimum: level, stations, pi
            if (writeOutC) write(fu_in1,*) 'before Allreduce: ', &
                                           'st1Id,st2Id,llocId,minDId',&
                                            st1Id,st2Id,llocId,minDId
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
            minDin(1) = minDId
            minDin(2) = myid
            call MPI_Allreduce(minDin, minDout, 1, &
                               MPI_2DOUBLE_PRECISION,MPI_MINLOC, & ! 0, &
                               MPI_COMM_WORLD, ierr)
            valTmp = minDout(1)
            piId = minDout(2)
!
            minDin(1) = minDId
            minDin(2) = st1Id
            call MPI_Allreduce(minDin, minDout, 1, &
                               MPI_2DOUBLE_PRECISION,MPI_MINLOC, & ! 0, &
                               MPI_COMM_WORLD, ierr)
            st1Id = minDout(2)
!
            minDin(1) = minDId
            minDin(2) = st2Id
            call MPI_Allreduce(minDin, minDout, 1, &
                               MPI_2DOUBLE_PRECISION,MPI_MINLOC, & ! 0, &
                               MPI_COMM_WORLD, ierr)
            st2Id = minDout(2)
! 
            if (writeOutC) write(fu_in1,*) 'after Allreduce: ', &
                                           'st1Id,st2Id,llocId,minDId',&
                                            st1Id,st2Id,llocId,valTmp
!
! flag the pair that has been clustered
            if (myid == piId) w(llocId) = .true.  
!
4444  format(a77,2(1x,i2.2),1x,f8.5,1x,i2.2)
!
! check to which PI the information should go to
            iclId = iclId+1
            if (myid == 0) then 
              iclTmp = 1
            else
!              write(*,*) 'id4cl,myid',id4cl,myid
              iclTmp = id4cl(myid)+1
!               write(*,*) 'myid,id4cl(myid),id4cl(myid+1),iclTmp', &
!                           myid,id4cl(myid),id4cl(myid+1),iclTmp
            end if
!            write(*,*) 'myid,iclId,iclTmp,id4cl: ', &
!                        myid,iclId,iclTmp,id4cl(myid+1)
            if (iclId >= iclTmp .and. iclId <= id4cl(myid+1)) then
              iclIdpi = iclIdpi+1
              xx(iclIdpi) = valTmp
              lki(iclIdpi) = st2Id
              lkj(iclIdpi) = st1Id
              if (writeOutC) write(fu_in2,*) icl,lki(iclIdpi),lkj(iclIdpi), &
                                             wj(st2Id),xx(iclIdpi) 
            end if
!
!  flag station2 of the clustered pair
!            write(*,*) 'myid,icl,st2Id',myid,icl,st2Id
            wj(st2Id) = .true.
!
!       Recompute the metric based on the pairs of stations 
!       that clustered last with all the other stations available
            !!$OMP DO PRIVATE(igc,jpl1,jpl2,i,j)
!
            do igc = 1,nst
              call MPI_Barrier(MPI_COMM_WORLD, ierr)
              if (writeOutC) write(fu_in2,*)  &
                             'myid,st1Id,st2Id,igc,wj(igc)', &
                              myid,st1Id,st2Id,igc,wj(igc), &
                             '################'
              if (wj(igc)) cycle
              if (igc == st1Id .or. igc == st2Id) cycle
!
!   get station for a possible pair based on station st2Id
              jpl1 = min(st2Id,igc)
              jpl2 = max(st2Id,igc)
              i = map(nst,jpl1,jpl2) ! i< j  i = jpl2 + (jpl1-1)*nst - (jpl1*(jpl1+1))/2
              if (i >= strS(myid+1) .and. i <= endS(myid+1) &
                                    .and. i <= npairsId) then 
                il = i
                dissM_i = dissMetric(il)
!           write(fu_in2,*) '1i is in the current PI,myid, i,il, dissM_i',&
!                           myid,i,il, dissM_i
                pi_i = myid
              elseif (i >= strS(myid+1) .and. i <= endS(myid+1) &
                                        .and. i > npairsId) then
                il = i - strS(myid+1)+ 1 
                dissM_i = dissMetric(il)
!           write(fu_in2,*) '2i is in the current PI,myid, i,il, dissM_i',&
!                           myid,i,il, dissM_i
                pi_i = myid
 
              else 
                do pi = 1, myprocs
                  if (i >= strS(pi) .and. i <= endS(pi)) then
                    il = i - strS(pi)+ 1 
                    pi_i = pi - 1
                    exit
                  end if
                end do 
                dissM_i = dissMetric(il)
!         write(fu_in2,*) '3i is not in the current PI,myid,i,il,dissM_i',&
!                           myid,i,il, dissM_i 
              end if  
              call MPI_Gather(dissM_i, 1, MPI_DOUBLE_PRECISION, &
                               rbuf_i, 1, MPI_DOUBLE_PRECISION,0, &
                               MPI_COMM_WORLD, ierr)
!
!   get station for a possible pair based on station st1Id
              jpl1 = min(st1Id,igc)
              jpl2 = max(st1Id,igc)
              j = map(nst,jpl1,jpl2) ! j = jpl2 + (jpl1-1)*nst - (jpl1*(jpl1+1))/2
              if (j >= strS(myid+1) .and. j <= endS(myid+1) &
                                    .and. j <= npairsId) then 
                jl = j
                dissM_j = dissMetric(jl)
!                write(fu_in2,*) '1j is in the current PI, ', &
!                                'myid, j,jl, dissM_j',&
!                                 myid,j,jl, dissM_j
                pi_j = myid
              elseif (j >= strS(myid+1) .and. j <= endS(myid+1) &
                                        .and. j > npairsId) then
                jl = j - strS(myid+1)+ 1
                dissM_j = dissMetric(jl)
!                write(fu_in2,*) '2j is in the current PI, ', &
!                                'myid, j,jl, dissM_j',&
!                                 myid,j,jl, dissM_j
                pi_j = myid
              else 
                do pi = 1, myprocs
!                  write(*,*) 'pi,strS(pi),endS(pi),j', &
!                             pi,strS(pi),endS(pi), j, j - strS(pi)+ 1
                  if (j >= strS(pi) .and. j <= endS(pi)) then
                    jl = j - strS(pi)+ 1
                    pi_j = pi-1
                    exit
                  end if
!                  write(*,*) 'pi,strS(pi),endS(pi),j', &
!                             pi,strS(pi),endS(pi), j, j - strS(pi)+ 1
                end do
                dissM_j = dissMetric(jl)
!                write(*,*) '3j is not in the current PI, ', &
!                           'myid, j,jl, dissM_j',&
!                           myid,j,jl, dissM_j 
              end if  
              call MPI_Gather(dissM_j, 1, MPI_DOUBLE_PRECISION, &
                              rbuf_j, 1, MPI_DOUBLE_PRECISION,0, &
                              MPI_COMM_WORLD, ierr)
!
! the clustering load is going to be on the root
              if (myid == 0) then
!                write(fu_in2,*) 'rbuf_i',rbuf_i,'rbuf_j',rbuf_j
!                write(fu_in2,*) 'i,il,j,jl,pi_i,pi_j',i,il,j,jl,pi_i,pi_j
!
!     Get pairs with st2Id, that haven't been clustered before 
                if (i == il) then
                  dissM_i = dissMetric(i)
                else
                  dissM_i = rbuf_i(pi_i+1)
                end if
!
!   Get pairs with st1Id, that haven't been clustered before
                if (j == jl) then
                  dissM_j = dissMetric(j)
                  dissMetric_tmp = 0.5*(dissM_i+dissM_j)
!                  write(fu_in2,*) 'j,dissM_i,dissM_j,dissMetric', &
!                              j,dissM_i,dissM_j,dissMetric_tmp
                else
                  dissM_j = rbuf_j(pi_j+1)
                  dissMetric_tmp =  0.5*(dissM_i+dissM_j)
!                  write(fu_in2,*) 'j,dissM_i,dissM_j,dissMetric_tmp', &
!                              j,dissM_i,dissM_j,dissMetric_tmp
                end if
!
              end if ! if (myid=0)
              call MPI_Barrier(MPI_COMM_WORLD, ierr)
              call MPI_Bcast(dissMetric_tmp, 1, MPI_DOUBLE_PRECISION, root, &
                             MPI_COMM_WORLD, ierr)
!
!   ndex j is the location where the value in the dissMetric will change
!              write(fu_in2,*) 'myid,pi_j,jl,dissMetric_tmp', &
!                         myid,pi_j,jl,dissMetric_tmp

              if (myid == pi_j) dissMetric(jl) = dissMetric_tmp
               
            end do
!            write(fu_in2,*) 'A: dissMetric',dissMetric
!
            call MPI_Barrier(MPI_COMM_WORLD, ierr)
!
          end do ! end clustering cycle
!
!          if (writeOutC) then
            do i = 1,inc
              write(fu_in,*) xx(i),lki(i),lkj(i)
            end do
            close(fu_in)
            close(fu_in1)
            close(fu_in2)
!          end if
!
          deallocate(minDloc, lloc)
          deallocate(dissMetric)
          deallocate(w, wj, lki, lkj, xx)
          deallocate(station1,station2)
          deallocate(rbuf_i,rbuf_j,rbuf_1,rbuf_2)
          deallocate(id4cl)
!                    
        end if IF_CLUSTERS_NAVAIL
        end_taskS = MPI_Wtime()
        write(*,*) 'PI',myid,'compute clustering running time (s)' , &
                    end_taskS-start_taskS
!       
        deallocate(strS, endS)

      end do !cycling over number of species    
!
      write(6,*)'pab finished pid# ',myid
        call MPI_FINALIZE(ierr)

              
      contains

!
!##########################################################################################
!
      subroutine d2j(dat,julian,ierr)
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! * Author:    John S. Urban
      ! * Version:   1.0 2015-12-21
      ! * Reference: From Wikipedia, the free encyclopedia 2015-12-19
      ! * There is no year zero
      ! * Julian Day must be non-negative
      ! * Julian Day starts at noon; while Civil Calendar date starts at midnight
      !-----------------------------------------------------------------------------------------------------------------------------------
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
      !-----------------------------------------------------------------------------------------------------------------------------------
      ! * Author:    John S. Urban
      ! * Version:   1.0 2015-12-21
      ! * Reference: From Wikipedia, the free encyclopedia 2015-12-19
      !-----------------------------------------------------------------------------------------------------------------------------------

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

      !call date_and_time(values=timezone)
      !tz = timezone(4)     

      ijul=int(julian)                           ! Integral Julian Day
      !second=real(julian-dble(ijul)*secday)      ! Seconds from beginning of Jul. Day
      !second=second+(tz*60)

      !if(second.ge.(secday/2.0d0)) then            ! In next calendar day
      !  ijul=ijul+1
      !  second=second-(secday/2.0d0)              ! Adjust from noon to midnight
      !else                                         ! In same calendar day
      !  second=second+(secday/2.0d0)              ! Adjust from noon to midnight
      !endif

      !if(second.ge.secday) then                    ! Final check to prevent time 24:00:00
      !  ijul=ijul+1
      !  second=second-secday
      !endif

      !minute=int(second/60.0)                      ! Integral minutes from beginning of day
      !second=second-float(minute*60)               ! Seconds from beginning of minute
      !hour=minute/60                               ! Integral hours from beginning of day
      !minute=minute-hour*60                        ! Integral minutes from beginning of hour

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
      !datAr(4)=tz
      !datAr(5)=hour
      !datAr(6)=minute
      !datAr(7)=int(second)
      !datAr(8)=int((second-int(second))*1000.0)
      dat = datAr(1)*10000+datAr(2)*100+datAr(3)
      ierr=0

      return
      end subroutine j2d

!
!************************************************************************************
!
      function map(n,i,j)
       !map row i and column j of upper half diagonal symmetric matrix onto vector
        integer(i8),intent(in)  :: n,i,j
        integer,parameter   :: i8 = selected_int_kind(15)
        integer(i8)         :: map

        map = (j + (i-1_i8)*n - (i*(i+1_i8))/2_i8)

      end function map
!
!************************************************************************************
!
      end program kzfilter_hc
