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
! mtype = 1
!Select standardization by range:
! 0 = no
! 1 = yes
! istd = 0
!select proximity according to dis or dimilarity
! 1 => dissimilarity
! -1 ==> similarity
! is = -1, if mtype=2,4 the is=1
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


      implicit none
       
       integer                                 :: mtype,itmp,is,istd, &
                                                  betalw,beta,itest, &
                                                  iclus1,iclus2,ncol, &
                                                  ir,inf,jr,i,j,l,k, &
                                                  nts,nst,nst2,w0, &
                                                  i1,nCells,kz,mkz,pkz, & 
                                                  nsum,fu_bg,key,ier, &
                                                  igc,ict,iIni,iEnd, &
                                                  jIni,jEnd,fcst_hr, &
                                                  iday,ithr
       integer                                 :: start_date, end_date

       real(kind=4)                            :: freq,period,avg, &
                                                  datmiss
       character(len=140)                      :: header
       character(len=40)                       :: headername
       character(len=14)                       :: cdum
       character(len=35)                       :: filein,obsheader
       integer,dimension(:),allocatable        :: id !,lk,lki,lkj 
       !real(kind=4),dimension(:),allocatable   :: lat,lon,xx,dissMatrix
       !integer(kind=4),dimension(:),allocatable:: stId1,stId2
       !integer(kind=4),dimension(:,:),allocatable :: cluster_table
       real(kind=4),allocatable,dimension(:)   :: i_array,j_array
       real(kind=4),allocatable,dimension(:,:) :: data_st,data_stT,varbg
       !areal(kind=4),dimension(:,:),allocatable :: kzout,data_in
       character(len=70) :: obs_file, outfile
!
       real(kind=4)                            :: i_stack,j_stack, &
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
! 4 FST FILES
       integer                     :: bggrid_npts
       integer                     :: dateo_bidon, deet_bidon, npas_bidon
       integer                     :: bggrid_ni, bggrid_nj, &
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
      real               :: start,finish,end_task(2), &
                            start_task(2)
      character(len=150) :: timeFNm

      logical, parameter :: ascii_out = .false.
      logical, parameter :: centroids = .false.
      character(len=1)   :: dataAvail,distMavail,clustersAvai,cntroidsAvail
      logical            :: data_not_available,distMatrix_not_available, &
                            clusters_not_available,centroids_not_available
       
      integer :: narg,stopt
      character(len=10) :: nstS,stOptS,lonS,latS
      character(len=1) :: metricS,istdS,itmpS
     
      integer(kind=4),dimension(8) :: date2julianS,date2julianE
      real                         :: jDayS,jDayE
      integer                      :: err
      character(len=8)             :: start_dateS,start_dateE,metric

!  Model surface ip1 values.
!
      !integer, parameter :: ip1g_surface = 93423264
      !integer, parameter :: ip1t_surface = 76696048
      !integer, parameter :: ip1m_surface = 75597472
      integer :: ip1g_surface,ip1t_surface,ip1m_surface
!
! parallel 
      integer, external :: omp_get_max_threads,omp_get_thread_num, &
                           omp_get_num_threads
!
! FST functions
      integer,external ::  fstinf, fstprm, fstlir, ezqkdef, gdrls, newdate, &
                           ezsetopt, fnom, fstouv, fstfrm, fclos, gdxyfll, &
                           gdllsval, gdllfxy, fstopc
!    
!  general settings
! 
      ip1g_surface = 93423264
      ip1t_surface = 76696048
      ip1m_surface = 75597472
!
      narg=iargc()
      call getarg(1,nstS)
      read(nstS,'(i)') nst
      if (mod((sqrt(nst*1.0)-1),2.0) /= 0) then
        write(6,*) 'The remainder of sqrt(nst)-1)/2 must be 0,',&
                   'check nst',nst
        write(6,*) 'exiting program'
        stop
      end if
      call getarg(2,metricS)
      read(metricS,'(i)') mtype
      !call getarg(3,istdS)
      !read(istdS,'(i)') istd
      is = 1                  ! similarity vs dissimilarity
      if (mtype == 3) is = -1
      call getarg(3,itmpS)
      read(itmpS,'(i)') itmp  ! clustering method
      if (itmp < 8) then
        betalw = 1
        beta = 1
      else 
        betalw = 0
        beta = 0 ! values needs to be between -1. and 1.
      end if
      !lower and upper bound of the clusering
      iclus1 = 2
      iclus2 = nst-1
      ! how many clusters/stations for dissimilarity output
      call getarg(4,stOptS)
      read(stOptS,'(i)') stOpt
     ! central-point of the domain
      call getarg(5,lonS)
      read(lonS,'(f10.5)') lonPS
      call getarg(6,latS)
      read(latS,'(f10.5)') latPS
     ! data existing or not
      call getarg(7,dataAvail)
      read(dataAvail,'(L1)') data_not_available
      call getarg(8,distMavail)
      read(distMavail,'(L1)') distMatrix_not_available
      call getarg(9,clustersAvai)
      read(clustersAvai,'(L1)') clusters_not_available
      call getarg(10,cntroidsAvail)
      write(6,*) 'cntroidsAvail',cntroidsAvail
      read(cntroidsAvail,'(L1)') centroids_not_available
      write(6,*) 'checks',data_not_available,distMatrix_not_available,&
                  clusters_not_available,centroids_not_available
!
!  GEM-MACH output
!
!      inDir = '/fs/site2/dev/eccc/aq/r1/jso001/'//&
!              'project/network_analysis/data/annual_run_extracts/'
      inDir = '/space/hall4/sitestore/eccc/aq/r1/jso001/'//&
           'project/network_analysis/data/annual_run_extracts/'
!'/users/tor/arqp/dej/newcetus3/annual_run_extracts/'
      fcst_hr = 18
      header = 'GEM-MACH Oil sans run between 09.2013-08.2014'
      date2julianS = (/2013,08,1,fcst_hr,00,00,00,00/)
      date2julianE = (/2013,08,29,fcst_hr,00,00,00,00/)
      call d2j(date2julianS,jDayS,err)
      call d2j(date2julianE,jDayE,err)
      nts = int((jDayE-jDayS+1)*24)
      start_date = date2julianS(1)*10000+date2julianS(2)*100+date2julianS(3)
      end_date = date2julianE(1)*10000+date2julianE(2)*100+date2julianE(3)
      write (dateIniStr,'(i8)') start_date ! converting integer to string using a 'internal file'
      write (dateEndStr,'(i8)') end_date ! converting integer to string using a 'internal file'
      write(fhr,'(i2.2)') fcst_hr
! output directory and parameters to compute HC   
!      outDir ='/home/ords/aq/jso001/project/network_analysis/'//& 
!              'phase4/'
      outDir = 'output/'
!
! ################  starting ############################ 
!
      write(6,*) '!-------------------------------------!'
      write(6,*) 'Hi world! Network analysis in progress.'
      write(6,*) '!-------------------------------------!'

      do l = 1,nbg
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
        write(6,*) 'start_date: ',dateIniStr,', end_date: ',dateEndStr
        write(6,'(a25,1x,a5,1x,a6)') ' central point (lat,lon):',latS,lonS
        write(6,*) 'species: ',trim(varbg_name(l))
        write(6,*) 'time series length:',nts, 'no. stations:', nst
        write(6,*) 'metric(1-5),standardization(0-no,1-yes):',mtype,istd
        write(6,*) 'clustering method (1-9):,',itmp
        write(6,*) 'no. clusters for optimization:',stOpt
        headername = '_modelPhase4_v23_'
        outFNmData = trim(outDir)//trim(varbg_name2(l))//&
                 trim(headername)
        outFNm = trim(outDir)//trim(varbg_name2(l))//&
                 trim(headername)//trim(metric)//'_'

        write(6,*) 'filename:',trim(headername)
        write(6,*) 'output filename: ',trim(outFNm)
!  Allocate global resources
        nCells = map(nst,nst-1,nst) 
! Get information about the model data
!       ### turn off FST information messages
        ier = fstopc('MSGLVL', 'SYSTEM', 0)
        fu_bg = 15 
        write(model_date, '(i8.8)') start_date
        write(tstep, '(i6.6)') 60
        bg_filename = model_date//trim(fhr)//'_'//trim(tstep)//'p'
        bg_filename = trim(inDir)//model_date//trim(fhr)//&
                      '/'//trim(bg_filename)
!    
!  Reading data if the GEM-MACH data is not available in a binary format 

       if (data_not_available) then
!   Open input files for extracting information not dependent of time step,
         ier = fnom(fu_bg, trim(bg_filename), 'RND+OLD+R/O', 0)
         if (ier < 0) then
           print *,'Ending,error while opening ',trim(bg_filename)
           stop
         end if
         ier = fstouv(fu_bg, 'RND')
!   Reading grid parameters
         bggrid_def_field = trim(varbg_name(l))
         key = fstinf(fu_bg,bggrid_ni,bggrid_nj, nk_bidon, -1, ' ',&
                     -1, -1, -1, ' ', bggrid_def_field)
         if (key < 0) then
           write(6,*) 'Key value negative after opening file',&
                      ' with grid definitions'
           stop
         end if
         !nbggrid_npts = bggrid_ni * bggrid_nj
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
        WRITE(*,*) 'Initializing grid at ', bggrid_ni, 'x', bggrid_nj
!    Initialize the bg grid for ezscint functions
         ier = ezsetopt('VERBOSE', 'NO')
         bggrid_id = ezqkdef(bggrid_ni,bggrid_nj,bggrid_type,bggrid_ig1,&
                             bggrid_ig2, bggrid_ig3, bggrid_ig4, fu_bg)
         ier = gdxyfll(bggrid_id,i_stack,j_stack,latPS,lonPS + 360.0,1) !!! REMOVE
         !
         write(*,*) 'gdxyfll returned ', i_stack, ', ', j_stack
         allocate(varbg(bggrid_ni,bggrid_nj))          
         ier = fstlir(varbg(:,:),fu_bg,ni_bidon,nj_bidon,nk_bidon , &
                      -1, '',ip1t_surface, -1, -1, '', bggrid_def_field)
         deallocate(varbg)
         igc = 0
         nst2 = floor(sqrt(real(nst))*0.5)
         iIni = i_stack-nst2
         iEnd = i_stack+nst2
         jIni = j_stack-nst2
         jEnd = j_stack+nst2
         allocate(id(nst)) !i_array(nst),j_array(nst))
         ! write(6,*) 'id(igc),i_array(igc),j_array(igc)'
!         do j = jIni,jEnd
!          do i = iIni,iEnd
!            igc=igc+1
!            id(igc) = igc
!            i_array(igc) = real(i)
!            j_array(igc) = real(j)
!            write(6,*) id(igc),i_array(igc),j_array(igc)
!          end do
!        end do          
         do j = -nst2,nst2
           do i = -nst2,nst2
             igc = igc+1
               id(igc) = igc
               !i_array(igc) = (i_stack+i) * 1.0
               !j_array(igc) = (j_stack+j) * 1.0
           end do
         end do
        ! ier = gdllfxy(bggrid_id,lat,lon,i_array,j_array,nst)
         !lon = lon-360.0
         ier = fstfrm(fu_bg)
         ier = fclos(fu_bg)
!      ### save data
         outFNmDataIn = trim(outFNmData)//'id_lat_lon_'// &
                        trim(nstS)//'gc.bin'
         open(unit=31,file=trim(outFNmDataIn),form='unformatted', &
                                              action='write')
         do igc = 1,nst
           write(31) id(igc) !,lat(igc),lon(igc)
         end do
         close(31)
         deallocate(id) !,lon,lat)
         !deallocate(i_array,j_array)
         call dtime(start_task)
         allocate(data_st(nst,nts))
         write(6,*) 'reading GEM-MACH data'
!   model data for each time step
         iday = 0
         cur_date = start_date
         do while (cur_date <= end_date)
           write(model_date, '(i8.8)') cur_date
           iday = iday+1 !#counter for the day
           do process_step = 1, 24
             model_tstep = process_step
             ict = 24*(iday-1)+process_step
             write(tstep, '(i6.6)') model_tstep * 60
             bg_filename = model_date//trim(fhr)//'_'//trim(tstep)//'p'
             bg_filename = trim(inDir)//model_date//trim(fhr)//&
                                '/'//trim(bg_filename)
             fu_bg = 10+process_step
!           Open input files
             ier = fnom(fu_bg, trim(bg_filename), 'RND+OLD+R/O', 0)
             if (ier < 0) then
               print *,'Ending,error while opening ',trim(bg_filename)
               stop
             end if
             ier = fstouv(fu_bg, 'RND')
             allocate(varbg(bggrid_ni,bggrid_nj))          
!   Read in the requested field for this model timestep
             ier = fstlir(varbg(:,:), fu_bg, ni_bidon, nj_bidon, nk_bidon, &
                          -1, '', ip1t_surface, -1, -1, '', varbg_name)
             igc = 0
             do j = -nst2,nst2 !jIni,jEnd
               do i = -nst2,nst2 !iIni,iEnd
                 igc=igc+1
                 data_st(igc,ict) = varbg(i_stack+i,j_stack+j) !varbg(i,j)
!                 write(6,*) 'igc,ict,varbg(i,j)',igc,ict,varbg(i,j)
               end do
             end do
             deallocate(varbg)
!
             ier = fstfrm(fu_bg)
             ier = fclos(fu_bg)
           end do !# process_step, read hourly files
!
           call incdat(idate_next, dateo_bidon, 24*iday)
           ier = newdate(idate_next, next_date, dummy_hour, -3)
           cur_date = next_date
         end do !# do while, period requested
! Transpose the matrix for metric calculation (code optimization) 
         allocate(data_stT(nts,nst))
         WRITE(*,*) ' transposing data_st'
        !$OMP WORKSHARE
         data_stT = transpose(data_st) 
        !$OMP END WORKSHARE
         WRITE(*,*) ' done transposing data_st'
         deallocate(data_st)
         call dtime(end_task)
         write(6,*) 'time to read data:',end_task(1)-start_task(1),end_task(2)-start_task(2)
!        ### save data
         outFNmDataIn = trim(outFNmData)//trim(nstS)//'gc.bin'
         open(unit=31,file=trim(outFNmDataIn),form='unformatted', &
                                              action='write')
         do igc = 1,nst
           do ict = 1,nts
             write(31) data_stT(ict,igc)         
           end do
         end do 
         close(31)
         deallocate(data_stT)
       end if ! if data_stT is not available
!!  Done reading hourly data. Start clustering!
!
       do kz=1,nkz
!  Reading time series as is
         if (kzname(kz) == 'kz_001_001') then
           !write(6,*) 'centroids_not_available',centroids_not_available
           if (centroids_not_available) then  !! needs input from sub clustering
!              write(6,*) 'clusters_not_available',clusters_not_available
             if (clusters_not_available) then !! needs input from sub compute_metric
               if (distMatrix_not_available) then 
!!  Compute dissimilarity matrix and save the matrix and stId1,stId2
                 call dtime(start_task)
                 if (mtype == 4) then
                   write(6,*) 'computing distance matrix'
                   call compute_metric_1_R(nst,nts,nCells,is, &
                                           outFNmData,outFNm,nstS)
                 end if
                 call dtime(end_task)
                 write(6,*) 'compute distance matrix running time', &
                            end_task(1)-start_task(1),end_task(2)-start_task(2)
               end if ! distMatrix_not_available
!!  Cluster and recompute the metric; save the parameters lkj,lki,lk,xx          
               call dtime(start_task)
               write(6,*) 'clustering'
               ! if (itmp == 4) then
                 call clustering_ave(nst,nCells,outFNm,nstS)
               ! end if
               call dtime(end_task)
               write(6,*) 'compute clustering running time', &
                          end_task(1)-start_task(1),end_task(2)-start_task(2)
             end if ! clusters_not_available
!!  Compute centroids if not available, save the parameters cluster_table      
             call dtime(start_task)
             write(6,*) 'computing centroids'
             call compute_centroids(nst,nCells,&
                                    iclus1,iclus2,outFNm,nstS)
             call dtime(end_task)
             write(6,*) 'compute centroids running time', &
                         end_task(1)-start_task(1),end_task(2)-start_task(2)
           end if !# if centroid_available
!! Write the output          
           call dtime(start_task)
           write(6,*) 'writing output'
           outFNmDataIn = trim(outFNm)//'id_lat_lon_'// &
                          trim(nstS)//'gc.bin'
           call write_output(nst,outFNmDataIn,outFNm, &
                             varbg_name(l), &
                             mtype,stOpt,inDir,start_date, &
                             latPS,lonPS,ip1t_surface,fhr)
           call dtime(end_task)
           write(6,*) 'write output running time', &
                        end_task(1)-start_task(1),end_task(2)-start_task(2)
         else ! kz > 001
!! time series will be fitered and then clustered
! Filter out (100% removal) all time signals with < 1 day duration, 50% pass on 2.75 day duration
! 100% pass on 42 day duration (the frequency returned is for the 50% level; multiply by 2.82 to get the 0%
! high frequency cut-off, and by 0.0657 to get the low frequency 100% pass levels):
           !allocate(kzout(nst,nts))
           !! Filter out all periods less than 1 day:
           if (kzname(kz) == 'kz_017_003') then 
             mkz = 17
             pkz = 3
           !! Filter out all periods less than 1 week:
           elseif (kzname(kz) == 'kz_095_005') then
             mkz = 95
             pkz = 5
           !! Filter out all periods less than 30 days :
           elseif (kzname(kz) == 'kz_523_003') then
             mkz = 523
             pkz = 3
           end if
!            w0=795
           i1= 8000 ! remove this eventually
           do k = 1,nst
             do j = 1,w0
              ! kzout(j,k) = 0.0
             end do
             !call kzf(data_in(1,k),nts,mkz,pkz,w0,i1,kzout(1,k),freq,period)
           end do
!!!!!!!!!!!!!!! STUFF MISSING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !deallocate(kzout)
         end if ! if kz > 1
!
        end do ! ### ikz 
!
      end do !### l, for each species 
       
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
      integer,intent(out)        :: ierr     ! Error return, 0 for successful execution,-1=invalid year,-2=invalid month,-3=invalid day,
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
       subroutine compute_metric_1_R(n,nt,nv,is,inFNm,outFNm,nstS)

       integer,intent(in)                         :: nt,is       ! no. records, (di)similarity
       integer(kind=4),intent(in)                 :: n,nv        ! no points,hours
       character(len=256),intent(in)              :: inFNm,outFNm
       character(len=10),intent(in)               :: nstS
       !local
       integer                          :: ihr,idThr,numThr,nThr, &    ! parallel
                                           omp_get_num_threads, &  
                                           omp_get_thread_num
       integer(kind=4), &
             allocatable,dimension(:)   :: station1,station2 ! station pairs ids
       integer(kind=4)                  :: gc,i,j,ipc 
       real(kind=4)                     :: ntr,isr,sxy,sxx,syy, &      ! temp vars
                                           data_st1,data_st2, &
                                           data_stTmp,sum_st1st2Tmp, &
                                           sum_stTmp,sum_stSQRTmp, &
                                           varTmp1,varTmp2,sumTmp, &
                                           start_task,end_task
       real(kind=4), &
             allocatable,dimension(:)   :: sum_st,sum_stSQR,metric ! precomputions 4 metric
       real(kind=4), &
             allocatable,dimension(:,:) :: data_stT                ! time series
       character(len=256)               ::fnm
!
! Get time series for each grid cell
       fnm = trim(inFNm)//trim(nstS)//'gc.bin'
       open(unit=31,file=trim(fnm),form='unformatted',action='read')
       allocate(data_stT(nt,n))
      !$OMP PARALLEL DO PRIVATE(ict,igc)
       do ict = 1,nt
         do igc = 1,n
           read(31) data_stT(ict,igc)
         end do
       end do 
      !$OMP END PARALLEL DO 
       close(31)   
!  compute the summation of a single timeseries for facilitating the metric computation later on
       allocate(sum_st(n),sum_stSQR(n))
       sum_st =  0.
       sum_stSQR = 0.
      !$OMP PARALLEL DO PRIVATE(igc,sum_stTmp,sum_stSQRTmp,ihr,data_stTmp)
       do igc = 1,n
         sum_stTmp = 0.
         sum_stSQRTmp = 0.
         do ihr = 1,nt,12
           data_stTmp = data_stT(ihr,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+1,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+2,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+3,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+4,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+5,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+6,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+7,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+8,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+9,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+10,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
!
           data_stTmp = data_stT(ihr+11,igc)
           sum_stTmp = sum_stTmp+data_stTmp
           sum_stSQRTmp = sum_stSQRTmp+data_stTmp*data_stTmp
         end do
         sum_st(igc) = sum_stTmp
         sum_stSQR(igc) = sum_stSQRTmp
       end do
      !$OMP END PARALLEL DO
!  
!   Trying to map the squared matrix: get all station pairs available
       allocate(station1(nv),station2(nv))
      !$OMP PARALLEL DO PRIVATE(i,j,ipc)
       do i = 1,n-1
         do j = i+1,n
           ipc = j + (i-1)*n - (i*(i+1))/2
           station1(ipc) = i
           station2(ipc) = j
         end do
       end do
      !$OMP END PARALLEL DO 

! Compute intermediate sums for computing metric and the metric
       ntr = 1/real(nt)
       isr = real(is)
       allocate(metric(nv))
       metric = 0.
      !$OMP PARALLEL DO PRIVATE(ipc,sum_st1st2Tmp,varTmp1,varTmp2,ihr) 
       do ipc = 1,nv
! Compute intermediate sums for computing metric below
         sum_st1st2Tmp = 0.
         varTmp1 = station1(ipc)
         varTmp2 = station2(ipc)
         do ihr = 1,nt,12
           sum_st1st2Tmp = sum_st1st2Tmp + &
                           data_stT(ihr,varTmp1) * &
                           data_stT(ihr,varTmp2) + &
                           data_stT(ihr+1,varTmp1) * &
                           data_stT(ihr+1,varTmp2) + &
                           data_stT(ihr+2,varTmp1) * &
                           data_stT(ihr+2,varTmp2) + &
                           data_stT(ihr+3,varTmp1) * &
                           data_stT(ihr+3,varTmp2) + &
                           data_stT(ihr+4,varTmp1) * &
                           data_stT(ihr+4,varTmp2) + &
                           data_stT(ihr+5,varTmp1) * &
                           data_stT(ihr+5,varTmp2) + &
                           data_stT(ihr+6,varTmp1) * &
                           data_stT(ihr+6,varTmp2) + &
                           data_stT(ihr+7,varTmp1) * &
                           data_stT(ihr+7,varTmp2) + &
                           data_stT(ihr+8,varTmp1) * &
                           data_stT(ihr+8,varTmp2) + &
                           data_stT(ihr+9,varTmp1) * &
                           data_stT(ihr+9,varTmp2) + &
                           data_stT(ihr+10,varTmp1) * &
                           data_stT(ihr+10,varTmp2) + &
                           data_stT(ihr+11,varTmp1) * &
                           data_stT(ihr+11,varTmp2)
         end do
         data_st1 = sum_st(varTmp1)
         data_st2 = sum_st(varTmp2)
         SXY = sum_st1st2Tmp - (data_st1*data_st2*ntr)
         SXX = sum_stSQR(varTmp1)-(data_st1*data_st1*ntr)
         SYY = sum_stSQR(varTmp2)-(data_st2*data_st2*ntr)
         metric(ipc) = (1.0-SXY/SQRT(SXX*SYY))*isr         
       end do
      !$OMP END PARALLEL DO
       deallocate(data_stT,sum_st,sum_stSQR) 
       write(6,*) 'metric',metric
!
!   ### save data 
       fnm = trim(outFNm)//trim(nstS)//'gc_distMatrix.bin'
       open(unit=31,file=trim(fnm),form='unformatted',action='write')
       do igc = 1,nv
         write(31) metric(igc),station1(igc),station2(igc)
       end do
       close(31)
       deallocate(metric,station1,station2)
!       
       return 

      end subroutine compute_metric_1_R

!
!#######################################################################
!
      subroutine clustering_ave(n,nv,outDir,nstS)

      integer(kind=4),intent(in)         :: n,nv      ! parameters
      character(len=256),intent(in)      :: outDir
      character(len=10),intent(in)       :: nstS

      ! local
      integer                            :: numThr,nThr,idThr         ! parallel cycles
      integer(kind=4)                    :: i,l,ltmp,icl,icl2,ipc,  & ! counters, temp vars 
                                            ki,ii,kk,ik,ji,j,jj,igc, &
                                            st1Tmp,st2Tmp,jpl1,jpl2,&
                                            inc   
      integer(kind=4),allocatable, &
                   dimension(:)          :: strT,endT,lloc            !interval, local l, parallel
      real(kind=4),allocatable, &
              dimension(:)               :: minDloc                   ! local min, parallel
      integer(kind=4),allocatable, &
              dimension(:)               :: st1,st2                   ! paired stations ids
      real(kind=4),allocatable, &
              dimension(:)               :: D                         ! paired stations metric
      integer(kind=4),allocatable, &
              dimension(:)               :: lki,lkj,lk                ! tracking clusters
      real(kind=4),allocatable, &
              dimension(:)               :: xx                        ! nodes
      real(kind=4)                       :: x,minD,perR               ! temp var, percentile
      real(kind=4)                       :: start_task,end_task, &    ! timming
                                            start_taskG,end_taskG
      logical,allocatable, &
              dimension(:)               :: w,wj  ! tracking clustered pairs, station index
      character(len=256)                 :: fnm

!
!
! Reading station pairs and its metric
      fnm = trim(outDir)//trim(nstS)//'gc_distMatrix.bin'
      allocate(D(nv),st1(nv),st2(nv))
      open(unit=31,file=trim(fnm),form='unformatted',action='read')
      do igc = 1,nv
        read(31) D(igc),st1(igc),st2(igc)
      end do
      close(31)
      call percentiler(D,nv,0.9,perR)
!
! Define interval for open MP computations 
     !$OMP PARALLEL SHARED(numThr)
       numThr = omp_get_num_threads()
     !$OMP END PARALLEL
      inc = nv/numThr
      allocate(minDloc(numThr),lloc(numThr))
      allocate(strT(numThr),endT(numThr))
      strT = 1
      endT = nv
      do idThr = 1,numThr-1
        endT(idThr) = strT(idThr) + inc - 1.0
        strT(idThr+1) = strT(idThr) + inc
      end do
!
      allocate(xx(nv),w(nv))
      allocate(wj(n),lki(n),lkj(n))
      XX = 0.
      w = .false.
      wj = .false.
      lki = 0
      lkj = 0
!   Clustering: 1cl=1 (avoid if statements) 
!               icl =2,ntmp(to where D>=P), icl2 = ntmp,n-1 
      icl = 1
!   Get the minimum value for the metric, keeping track if the pair
!   has been clustered or not via  w(l)
      minDloc = 999999999.
      lloc = -999999999
     !$OMP PARALLEL PRIVATE(idThr,l)
      idThr = omp_get_thread_num()+1
      do l = strT(idThr),endT(idThr)
        if (w(l)) cycle                        ! doesn't read clustered pairs
        if (D(l) >= perR) cycle                ! takes only values lower than perR
        if (wj(st1(l)) .or. wj(st2(l))) cycle  ! does't read if satation was already clustered
        if (D(l) > minDloc(idThr)) cycle
        minDloc(idThr) = D(l)
        lloc(idThr) = l
      end do
     !$OMP END PARALLEL
      minD = minval(minDloc)
      do idThr = 1,numThr
        if(minDloc(idThr) == minD) ltmp = lloc(idThr)
      end do
      XX(icl) = minD
      st1Tmp = st1(lTmp) ! j
      st2Tmp = st2(lTmp) ! i
      lki(icl) = st2Tmp
      lkj(icl) = st1Tmp
      wj(st1Tmp) = .true.
      w(lTmp) = .true.
!   Taking only the pairs of stations that cluster with all the other stations available
     !$OMP PARALLEL DO PRIVATE(igc,jpl1,jpl2,i,j)
      do igc = 1,n
        if (wj(igc)) cycle
        if (igc == st2Tmp) cycle
        jpl1 = min(st1Tmp,igc)
        jpl2 = max(st1Tmp,igc)
        !i = map(n,jpl1,jpl2) !map(n,i,j)=j + (i-1)*n - (i*(i+1))/2 , i<j
        i = jpl2 + (jpl1-1)*n - (jpl1*(jpl1+1))/2
        jpl1 = min(st2Tmp,igc)
        jpl2 = max(st2Tmp,igc)
        j = jpl2 + (jpl1-1)*n - (jpl1*(jpl1+1))/2
        ! D(j) = 0.5*(D(i)+D(j))
        D(j) = min(D(i),D(j))
      end do
     !$OMP END PARALLEL DO 
!
      do icl = 2,n-1
        minDloc = 999999999.
        lloc = -999999999
       !$OMP PARALLEL PRIVATE(idThr,l)
        idThr = omp_get_thread_num()+1
        do l = strT(idThr),endT(idThr)
          if (w(l)) cycle                     
          if (D(l) >= perR) cycle
          if (wj(st1(l)) .or. wj(st2(l))) cycle
          if (D(l) > minDloc(idThr)) cycle
          minDloc(idThr) = D(l)
          lloc(idThr) = l
        end do
       !$OMP END PARALLEL
        minD = minval(minDloc)
        do idThr = 1,numThr
          if(minDloc(idThr) == minD) ltmp = lloc(idThr)
        end do
        if (minD == 999999999.) exit
        XX(icl) = minD
        st1Tmp = st1(lTmp) ! j
        st2Tmp = st2(lTmp) ! i
        lki(icl) = st2Tmp
        lkj(icl) = st1Tmp
        wj(st1Tmp) = .true. 
        w(lTmp) = .true.    
       !$OMP PARALLEL DO PRIVATE(igc,jpl1,jpl2,i,j)
        do igc = 1,n
          if (wj(igc)) cycle
          if (igc == st2Tmp) cycle
          jpl1 = min(st1Tmp,igc)
          jpl2 = max(st1Tmp,igc)
          !i = map(n,jpl1,jpl2) !map(n,i,j)=j + (i-1)*n - (i*(i+1))/2 , i<j
          i = jpl2 + (jpl1-1)*n - (jpl1*(jpl1+1))/2
          jpl1 = min(st2Tmp,igc)
          jpl2 = max(st2Tmp,igc)
          j = jpl2 + (jpl1-1)*n - (jpl1*(jpl1+1))/2
          ! D(j) = 0.5*(D(i)+D(j))
          D(j) = min(D(i),D(j))
        end do
       !$OMP END PARALLEL DO
!
!      ### save data
       if (mod(icl,50000) == 0) then
         write(6,*) 'Storing data due to time constraints, icl=',icl
         fnm = trim(outFNm)//'lkj_lki_lk_xx_wj_tmp.bin'
         open(unit=31,file=trim(fnm),form='unformatted', &
                                            action='write')
         do igc = 1,n
           write(31) lkj(igc),lki(igc),lk(igc),xx(igc),wj(igc)
         end do
         close(31)
         fnm = trim(outFNm)//'w_D_tmp.bin'
         open(unit=31,file=trim(fnm),form='unformatted', &
                                     action='write')
         do igc = 1,nv
           write(31) w(igc),D(igc)
         end do
         close(31)
       end if  
      end do
!
      do icl2 = icl,n-1
        minDloc = 999999999.
        lloc = -999999999
       !$OMP PARALLEL PRIVATE(idThr,l)
        idThr = omp_get_thread_num()+1
!       !$OMP DO 
        do l = strT(idThr),endT(idThr)
          if (w(l)) cycle                        !ensures that it doesn't read clustered pairs
          if (wj(st1(l)) .or. wj(st2(l))) cycle
          if (D(l) > minDloc(idThr)) cycle
          minDloc(idThr) = D(l)
          lloc(idThr) = l
        end do
!       !$OMP END DO
        !write(6,*) 'icl,idThr,minDloc',icl,idThr,minDloc(idThr)
       !$OMP END PARALLEL
        minD = minval(minDloc)
        do idThr = 1,numThr
          if(minDloc(idThr) == minD) ltmp = lloc(idThr)
        end do
        XX(icl2) = minD
        st1Tmp = st1(lTmp) ! j
        st2Tmp = st2(lTmp) ! i
        lki(icl2) = st2Tmp
        lkj(icl2) = st1Tmp
        wj(st1Tmp) = .true.
        w(lTmp) = .true.
       !$OMP PARALLEL DO PRIVATE(igc,jpl1,jpl2,i,j)
        do igc = 1,n
          if (wj(igc)) cycle
          if (igc == st2Tmp) cycle
          jpl1 = min(st1Tmp,igc)
          jpl2 = max(st1Tmp,igc)
          i = jpl2 + (jpl1-1)*n - (jpl1*(jpl1+1))/2
          jpl1 = min(st2Tmp,igc)
          jpl2 = max(st2Tmp,igc)
          j = jpl2 + (jpl1-1)*n - (jpl1*(jpl1+1))/2
          ! D(j) = 0.5*(D(i)+D(j))
          D(j) = min(D(i), D(j))
        end do
       !$OMP END PARALLEL DO 
!      ### save data
       if (mod(icl,50000) == 0) then
         write(6,*) 'Storing data due to time constraints, icl=',icl
         fnm = trim(outFNm)//'lkj_lki_lk_xx_wj_tmp.bin'
         open(unit=31,file=trim(fnm),form='unformatted', &
                                            action='write')
         do igc = 1,n
           write(31) lkj(igc),lki(igc),lk(igc),xx(igc),wj(igc)
         end do
         close(31)
         outFNmData = trim(outFNm)//'w_D_tmp.bin'
         open(unit=31,file=trim(outFNm),form='unformatted', &
                                            action='write')
         do igc = 1,nv
           write(31) w(igc),D(igc)
         end do
         close(31)
        end if
      end do
      xx(n) = xx(n-1)
      deallocate(D,st1,st2,w,wj)
!
!  Ordering the stations according to the clusters with LK being used for labeling the stations  
      allocate(lk(n))
      lk(1) = lkj(n-1)
      lk(2) = lki(n-1)
      do ii = 2, n-1
        do kk = 1, ii
          IF (lki(n-ii) == lk(kk)) exit
        end do
        do ik = kk,ii
          ki = ii+kk-ik
          lk(ki+1) = lk(ki)
        end do
        lk(kk) = lkj(n-ii)
      end do
      do  i = 1, n-1
        do ji = 1, n
          if (lk(ji) == lki(i)) exit
        end do
        do JJ = 1, n
          if (lk(jj) == lkj(i)) exit
        end do
        if (ji > n) ji = n
        if (jj > n) jj = n
        lki(i) = ji
        lkj(i) = jj
      end do
      do i = 2, n-2
        do j = 1, i-1
          if (lkj(i) /= lki(i-j)) cycle
          lkj(i) = lkj(i-j)
          exit
        end do
      end do
!
!   ### save data
      fnm = trim(outFNm)//trim(nstS)//'gc_clusters.bin'
      open(unit=31,file=trim(fnm),form='unformatted',action='write')
      do igc = 1,n
        write(31) lkj(igc),lki(igc),lk(igc),xx(igc)
      end do
      close(31)
!
     ! write(6,*) 'lki',lki
     ! write(6,*) 'lkj',lkj
      ! write(6,*) 'lk',lk
      write(*,*) 'lkj, lki, lk, 1-R'
      do igc = 1,n
         write(*,*) lkj(igc),lki(igc),lk(igc),xx(igc)
      end do
     deallocate(lkj,lki,lk,xx)
      
     return

      end subroutine clustering_ave
!
!#######################################################################
!
      subroutine compute_centroids(n,nv,iclus1,iclus2,outFNm,nstS)
    
       integer,intent(in)                    :: n,nv,iclus1,iclus2 
       character(len=256), intent(in)        :: outFNm
       character(len=10),intent(in)          :: nstS

       !local
       integer(kind=4)                   :: iclus,jclus,i,ii,j,jj,ij, &
                                            jpl1,jpl2,ncr
       integer(kind=4),allocatable, &
                       dimension(:)      :: ip,wd      ! tracking 
       integer(kind=4),allocatable, &
                       dimension(:)      :: lk,lki,lkj ! tracking clusters
       integer,allocatable,&
               dimension(:,:)            :: cluster_table  !order of clustering
       real(kind=4),allocatable, &
                       dimension(:)      :: xx         ! metric
       character(len=256)                :: fnm
!       
! Get data needed for computing centroids
       allocate(lkj(n),lki(n),lk(n),xx(n))
       fnm = trim(outFNm)//trim(nstS)//'gc_clusters.bin'
       open(unit=31,file=trim(fnm),form='unformatted',action='read')
       do igc = 1,n
         read(31) lkj(igc),lki(igc),lk(igc),xx(igc)
       end do
       close(31)
!
       allocate(wd(nv))
       wd = 1 ! false
       allocate(cluster_table(n,n))
       cluster_table = 0
       do i = 1, n-2
!        !$OMP PARALLEL PRIVATE(ij,ii,jpl1,jpl2,j)
         do ij = lkj(i),lki(i)-1
          !$OMP PARALLEL DO PRIVATE(ii,jpl1,jpl2,j)
           do ii = ij+1,lki(i)
             jpl1 = min(lk(ii),lk(ij))
             jpl2 = max(lk(ii),lk(ij))
             !j = map(n,jpl1,jpl2)  !map(n,i,j)=j + (i-1)*n - (i*(i+1))/2 , i<j
             j = jpl2 + (jpl1-1)*n - (jpl1*(jpl1+1))/2
             wd(j) = 0
           end do
         !$OMP END PARALLEL DO
         end do
!       !$OMP END PARALLEL
         if (i >= n-iclus2 .or. i <= n-iclus1) then
           allocate(ip(n))
           ip = 0
           ip(1) = 1
           ncr = 1
           ii = 1
           jj = 2
           call cluster_analysis(jj,n,nv,ncr,ii,i,wd,ip)
          !$OMP PARALLEL DO PRIVATE(j)
           do j = 1,n
             cluster_table(n-i,j) = ip(j)
           end do         
          !$OMP END PARALLEL DO
          deallocate(ip)
         end if 
       end do
!
!   ### save data
       fnm = trim(outFNm)//trim(nstS)//'gc_centroids.bin'
       open(unit=31,file=trim(fnm),form='unformatted',action='write')
       do i = 1,n
         do j = 1,n
           write(31) cluster_table(j,i)
         end do
       end do
       close(31)

       return

      end subroutine compute_centroids 
!
!************************************************************************
!
      recursive subroutine cluster_analysis(j,n,nv,ncr,ii,iipx1,wd,ip)

        integer,intent(in)                   :: nv,n,iipx1 
        integer,intent(in)                   :: wd(nv)     ! tracking clustering
        integer,intent(inout)                :: j,ii,ncr
        integer,dimension(n),intent(inout)   :: ip           

        !local
        integer                 :: jj,ij,jpl1,jpl2  ! counters, temp vars
!
        do jj = j, n
          jpl1 = min(ii,jj)
          jpl2 = max(ii,jj)
          !ij = map(n,jpl1,jpl2)   !map(n,i,j)=j + (i-1)*n - (i*(i+1))/2 , i<j
          ij = jpl2 + (jpl1-1)*n - (jpl1*(jpl1+1))/2
          if (wd(ij) == 1) cycle
          ip(jj) = ncr
          exit
        end do  
        ii = jj
        j = ii+1
        do while (ii >= n)
          ncr = ncr+1
          if (ncr > n-iipx1) return
          do JJ = 1, N
            if (IP(JJ) == 0) exit 
          end do
          if (jj > n) exit
          ip(jj) = ncr
          ii = jj
          j = ii+1
        end do
!
        call cluster_analysis(j,n,nv,ncr,ii,iipx1,wd,ip)
        
        return
        
      end subroutine cluster_analysis

!
!##################################################################################
!
       subroutine kzf(s,nt,m,n,w0,i1,kz,freq,period)
!  Calculates KZ(m,n)
         implicit none
         integer(kind=4), intent(in) :: nt,m,n, w0,i1
         real(kind=4), dimension(nt), intent(in) :: s
         real(kind=4), dimension(nt), intent(out) :: kz
         real(kind=4), intent(out) :: freq,period
!  local arrays, variables
         integer(kind=4) :: j,k, iter,ien,ist,jb,ju
         real(kind=4), dimension(nt) :: tmp,wrkz
         real(kind=4) :: sum, pi
!
         pi = acos(-1.0)
! Calculate frequency and period of filter in hours-1 and hours, respectively
         freq = sqrt(6.0)/pi * sqrt((1.0 - 0.5**(1.0/(2.*real(n))))/ &
                (real(m*m) - 0.5**(1.0/(2.*real(n)))))
         period = 1.0 / freq
!
! Fill work array with first pass values
!
         do j = 1, nt
           wrkz(j)   = s(j)
         end do
!
! Carry out iterations
!  Note that the span of the iterations decreases with each 
!  iteration, in order to allow the KZ filter to work on a 
!  previous running average. i1 is the desired starting hour
!  of the final output array, and w0 is the number of elements 
!  in that array.
!
!  First, do a check on the required bounds of the input data
!  in order to prevent attempts to read beyond the bounds in
!  creating the running averages.
! 
!
        if(i1 - n*(m-1)/2 < 1) then
          write(6,*) "Requested values of m and n for KZ are larger"
          write(6,*) "than can be accomodated with the input data array",&
                     " length included here"
          write(6,*) "Starting location requested in the data array", i1
          write(6,*) "Values of m and n requested were : ",m,n
          write(6,*) "Resulting initial array element in the data array"
          write(6,*) "that would need to be accessed would be :", &
                      (i1 - n*(m-1)/2)," which is out of bounds."
          write(6,*) "Stopping code..."
           stop
        end if
        if(i1 + w0 - 1 + n * (m-1) /2 > nt) then
          write(6,*) "Requested values of m and n for KZmn are larger"
          write(6,*) "than can be accomodated with the input data array",&
                     "  length included here"
          write(6,*) "Starting location requested in the data array: ",i1
          write(6,*) "Requested length of final array", w0
          write(6,*) "Values of m and n requested were : ",m,n
          write(6,*) "Resulting final array element in the data array"
          write(6,*) "that would  need to be accessed would be :", &
                      (i1+w0-1+n*(m-1)/2)," which is out of bounds"
          write(6,*) " due to maximum array length ",nt
          write(6,*) "Stopping code..."
          stop
        end if
!  Continue with the KZ filter
         do iter = n,1,-1
!  First, fill the outside-of-original-array values with values from the main array, assuming periodicity
! Start and end time index values of current iteration:
           ist = i1 - (iter - 1)*(m-1)/2
           ien = i1 + W0 -1 + (iter - 1)*(m - 1) / 2
!
!
!      write(6,*) 'KZ(',m,n,') iteration: ',iter
           do k = ist, ien
! Start and end intex values of the current running average:
             jb = k - (m - 1)/2
             ju = k + (m - 1)/2
             sum = 0.
             do j = jb,ju
                sum = sum + wrkz(j)
             end do
             tmp(k)=sum/real(m)
           end do
! Update the part of the work array for which the running average has been created.
! Not that the values of ist and ien will be smaller on the next iteration.I
           do k = ist,ien
              wrkz(k) = tmp(k)
           end do
         end do
!
!  Assign to output array
!
         do k = 1,w0
           kz(k) = wrkz(k+i1-1)
         end do
!
       return
       end subroutine kzf

!
!#######################################################################
!
      subroutine write_output(n,outFNmDataIn,outFNm,&
                              varbg_name, &
                              !mtype,diss,dirM,startDate, &
                              mtype,nstOpt,dirM,startDate, &
                              lat_mod,lon_mod,ip1,fhr)

       integer,intent(in)                   :: n,mtype,startDate,ip1,nstOpt
       character(len=2),intent(in)          :: fhr
       character(len=256),intent(in)        :: outFNmDataIn,outFNm,dirM
       character(len=4),intent(in)          :: varbg_name
       real(kind=4),intent(in)              :: lat_mod,lon_mod
       !real,intent(in)                     :: diss

       !local
       integer                        :: key,id_grid,&
                                         ni_grid,nj_grid,nk_grid, &
                                         ni_grid2,nj_grid2,nk_grid2
       integer                        :: swa_grid,lng_grid, dltf_grid,&
                                         ubc_grid,e1_grid,e2_grid,e3_grid
       integer                        :: dateo_grid,datyp_grid, &
                                         deet_grid,npas_grid,nbits_grid,&
                                         dateo_grid2,datyp_grid2, &
                                         deet_grid2,npas_grid2,nbits_grid2
       integer                        :: ip1_grid,ip2_grid,ip3_grid,ip4_grid, &
                                         ip1_grid2,ip2_grid2,ip3_grid2,ip4_grid2, &
                                         ig1_grid,ig2_grid,ig3_grid,ig4_grid, &
                                         ig1_grid2, ig2_grid2, ig3_grid2, ig4_grid2
       character(len=2)               :: grtyp_grid,typvar_grid, &
                                         grtyp_grid2,typvar_grid2
       character(len=12)              :: etiket_grid,etiket_grid2

       real(kind=4)                   :: istart,jstart, iend, jend, &
                                         i_stack, j_stack, &
                                         lat_stack,lon_stack
       integer(kind=4),allocatable, &
                       dimension(:)   :: lki,lkj,lk
       integer(kind=4),allocatable, &
                       dimension(:,:) :: cluster_table,cluster_matrix
       real(kind=4),allocatable, &
                       dimension(:)   :: xx
       real(kind=4),allocatable, &
                       dimension(:,:) :: rr_matrix,rr_work

       integer                              :: ier,fu_bg,fu_out,ii,jj,&
                                               nCells, i, j
       character(len=256)                   :: outFst
       character(len=7)                     :: mtypeStr
       character(len=10)                    :: dissStr,nstS

       integer                        :: ncl,jcl, kk,icl,iicl,jjcl,il,&
                                         iclus
       real                           :: rcold,rcoldi,rcoldj,dx,x1,x2,x3,x4,ykey,&
                                         xt,x0,y0,y2,y3
       logical                        :: moved
       character(len=3)               :: dig3
       character(len=4)               :: bggrid_def_field, nomvar_grid
       character(len=6)               :: tstep
       character(len=8)               :: model_date
       character(len=256)             :: obsfilename,dendrofilename,&
                                         dendronamefile,stationmerge,&
                                         airshedlocation,typelocation,&
                                         bg_filename
       character(len=60)              :: infoline
       character(len=60),dimension(n) :: dataname
       integer(kind=4),allocatable, &
               dimension(:)           :: station_id, &
                                         in_list,ki,kj,cm,cmnew, &
                                         alstatid,statid2, statid3,clt,clt2
       !integer,dimension(1,n)         :: ind, indnew                    
       integer,dimension(n,n)         :: ind, indnew
       real(kind=4),allocatable, &
               dimension(:)           :: longitude,latitude, &
                                         xx2, xx3 ,rc,rcnew
     
       logical                        :: KI_NOTINLIST,kj_notinlist
       real(kind=4),dimension(n)  :: longit2,latitu2, longit3,latitu3
       real(kind=4),allocatable,&
                     dimension(:) :: lon_bg, lat_bg

!  External functions
       integer  fstecr, gdll, fstinf, fstprm, fstlir, ezqkdef, gdrls, newdate, &
                ezsetopt, fnom, fstouv, fstfrm, fclos, gdxyfll, gdllsval, gdllfxy, fstopc
       external  fstec, gdll, fstinf, fstprm, fstlir,nezqkdef, gdrls, newdate, &
                ezsetopt, fnom, fstouv, fstfrm, fclos, gdxyfll, gdllsval, gdllfxy, fstopc

!
      allocate(station_id(nst),latitude(nst),longitude(nst))
      allocate(statid2(nst),statid3(nst), clt(nst), clt2(nst))
      open(unit=31,file=trim(outFNmDataIn),form='unformatted', &
                                           action='write')
      do igc = 1,nst
        write(31) station_id(igc),latitude(igc),longitude(igc)
      end do
      close(31)
      !
      write(nstS,'(i10)') nst
      nstS = adjustl(nstS)
     allocate(lkj(nst),lki(nst),lk(nst),xx(nst))
     allocate(in_list(nst),ki(nst),kj(nst))
     allocate(xx2(nst), xx3(nst))
     outFNmData = trim(outFNm)//trim(nstS)//'gc_clusters.bin'
     open(unit=31,file=trim(outFNmData),form='unformatted', &
                                                 action='read')
              do igc = 1,nst
                read(31) lkj(igc),lki(igc),lk(igc),xx(igc)
              end do
              close(31)
!
              allocate(cluster_table(nst,nst))
              outFNmData = trim(outFNm)//trim(nstS)//'gc_centroids.bin'
              open(unit=31,file=trim(outFNmData),form='unformatted', &
                                                 action='read')
              do i = 1,nst
                do j = 1,nst
                  read(31) cluster_table(j,i)
                end do
              end do
              close(31)

! Check if sorted locations lki(i) and lkj(i) have appeared already in the list
      ncl = 2
      in_list(1) = lki(1)
      in_list(2) = lkj(1)
      ki(1) = -1
      kj(1) = -1
      do i = 2,n-1
         ki_notinlist = .true.
         kj_notinlist = .true.
         do ii = 1,ncl
           if(in_list(ii) == lki(i)) ki_notinlist = .false.
           if(in_list(ii) == lkj(i)) kj_notinlist = .false.
         end do
         if(ki_notinlist) then
            ncl = ncl + 1
            in_list(ncl) = lki(i)
            ki(i) = -1
         else
            ki(i) = 1
         end if
         if(kj_notinlist) then
            ncl = ncl + 1
            in_list(ncl) = lkj(i)
            kj(i) = -1
         else
            kj(i) = 1
         end if
      end do
!
     !$OMP PARALLEL DO PRIVATE(i)
      do i = 1,n-1
         ki(i) = lki(i)*ki(i)
         kj(i) = lkj(i)*kj(i)
      end do
     !$OMP END PARALLEL DO
!
! Getting clusters based on dissimilarity or number of clusters 
      iclus=0
      !do i = n,1,-1
      !  iclus=iclus+1
      !  if (xx(i) < diss) exit
      !end do
      do i = 1,n
        if (maxval(cluster_table(i,:))>nstOpt) exit
        iclus=iclus+1
      end do
      if (iclus == 0) then
        write(6,*) 'WARNING: no. of clusters requested is too high'
        write(6,*) 'WARNING: no. of clusters set to 1'
        iclus=1
      end if
!
!  The first time ki or kj is negative corresponds to the dissimilarity when that
!  station merges with another cluster.  Determine the numerical R values of the 
!  stations when this happens:
      il = 0
      do i = 1,n-1
        if(ki(i) < 0) then
          il = il + 1
          ii = lk(-ki(i))
          xx2(il) = xx(i)
          statid2(il) = station_id(ii)
          longit2(il) = longitude(ii)
          latitu2(il) = latitude(ii)
          if (iclus /= 1 .and. iclus /= n) then
            clt(il) =  cluster_table(iclus,ii) ! cluster_table order stations from 1:n
          else
            if (iclus == 1) clt(il) = 1
            if (iclus == n) clt(il) = il 
         end if           
        end if
        if(kj(i) < 0) then
          il = il + 1
          ii = lk(-kj(i))
          xx2(il) = xx(i)
          statid2(il) = station_id(ii)
          longit2(il) = longitude(ii)
          latitu2(il) = latitude(ii)
          if (iclus /= 1 .and. iclus /= n) then
            clt(il) = cluster_table(iclus,ii)
          else 
            if (iclus == 1) clt(il) = 1
            if (iclus == n) clt(il) = il
          end if
        end if
      end do
!
! order stations for fst file
     !$OMP PARALLEL PRIVATE(i,j)
      do i = 1,n
        do j = 1,n
          if (statid2(j) /= i) cycle
          statid3(i) = statid2(j)
          longit3(i) = longit2(j)
          latitu3(i) = latitu2(j)
          xx3(i) =  xx2(j)
          clt2(i) =  clt(j)
          exit
        end do 
      end do
     !$OMP END PARALLEL 
!
!  write binary files 
!
      fu_out = 100
      outfilename = TRIM(outFNm)//'cluster_table.bin'
      open(unit=fu_out,file=outfilename,form='unformatted')
      write(fu_out) cluster_table
      close(fu_out)
!
! colect data for the fst
!
!! infomation abot the model domain for the fst files 
      !   Turn off FST information messages
      ier = fstopc('MSGLVL', 'SYSTEM', 0)
      fu_bg = 15
      !!! get all the information needed from 1st time step
      write(model_date, '(i8.8)') startDate
      write(tstep, '(i6.6)') 1 * 60
      bg_filename = trim(dirM)//model_date//trim(fhr)//&
                    '/'// model_date//trim(fhr)//'_'//&
                    trim(tstep)//'p'
!		    
! Opening input files
      ier = fnom(fu_bg, trim(bg_filename), 'RND+OLD+R/O', 0)
      if (ier < 0) print *, 'error while opening ', trim(bg_filename)
      ier = fstouv(fu_bg, 'RND')
      nomvar_grid = trim(varbg_name)
      key = fstinf(fu_bg,ni_grid,nj_grid, nk_grid, -1, ' ',&
                               -1, -1, -1, ' ',nomvar_grid)
      if (key < 0) then
        write(6,*) 'Key value negative after opening file',&
                   ' with grid definitions'
       stop
      end if
!  Getting all the description information of the record given the handle(key)
      ier = fstprm(key, dateo_grid, deet_grid, npas_grid,     &
                   ni_grid, nj_grid, nk_grid, nbits_grid,     &
                   datyp_grid, ip1_grid, ip2_grid, ip3_grid,  &
                   typvar_grid, nomvar_grid, etiket_grid,     &
                   grtyp_grid, &
                   ig1_grid, ig2_grid, ig3_grid, ig4_grid,    &
                   swa_grid, lng_grid, dltf_grid, ubc_grid,   &
                   e1_grid, e2_grid, e3_grid)
      if (ier < 0) then
        write(6,*) 'ier value negative after attempting to ', &
                   'read grid parameters'
        stop
      end if
      ier = ezsetopt('VERBOSE', 'NO')
      id_grid = ezqkdef(ni_grid,nj_grid,grtyp_grid,ig1_grid,&
                        ig2_grid,ig3_grid,ig4_grid,fu_bg)
      allocate(rr_matrix(ni_grid,nj_grid))
      allocate(cluster_matrix(ni_grid,nj_grid))
      rr_matrix = 0
      cluster_matrix = 0
      nCells = floor(sqrt(real(n))*0.5)
      ier = gdxyfll(id_grid, i_stack, j_stack, lat_mod,lon_mod, 1)
      i=0
      do jj = -nCells,nCells
        do ii = -nCells,nCells
          i=i+1
          rr_matrix(i_stack+ii,j_stack+jj) = xx3(i)
          cluster_matrix(i_stack+ii,j_stack+jj) = clt2(i)
        end do
      end do
!
! write fst files
!
      write(nstS,'(i10.10)') n
      outFst = trim(outFnm)//'_'//trim(nstS)//'.fst'
      fu_out = 100
      if (mtype == 1) mtypeStr = 'EuD'
      if (mtype == 2) mtypeStr = 'EuDN'
      if (mtype == 3) mtypeStr = 'R'
      if (mtype == 4) mtypeStr = '1-R'
      if (mtype == 5) mtypeStr = '1-RxEuD'
      if (mtype == 6) mtypeStr = 'variance'
     ! write(dissStr,'(a5,f5.2)') 'diss_',diss
      write(dissStr,'(a5,f5.2)') 'diss_',xx(n-iclus)
      allocate(rr_work(ni_grid,nj_grid))
      call em_open_fstfile(fu_out, trim(outFst),'STD+RND', 'RND')
!  Locating the next record that matches the search keys (datev,etiket,ip1,ip2,ip3,typvar,nomvar) 
!  Reading grid parameters
      key = fstinf(fu_bg, ni_grid2, nj_grid2, nk_grid2, -1,  &
                   ' ', 36776, 89595, -1, 'X', '>>')
      if (key < 0) then
        write(6,*) 'Key value for tic is negative after',&
                   ' opening file with grid definitions'
        stop
      end if
      allocate(lon_bg(ni_grid2))
      key = fstinf(fu_bg, ni_grid2, nj_grid2, nk_grid2, -1,  &
                   ' ', 36776, 89595, -1, 'X', '^^')
      if (key < 0) then
        write(6,*) 'Key value for  tac is negative after',&
                   ' opening file with grid definitions'
        stop
      end if
      allocate(lat_bg(nj_grid2))
! Geting all the description information of the record given the handle 
      ier = fstprm(key, dateo_grid2, deet_grid2, npas_grid2,     &
                   ni_grid2, nj_grid2, nk_grid2, nbits_grid2,     &
                   datyp_grid2, ip1_grid2, ip2_grid2, ip3_grid2,  &
                   typvar_grid2, nomvar_grid, etiket_grid2,     &
                   grtyp_grid2, &
                   ig1_grid2, ig2_grid2, ig3_grid2, ig4_grid2,    &
                   swa_grid, lng_grid, dltf_grid, ubc_grid,   &
                   e1_grid, e2_grid, e3_grid)
! Locating and reading the next record that matches the search keys (datev,etiket,ip1,ip2,ip3,typvar,nomvar). 
      ier = fstlir(lon_bg, fu_bg,  ni_grid2, nj_grid2, nk_grid2, &
                   -1, etiket_grid2, ip1_grid2, ip2_grid2, -1,&
                   'X', '>>')
      if (ier < 0) then
        write(6,*) 'ier value negative after attempting to ', &
                   'set output grid parameters^^'
        stop
      end if
!  Storing lat and lon fields
      ier = fstecr(lon_bg, rr_work, -32, fu_out, dateo_grid2, 0, 0, &
                   ni_grid2, 1, 1, ip1_grid2, ip2_grid2, ip3_grid2, &
                   'X', '>>', etiket_grid2, grtyp_grid2, &
                   ig1_grid2, ig2_grid2, ig3_grid2, ig4_grid2, &
                   datyp_grid2, .true.)
!
      ier = fstlir(lat_bg, fu_bg, ni_grid2, nj_grid2, nk_grid2, &
                   -1, etiket_grid2, ip1_grid2, ip2_grid2, -1,&
                   'X', '^^')
      if (ier < 0) then
        write(6,*) 'ier value negative after attempting to ', &
                   'set output grid parameters:>>'
        stop
      end if
      ier = fstecr(lat_bg, rr_work, -32, fu_out, dateo_grid2, 0, 0, &
                   1,nj_grid2, 1, ip1_grid2, ip2_grid2, ip3_grid2, &
                   'X', '^^', etiket_grid2, grtyp_grid2, &
                   ig1_grid2, ig2_grid2, ig3_grid2, ig4_grid2, &
                   datyp_grid2, .true.)
!  Storing metric and clusters fields
      ier = fstecr(rr_matrix, rr_work, -nbits_grid, fu_out,&
                   dateo_grid, deet_grid, npas_grid, &
                   ni_grid,nj_grid, 1, ip1_grid, ip2_grid, ip3_grid, &
                   typvar_grid, trim(mtypeStr), etiket_grid, grtyp_grid, &
                   ig1_grid, ig2_grid, ig3_grid, ig4_grid, datyp_grid, .true.)
      ier = fstecr(cluster_matrix, rr_work, -nbits_grid, fu_out,&
                   dateo_grid, deet_grid, npas_grid, &
                   ni_grid,nj_grid, 1, ip1_grid, ip2_grid, ip3_grid, &
                   typvar_grid, trim(dissStr), etiket_grid, grtyp_grid, &
                   ig1_grid, ig2_grid, ig3_grid, ig4_grid, datyp_grid, .true.)
!
      call em_close_fstfile(fu_out)
      call em_close_fstfile(fu_bg)
      
      deallocate(lon_bg,lat_bg)
      deallocate(rr_matrix,cluster_matrix,rr_work)

      return

      end subroutine write_output
!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sort(n,arr)

       integer,intent(in)              :: n
       integer,dimension(n),intent(inout) :: arr

       ! local
       integer(kind=4),dimension(n) :: indx

       call indexx(n,arr,indx)

       arr = arr(indx)

       return

      end subroutine sort
!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sortr(n,arr)

       integer,intent(in)              :: n
       real,dimension(n),intent(inout) :: arr

       ! local
       integer(kind=4),dimension(n) :: indx

       call indexxr(n,arr,indx)

       arr = arr(indx)

       return

      end subroutine sortr

!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sort2(n,arr,arrTmp)

       integer,intent(in)              :: n
       integer,dimension(n),intent(inout) :: arr,arrTmp

       ! local
       integer(kind=4),dimension(n) :: indx

       call indexx(n,arr,indx)

       arr = arr(indx)
       arrTmp = arrTmp(indx)

       return

      end subroutine sort2 
!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sort2_r_i(n,arr1,arr2)

       integer,intent(in)              :: n
       real,dimension(n),intent(inout) :: arr1
       integer,dimension(n),intent(inout) :: arr2

       ! local
       integer(kind=4),dimension(n) :: indx

       call indexxr(n,arr1,indx)

       arr1 = arr1(indx)
       arr2 = arr2(indx)

       return

      end subroutine sort2_r_i
!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sort3(n,arr,arrTmp1,arrTmp2)
  
       integer,intent(in)                 :: n
       integer,dimension(n),intent(inout) :: arrTmp1,arrTmp2
       integer,dimension(n),intent(inout)    :: arr
   
       ! local
       integer(kind=4),dimension(n) :: indx

       call indexx(n,arr,indx)
       
       arr = arr(indx)
       arrTmp1 = arrTmp1(indx)
       arrTmp2 = arrTmp2(indx)

       return

      end subroutine sort3
!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sort3_2i_1r(n,arr,arr1,arr2)

       integer,intent(in)                 :: n
       integer,dimension(n),intent(inout) :: arr,arr1
       real,dimension(n),intent(inout)    :: arr2

       ! local
       integer(kind=4),dimension(n) :: indx

       call indexx(n,arr,indx)

       arr = arr(indx)
       arr1 = arr1(indx)
       arr2 = arr2(indx)

       return

      end subroutine sort3_2i_1r

!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sortr_2i(n,arr,arr1,arr2)

       integer,intent(in)                 :: n
       integer,dimension(n),intent(inout) :: arr1,arr2
       real,dimension(n),intent(inout)    :: arr

       ! local
       integer(kind=4),dimension(n) :: indx

       call indexxr(n,arr,indx)

       arr = arr(indx)
       arr1 = arr1(indx)
       arr2 = arr2(indx)

       return

      end subroutine sortr_2i
!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sort4(n,arr,arrTmp1,arrTmp2,arrTmp3)

       integer,intent(in)                 :: n
       integer,dimension(n),intent(inout) :: arrTmp1,arrTmp2,arrTmp3
       integer,dimension(n),intent(inout)    :: arr

       ! local
       integer(kind=4),dimension(n) :: indx

       call indexx(n,arr,indx)

       arr = arr(indx)
       arrTmp1 = arrTmp1(indx)
       arrTmp2 = arrTmp2(indx)
       arrTmp3 = arrTmp3(indx)

       return

      end subroutine sort4
!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sort4_3i_1r(n,arr,arr1,arr2,arr3)

       integer,intent(in)                 :: n
       integer,dimension(n),intent(inout) :: arr,arr1,arr2
       real,dimension(n),intent(inout)    :: arr3

       ! local
       integer(kind=4),allocatable, &
                       dimension(:)    :: indx

       allocate(indx(n))
       call indexx(n,arr,indx)
       arr = arr(indx)
       arr1 = arr1(indx)
       arr2 = arr2(indx)
       arr3 = arr3(indx)
       deallocate(indx)

       return

      end subroutine sort4_3i_1r
!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine sort4_1r_3i(n,arr,arr1,arr2,arr3)

       integer,intent(in)                 :: n
       real,dimension(n),intent(inout)    :: arr
       integer(kind=4),dimension(n),intent(inout) :: arr3,arr1,arr2

       ! local
       integer,dimension(n) :: indx

       call indexxr(n,arr,indx)

       arr = arr(indx)
       arr1 = arr1(indx)
       arr2 = arr2(indx)
       arr3 = arr3(indx)

       return

      end subroutine sort4_1r_3i

!
!#######################################################################
!  Based on the numerical recipies for F90

      subroutine indexx(n,arr,indx)
       integer(kind=4),intent(in)               :: n
       integer(kind=4),dimension(n),intent(in)     :: arr
       integer(kind=4),dimension(n),intent(out)  :: indx
      
       ! local
       integer, parameter        :: nn = 15,nstack = 50 
       integer(kind=4)           :: indext,varTmp,a
       integer(kind=4)           :: i,ir,itemp,j,jstack,k, &
                                    l,r
       integer,dimension(nstack) :: istack
 
      indx=arth_i(1,1,n)
      jstack=0
      l=1
      r=n
      do
      if (r-l < NN) then
        do j=l+1,r
          indext=indx(j)
          a=arr(indext)
          do i=j-1,l,-1
            if (arr(indx(i)) <= a) exit
            indx(i+1)=indx(i)
          end do
          indx(i+1)=indext
        end do
        if (jstack == 0) RETURN
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+r)/2
        ! swap
        varTmp = indx(k)
        indx(k) = indx(l+1)
        indx(l+1) = varTmp
        if (arr(indx(r)) < arr(indx(l))) then
          varTmp=indx(l)
          indx(l) = indx(r)
          indx(r) = varTmp
        end if
        if (arr(indx(r)) < arr(indx(l+1))) then
          varTmp=indx(l+1)
          indx(l+1) = indx(r)
          indx(r) = varTmp
        end if
        if (arr(indx(l+1)) < arr(indx(l))) then
          varTmp=indx(l)
          indx(l) = indx(l+1)
          indx(l+1) = varTmp
        end if
        
        i=l+1
        j=r
        indext=indx(l+1)
        a=arr(indext)
        do
          do
            i=i+1
            if (arr(indx(i)) >= a) exit
          end do
          do
            j=j-1
            if (arr(indx(j)) <= a) exit
          end do
          if (j < i) exit
          !swap
          varTmp = indx(i)
          indx(i) = indx(j)
          indx(j) = varTmp
        end do
        indx(l+1)=indx(j)
        indx(j)=indext
        jstack=jstack+2
        if (jstack > NSTACK) write(6,*) ('indexx:NSTACK too small')
        if (r-i+1 >= j-l) then
          istack(jstack)=r
          istack(jstack-1)=i
          r=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        end if
      end if
      end do
 
      return

      end subroutine indexx 
!
! ###############################################
!
      subroutine indexxr(n,arr,indx)
       integer,intent(in)               :: n
       real,dimension(n),intent(in)     :: arr
       integer(kind=4),dimension(n),intent(out)  :: indx

       ! local
       integer, parameter        :: nn = 15,nstack = 50
       integer                   :: i,indext,ir,itemp,j,jstack,k, &
                                    l,varTmp,r
       integer,dimension(nstack) :: istack
       real                      :: a
!
!
      indx=arth_i(1,1,n)
      jstack=0
      l=1
      r=n
      do
      if (r-l < NN) then
        do j=l+1,r
          indext=indx(j)
          a=arr(indext)
          do i=j-1,l,-1
            if (arr(indx(i)) <= a) exit
            indx(i+1)=indx(i)
          end do
          indx(i+1)=indext
        end do
        if (jstack == 0) RETURN
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+r)/2
        ! swap
        varTmp = indx(k)
        indx(k) = indx(l+1)
        indx(l+1) = varTmp
        if (arr(indx(r)) < arr(indx(l))) then
          varTmp=indx(l)
          indx(l) = indx(r)
          indx(r) = varTmp
        end if
        if (arr(indx(r)) < arr(indx(l+1))) then
          varTmp=indx(l+1)
          indx(l+1) = indx(r)
          indx(r) = varTmp
        end if
        if (arr(indx(l+1)) < arr(indx(l))) then
          varTmp=indx(l)
          indx(l) = indx(l+1)
          indx(l+1) = varTmp
        end if

        i=l+1
        j=r
        indext=indx(l+1)
        a=arr(indext)
        do
          do
            i=i+1
            if (arr(indx(i)) >= a) exit
          end do
          do
            j=j-1
            if (arr(indx(j)) <= a) exit
          end do
          if (j < i) exit
          !swap
          varTmp = indx(i)
          indx(i) = indx(j)
          indx(j) = varTmp
        end do
        indx(l+1)=indx(j)
        indx(j)=indext
        jstack=jstack+2
        if (jstack > NSTACK) write(6,*) ('indexx:NSTACK too small')
        if (r-i+1 >= j-l) then
          istack(jstack)=r
          istack(jstack-1)=i
          r=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        end if
      end if
      end do

      return

      end subroutine indexxr

!
!######################################################################
!
      function arth_i(first,increment,n)
        
        integer,intent(in) :: first,increment,n

        ! local
        integer              :: k,k2,temp
        integer,dimension(n) :: arth_i
        integer,parameter :: NPAR_ARTH = 16, NPAR2_ARTH = 8

        if (n > 0) arth_i(1)=first
        if (n <= NPAR_ARTH) then
          do k=2,n
            arth_i(k)=arth_i(k-1)+increment
          end do
        else
          do k=2,NPAR2_ARTH
            arth_i(k)=arth_i(k-1)+increment
          end do
          temp=increment*NPAR2_ARTH
          k=NPAR2_ARTH
          do
            if (k >= n) exit
            k2=k+k
            arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
          end do
        end if
 
       end function arth_i
!
!######################################################################
!
      function map(n,i,j)
       !map row i and column j of upper half diagonal symmetric matrix onto vector
        integer,intent(in) :: n,i,j
        integer :: map

        map = j + (i-1)*n - (i*(i+1))/2


      end function map
!
!######################################################################
!
      subroutine percentiler(x, length, per, oout)

      integer(kind=4),intent(in)         :: length    ! elements # of x and per
      real(kind=4),dimension(length), &
                   intent(in)           :: x   ! input array, and its non-MISSING part
      real(kind=4),intent(in)           :: per !,oout  ! percentile level and results
      real(kind=4),intent(out)          :: oout  ! percentile level and results

      ! local
      real(kind=4),allocatable, &
                   dimension(:)      ::xtos    ! input array, and its non-MISSING part
      real(kind=4)    :: bb,cc ! temporary variable
      integer :: nn,i  ! loop index and count #

     allocate(xtos(length))
      nn=0
      do i=1, length
          nn=nn+1
          xtos(nn)=x(i)
      enddo
      bb=nn*per+per/3.+1./3.
      cc=real(int(bb))
      if(int(cc).ge.nn) then
        oout=xtos(nn)
      else
        oout=xtos(int(cc))+(bb-cc)*       &
            (xtos(int(cc)+1)-xtos(int(cc)))
      endif

     return

      end subroutine percentiler

!
!######################################################################
!
! Calculate percentile
!
      !subroutine percentile(x, length, nl, per, oout)
      subroutine percentile(x, length, per, oout)
      !use COMM
      !use functions
      integer :: length,nl    ! elements # of x and per
      !real(kind=4),dimension(length):: x,xtos    ! input array, and its non-MISSING part
      integer(kind=4),dimension(length):: x,xtos    ! input array, and its non-MISSING part
      !real(kind=4),dimension(nl)    :: per,oout  ! per! percentile level and results
      real(kind=4)    :: per,oout  ! percentile level and results
      real(kind=4)    :: bb,cc ! temporary variable
      integer :: nn,i  ! loop index and count #

! check if percentile level is out of bound.
      !if(maxval(per)>1..or.minval(per)<0.) & 
      !   stop 'ERROR in percentile: per out of range [0.,1.] !'
! choose non-MISSING input, but is it necessary? I think x is already chosen ......
      nn=0
      do i=1, length
      !  if(nomiss(x(i)))then
          nn=nn+1
          xtos(nn)=x(i)
      !  endif
      enddo

      !if(nn.eq.0) then
      !  oout=MISSING
      !else
        !call sortr(nn,xtos) ! ascending numerical order
        call sort(nn,xtos) ! ascending numerical order
        !do i=1,nl
          !bb=nn*per(i)+per(i)/3.+1./3.
          bb=nn*per+per/3.+1./3.
          cc=real(int(bb))
          if(int(cc).ge.nn) then
            !oout(i)=xtos(nn)
            oout=xtos(nn)
          else
            !oout(i)=xtos(int(cc))+(bb-cc)*       &
            oout=xtos(int(cc))+(bb-cc)*       &
               (xtos(int(cc)+1)-xtos(int(cc)))
          endif
        !enddo
     ! endif

      return
      end subroutine percentile

!
!######################################################################
!
!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : em_open_fstfile.ftn90
! Creation       : H. Landry - octobre 2006
! Description    : Open an fst file
!




! Extra info     :
!
!! Arguments: IN
!              file_name -> speaks for itself
!              file_type -> RPN type value for fnom
!                          (usually something like 'RND+OLD')
!              options   -> RPN options for fstouv
!                          (usually something like 'RND')
!            IN/OUT
!              file_unit -> the file unit to open (an initial value of 0 will
!              assign a unit automatically - Recommended)
!
!==============================================================================
subroutine em_open_fstfile(FILE_UNIT, file_name, file_type, options)
   implicit none
!  Arguments
   integer, intent(inout) :: file_unit
   character*(*), intent(in) :: file_name
   character*(*), intent(in) :: file_type
   character*(*), intent(in) :: options
!  External functions
   integer  fnom, fstouv
   external fnom, fstouv
   if (fnom (file_unit, file_name, file_type, 0) < 0) then
      write(0, *) '### Error in em_open_fstfile ###'
      write(0, *) '# fnom failed for:', file_name
      write(0, *) '###           ABORT          ###'
      call qqexit(1)
   endif
   if (fstouv(file_unit, options) < 0) then
      write(0, *) '### Error in em_open_fstfile ###'
      write(0, *) '# fstouv failed for:', file_name
      write(0, *) '###           ABORT          ###'
      call qqexit(1)
   endif
   print *, '   Opened ', trim(file_name), ' with unit:', file_unit
end subroutine em_open_fstfile
                                                            

!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : em_close_fstfile.ftn90
! Creation       : H. Landry, Oct 2006
! Description    : Close an opened fst file
!
! Extra info     :
!
! Arguments:  IN
!               file_unit -> the file unit to close
!
!============================================================================
subroutine em_close_fstfile(file_unit)
   implicit none
   integer, intent(in) :: file_unit
   integer err
   integer  fclos, fstfrm
   external fclos, fstfrm
   err = fstfrm(file_unit)
   err = fclos(file_unit)
   print *, '   Closed unit:', file_unit
end subroutine em_close_fstfile

      end program kzfilter_hc
