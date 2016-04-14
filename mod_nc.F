      module mod_nc
      use mod_xc        ! HYCOM communication interface
      use mod_za        ! HYCOM I/O interface
      use mod_cb_arrays ! HYCOM saved arrays

      implicit none
      include 'netcdf.inc'

!---  Module to output Netcdf variable (debug purposes) ---


      contains

      subroutine ncput_2d(var_nc2d,var_name,file_name)
      use mod_xc        ! HYCOM communication interface
      use mod_cb_arrays ! HYCOM saved arrays

!--- To output 2d variables 
      implicit none      
      include 'netcdf.inc'

      real*8, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy )::
     &  var_nc2d
      character(len=*) :: file_name,var_name
      CHARACTER*(NF_MAX_NAME) recname

      integer ::  ncid,i,j,nrec
      integer  , parameter :: ndims=3
      integer  , parameter :: nrecs=5
      character(len=*), parameter :: lvl_name  = 'layer'
      character(len=*), parameter :: lat_name  = 'latitude'
      character(len=*), parameter :: lon_name  = 'longitude'
      character(len=*), parameter :: rec_name  = 'time'
      character(len=*), parameter :: units     = 'units'
      character(len=*), parameter :: lat_units = 'degrees_north'
      character(len=*), parameter :: lon_units = 'degrees_east'
      character(len=*), parameter :: time_units = 'step'
      integer :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid,
     &           lon_varid, lat_varid, var_varid, retval,rec_varid
      integer, dimension(ndims) :: dimids, start, count
      integer, dimension(2)     :: cdimids
      real*8 :: date
      logical :: wcoord
                     
! --------------------------------------------------------------
!            OPEN NETCDF FILE
! --------------------------------------------------------------
      if (mnproc.eq.1) then              
!---  Create the file 
          retval = nf_create(file_name, nf_noclobber, ncid)

!---  Check if file exists or not      
          if (retval .eq. nf_eexist) then
              retval = nf_open(file_name, nf_write, ncid)
              retval = nf_inq_varid(ncid,var_name,var_varid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              wcoord = .false.
          else
!---  Define the dimensions 
              retval = nf_def_dim(ncid,lat_name,jtdm,lat_dimid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_def_dim(ncid,lon_name,itdm,lon_dimid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_def_dim(ncid,rec_name,nf_unlimited,rec_dimid)
              if (retval .ne. nf_noerr) call handle_err(retval)

!---  Define the coordinate variables
              cdimids(1)=lon_dimid
              cdimids(2)=lat_dimid
!              cdimids(3)=rec_dimid
              retval = nf_def_var(ncid, lat_name, nf_real, 2, cdimids, 
     &            lat_varid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_def_var(ncid, lon_name, nf_real, 2, cdimids, 
     &            lon_varid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              
!---  Assign units attributes to coordinate variables
              retval = nf_put_att_text(ncid, lat_varid, units,
     &            len(lat_units), lat_units)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_put_att_text(ncid, lon_varid, units,
     &            len(lon_units), lon_units)
              
!---  Assign dimensions of the netCDF variables 
              dimids(1) = lon_dimid
              dimids(2) = lat_dimid
              dimids(3) = rec_dimid

!---  Define the netCDF variables var_nc2d.
              retval = nf_def_var(ncid, var_name,nf_double,ndims,dimids, 
     &            var_varid)
              if (retval .ne. nf_noerr) call handle_err(retval)

!---  End define mode
              retval = nf_enddef(ncid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              wcoord = .true.
          endif  ! .nf_eexist. 

! --- check time status
          retval = nf_inq_dimid (ncid, rec_name, rec_varid)
          if (retval .ne. nf_noerr) call handle_err(retval)
          retval = nf_inq_dimlen(ncid, rec_varid, nrec)
          if (retval .ne. nf_noerr) call handle_err(retval)
      endif                     ! mnproc.eq.1

!---  Write the coordinate lon/lat
      do j=1,jj
        do i=1,ii
          util1(i,j) = plat(i,j)
        enddo
      enddo
      call xcaget(worknc1, util1, 1) !! get value from all procs
      if (mnproc.eq.1) then
        if (wcoord) then
          retval = nf_put_var_double(ncid, lat_varid, worknc1)
          if (retval .ne. nf_noerr) call handle_err(retval)
        endif ! wcoord
      endif  ! mnproc
      do j=1,jj
        do i=1,ii
          util1(i,j) = plon(i,j)
        enddo
      enddo      
      call xcaget(worknc1, util1, 1) !! get value from all procs
      if (mnproc.eq.1) then
        if (wcoord) then
          retval = nf_put_var_double(ncid, lon_varid, worknc1)
          if (retval .ne. nf_noerr) call handle_err(retval)
        endif  !wcoord
      endif ! mnproc

!---  Write one timestep of data 
      if (mnproc.eq.1) then 
          count(1) = itdm
          count(2) = jtdm
          count(3) = 1
          start(1) = 1
          start(2) = 1
          start(3) = nrec+1
      endif ! mnproc

!---  Write var_nc2d
      do j=1,jj
        do i=1,ii
          util1(i,j) = var_nc2d(i,j)
        enddo
      enddo
      call xcaget(worknc1, util1, 1) !! get value from all procs
      if (mnproc.eq.1) then
          retval = nf_put_vara_double(ncid, var_varid, start, count, 
     &        worknc1(1:itdm,1:jtdm))
          if (retval .ne. nf_noerr) call handle_err(retval)
      endif !mnproc

!---  Close the file
      if (mnproc.eq.1) then
          retval = nf_close(ncid)
          if (retval .ne. nf_noerr) call handle_err(retval)
      endif ! mnproc
      
      end subroutine ncput_2d


      subroutine ncput_3d(var_nc3d,var_name,file_name)
      use mod_xc        ! HYCOM communication interface
      use mod_cb_arrays ! HYCOM saved arrays

!--- To output 3d variables 
      implicit none      
      include 'netcdf.inc'

      real*8,dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm),intent(in)::
     &  var_nc3d
      character(len=*), intent(in) :: file_name,var_name
      CHARACTER*(NF_MAX_NAME) recname

      integer ::  ncid,i,j,k,nrec
      integer  , parameter :: ndims=4
      integer  , parameter :: nrecs=5
      character(len=*), parameter :: lvl_name  = 'layer'
      character(len=*), parameter :: lat_name  = 'latitude'
      character(len=*), parameter :: lon_name  = 'longitude'
      character(len=*), parameter :: rec_name  = 'time'
      character(len=*), parameter :: units     = 'units'
      character(len=*), parameter :: lat_units = 'degrees_north'
      character(len=*), parameter :: lon_units = 'degrees_east'
      integer :: lvl_dimid, lon_dimid, lat_dimid, rec_dimid,
     &           lon_varid, lat_varid, var_varid, retval
      integer, dimension(ndims) :: dimids, start, count,rec_varid
      integer, dimension(2)     :: cdimids
      logical :: wcoord

! --------------------------------------------------------------
!            OPEN NETCDF FILE
! --------------------------------------------------------------
      if (mnproc.eq.1) then              
!---  Create the file 
          retval = nf_create(file_name, nf_noclobber, ncid)

!---  Check if file exists or not      
          if (retval .eq. nf_eexist) then
              retval = nf_open(file_name, nf_write, ncid)
              retval = nf_inq_varid(ncid,var_name,var_varid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              wcoord = .false.
          else
!---  Check if file exists or not      
!          if (retval .eq. nf_eexist) then
!              retval= nf_open(file_name, nf_write, ncid)
!          else

!---  Define the dimensions 
              retval = nf_def_dim(ncid,lvl_name,kdm ,lvl_dimid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_def_dim(ncid,lat_name,jtdm,lat_dimid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_def_dim(ncid,lon_name,itdm,lon_dimid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_def_dim(ncid,rec_name,nf_unlimited,rec_dimid)
              if (retval .ne. nf_noerr) call handle_err(retval)

!---  Define the coordinate variables
              cdimids(1)=lon_dimid
              cdimids(2)=lat_dimid
              retval = nf_def_var(ncid, lat_name, nf_real, 2, cdimids, 
     &            lat_varid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_def_var(ncid, lon_name, nf_real, 2, cdimids, 
     &            lon_varid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              
!---  Assign units attributes to coordinate variables
              retval = nf_put_att_text(ncid, lat_varid, units,
     &            len(lat_units), lat_units)
              if (retval .ne. nf_noerr) call handle_err(retval)
              retval = nf_put_att_text(ncid, lon_varid, units,
     &            len(lon_units), lon_units)
              if (retval .ne. nf_noerr) call handle_err(retval)
              
!---  Assign dimensions of the netCDF variables 
              dimids(1) = lon_dimid
              dimids(2) = lat_dimid
              dimids(3) = lvl_dimid
              dimids(4) = rec_dimid

!---  Define the netCDF variables var_nc2d.
              retval = nf_def_var(ncid,var_name,nf_double,ndims,dimids, 
     &            var_varid)
              if (retval .ne. nf_noerr) call handle_err(retval)

!---  End define mode
              retval = nf_enddef(ncid)
              if (retval .ne. nf_noerr) call handle_err(retval)
              wcoord = .true.
          endif  ! .nf_exist. 
! --- check time status
          retval = nf_inq_dimid (ncid, rec_name, rec_varid)
          if (retval .ne. nf_noerr) call handle_err(retval)
          retval = nf_inq_dimlen(ncid, rec_varid, nrec)
          if (retval .ne. nf_noerr) call handle_err(retval)
      endif  ! mnproc.eq.1

!---  Write the coordinate lon/lat 
      do j=1,jj
        do i=1,ii
          util1(i,j) = plat(i,j)
        enddo
      enddo
      call xcaget(worknc1, util1, 1) !! get value from all procs
      if (mnproc.eq.1) then
        if (wcoord) then
          retval = nf_put_var_double(ncid, lat_varid, worknc1)
          if (retval .ne. nf_noerr) call handle_err(retval)
        endif ! wcoord
      endif  ! mnproc
      do j=1,jj
        do i=1,ii
          util1(i,j) = plon(i,j)
        enddo
      enddo      
      call xcaget(worknc1, util1, 1) !! get value from all procs
      if (mnproc.eq.1) then
        if (wcoord) then
          retval = nf_put_var_double(ncid, lon_varid, worknc1)
          if (retval .ne. nf_noerr) call handle_err(retval)
        endif ! wcoord
      endif ! mnproc

!---  Write one timestep of data 
      if (mnproc.eq.1) then 
          count(1) = itdm
          count(2) = jtdm
          count(3) = kdm
          count(4) = 1
          start(1) = 1
          start(2) = 1
          start(3) = 1
          start(4) = nrec+1
      endif ! mnproc


!---  Write var_nc3d 
      do k=1,kdm
        do j=1,jj
          do i=1,ii
            util1(i,j) = var_nc3d(i,j,k)
          enddo
        enddo
        call xcaget(worknc1, util1, 1) !! get value from all procs
        worknc2(:,:,k)=worknc1
      enddo 
      if (mnproc.eq.1) then
          retval = nf_put_vara_double(ncid, var_varid, start, count, 
     &        worknc2(1:itdm,1:jtdm,1:kdm))
          if (retval .ne. nf_noerr) call handle_err(retval)
      endif ! mnproc
     

!---  Close the file
      if (mnproc.eq.1) then
          retval = nf_close(ncid)
          if (retval .ne. nf_noerr) call handle_err(retval)
      endif ! mnproc

 
      end subroutine ncput_3d

      subroutine handle_err(errcode)
      implicit none

!--- To handle errors in the creation and filling of the netcdf 

      include 'netcdf.inc'

      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2

      end subroutine handle_err

      end module mod_nc
      

      
