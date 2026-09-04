module init_netcdf

  contains
  
    subroutine netcdf()
      use typy
      use netcdf
      use ncglobvars
      use datetime
      use globals
      use global_objs
      use nctools
      use debug_tools
      use ncdem
      use netcdfflux
      use core_tools
      use pde_objs
      use readtools
      use ncfluxarea
      use ncmesh
      use ncmap
      use geom_tools
            
      integer :: ierr, filetmp, fileconf
      integer(kind=ikind) :: i, bccnt, j
      logical :: success
      real(kind=rkind) :: q, segment_distance, best_segment_distance
      character(len=1024) :: errmsg
      integer(kind=ikind), dimension(3) :: datearray
      real(kind=rkind), dimension(2) :: xy, A, B, C
      real(kind=rkind), dimension(4,2) :: pts
      
      starttime%year = 2010
      starttime%month = 1
      starttime%day = 1
      starttime%hour = 0
      starttime%minute = 0
      starttime%second = 0
      
      
      
      
      open(newunit=fileconf, file="drutes.conf/netcdf/netcdf.conf", action="read", status="old")
      
!      read_int_array(r, fileid, ranges, errmsg, checklen, noexit)

      call fileread(datearray, fileconf)
      
      if (datearray(1) < 1900) then
        write(errmsg, *) "----------------------------------", new_line("a"), &
                        "W: the year of your simulation start is very old, check drutes.conf/netcdf/netcdf.conf", new_line("a"), &
                        "the year you defined is:", datearray(1), new_line("a"), "check your inputs! however simulation continues",&
                        new_line("a"), "----------------------------------"
        call write_log(cut(errmsg))
      end if
      
      starttime%year = datearray(1)
      
      if (datearray(2) < 1 .or. datearray(2) > 12) then
        write(errmsg, *) "incorrect month specified, month should be within 1-12, you specified", datearray(2)
        call file_error(fileconf, errmsg)
      end if
     
      starttime%month = datearray(2)
      
      starttime%day = datearray(3)
    
      
      call fileread(LSdisp, fileconf, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
                    errmsg="incorrect dispersivity definition in drutes.conf/netcdf/netcdf.conf")
      
      
      call fileread(Qmin, fileconf, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
                    errmsg="incorrect Qmin definition in drutes.conf/netcdf/netcdf.conf")
    
      call fileread(cinit_ls, fileconf, ranges=(/0.0_rkind, huge(0.0_rkind)/), &
                    errmsg="incorrect c initial definition in drutes.conf/netcdf/netcdf.conf")

      
      
      addedbc = maxval(nodes%edge) + 1
      
      
      ierr = nf90_open(path="drutes.conf/netcdf/mRM_Fluxes_States.nc", mode=nf90_nowrite, ncid=netcdfID)
      
      if (ierr /= nf90_noerr) ERROR STOP "Error opening NetCDF file: drutes.conf/netcdf/mRM_Fluxes_States.nc"
    
      call read_ncorigin(netcdfID, ncstart)
      

      
      call ncflux_init(success, errmsg)
   
      if (.not. success) then
        print *, cut(errmsg)
        print *, "unsupported structure of the file drutes.conf/netcdf/mRM_Fluxes_States.nc"
        print *, "is this correct output from mHM?"
        print *, "after succesfull mHM simulation you should copy "
        print *, "      mRM_Fluxes_States.nc -> drutes.conf/netcdf/mRM_Fluxes_States.nc"
        ERROR STOP
      end if
      
      call read_ncbounds(netcdfID, success, errmsg)
      
      if (.not. success) then
        print *, cut(errmsg)
        print *, "unsupported structure of the file drutes.conf/netcdf/mRM_Fluxes_States.nc"
        print *, "is this correct output from mHM?"
        print *, "after succesfull mHM simulation you should copy "
        print *, "      mRM_Fluxes_States.nc -> drutes.conf/netcdf/mRM_Fluxes_States.nc"
        ERROR STOP
      end if
    
      call getmeshalt()
      
      ora_di_ini = difftime(ncstart, starttime, "hrs")
      
    
      do i=1, nodes%kolik
      
        call ncflux_get_xy_cell(nodes%data(i,1), nodes%data(i,2), ora_di_ini, q, success, errmsg)
        if (.not. success) then
          nodes%edge(i) = addedbc
        else if (q < Qmin) then
          nodes%edge(i) = addedbc
        end if
        
!         if (q > 0.0) then
!           if (nint(nodealt(i)) == missing) then
!             write(errmsg, *) "W: your dem model doesn't conver the entire watershed, update drutes.conf/netcdf/dem.nc, node:", i, &
!               "will be deactivated"
!             call write_log(errmsg)
!             nodes%edge(i) = addedbc
!           end if
!         end if
        
        
      end do

      allocate(ncfluxdata%activeel(elements%kolik))
      
      do i=1, elements%kolik
        bccnt = 0
        do j=1, ubound(elements%data,2)
          if (nodes%edge(elements%data(i,j)) == addedbc) then
            bccnt = bccnt + 1
          end if
        end do
        if (bccnt == ubound(elements%data,2)) then
          ncfluxdata%activeel(i) = .false.
        else
          ncfluxdata%activeel(i) = .true.
        end if
      end do
      
!       ncfluxdata%activeel(:) = .true.
      
      allocate(ncfluxdata%cellarea(elements%kolik))
      allocate(ncfluxdata%fluxvct(elements%kolik,2))
      
      ncfluxdata%fluxvct = 0.0_rkind
      
      
      
      do i = 1, elements%kolik
        xy(1) = avgarr(nodes%data(elements%data(i,:),1))
        xy(2) = avgarr(nodes%data(elements%data(i,:),2))
        call ncflux_cell_area_xy(xy(1), xy(2), ncfluxdata%cellarea(i), success, errmsg)
      end do
      
      
      call getncmesh(netcdfID, ncnodes, ncelements, success, errmsg)

      if (.not. success) then
        print *, trim(errmsg)
        error stop
      end if
      
      call ncelslope()
      
      call mapel()
      
      call terrain_slopes()
      
      call readchannel()
      
      do i=1, elements%kolik
        if (ncfluxdata%activeel(i)) then 
          C(1) = mean_array(nodes%data(elements%data(i,:),1))
          C(2) = mean_array(nodes%data(elements%data(i,:),2))
          best_segment_distance = huge(1.0_rkind)
          channel: do j=1, channel_el%kolik
                    A = channel_nd%data(channel_el%data(j,1),:)
                    B = channel_nd%data(channel_el%data(j,2),:)
                    segment_distance = point_segment_distance(A, B, C)
                    if (segment_distance < best_segment_distance) then
                      best_segment_distance = segment_distance
                      ncfluxdata%fluxvct(i,:) = unit_vector(A,B)
                    end if
          end do channel
        end if
      end do
      
      
      allocate(ncelements%areas(ncelements%kolik))
      
      do i=1, ncelements%kolik
        do j=1,4
          pts(j,:) = ncnodes%data(ncelements%data(i,j),1:2)
        end do

        ncelements%areas(i) = quad_area(pts(1,:), pts(2,:), pts(3,:), pts(4,:))
      end do
      
      
      
      call readbcvals(unitW=fileconf, struct=pde(1)%bc, dimen=2_ikind, &
        dirname="drutes.conf/netcdf/")


      pde(1)%problem_name(1) = "ADE_in_watershed"
      pde(1)%problem_name(2) = "Advection-dispersion-reaction equation for watershed large scale"

      pde(1)%solution_name(1) = "solute_concentration" !nazev vystupnich souboru
      pde(1)%solution_name(2) = "c  [M/L^3]" !popisek grafu

      pde(1)%flux_name(1) = "conc_flux"  
      pde(1)%flux_name(2) = "concentration flux [M.L^{-2}.T^{-1}]"
      
      allocate(pde(1)%mass_name(1,2))

      pde(1)%mass_name(1,1) = "conc_in_river"
      pde(1)%mass_name(1,2) = "concetration [M/L^3]"
      
      pde(1)%print_mass = .true.      

    end subroutine netcdf
    
    
    subroutine read_ncorigin(ncid, ncstart)
      use netcdf
      use datetime
      implicit none

      integer, intent(in) :: ncid
      type(datetime_t), intent(out) :: ncstart

      integer :: ierr, time_varid, p
      character(len=256) :: units
      character(len=256) :: refstr

      ! find variable "time"
      ierr = nf90_inq_varid(ncid, "time", time_varid)
      if (ierr /= nf90_noerr) then
        print *, "Cannot find variable 'time': ", trim(nf90_strerror(ierr))
        stop
      end if

      ! read attribute units, e.g. "hours since 1950-01-01 00:00:00"
      ierr = nf90_get_att(ncid, time_varid, "units", units)
      if (ierr /= nf90_noerr) then
        print *, "Cannot read attribute 'units': ", trim(nf90_strerror(ierr))
        stop
      end if

      units = trim(adjustl(units))

      ! locate the word "since"
      p = index(units, "since")
      if (p <= 0) then
        print *, "Unsupported time units format: ", trim(units)
        stop
      end if

      ! take everything after "since"
      refstr = adjustl(units(p+5:))
      refstr = trim(refstr)

      ! initialize defaults
      ncstart%year   = 0
      ncstart%month  = 0
      ncstart%day    = 0
      ncstart%hour   = 0
      ncstart%minute = 0
      ncstart%second = 0

      ! support formats:
      !   YYYY-MM-DD HH:MM:SS
      !   YYYY-MM-DDTHH:MM:SSZ
      !   YYYY-MM-DD
      if (len_trim(refstr) >= 19) then
        read(refstr(1:4),  *) ncstart%year
        read(refstr(6:7),  *) ncstart%month
        read(refstr(9:10), *) ncstart%day

        if (refstr(11:11) == 'T' .or. refstr(11:11) == ' ') then
          read(refstr(12:13), *) ncstart%hour
          read(refstr(15:16), *) ncstart%minute
          read(refstr(18:19), *) ncstart%second
        end if

      else if (len_trim(refstr) >= 10) then
        read(refstr(1:4),  *) ncstart%year
        read(refstr(6:7),  *) ncstart%month
        read(refstr(9:10), *) ncstart%day
      else
        print *, "Reference date in units is too short: ", trim(refstr)
        stop
      end if

    end subroutine read_ncorigin
    
    
    subroutine read_ncbounds(ncid, ok, errmsg)
      use netcdf
      use typy
      use ncglobvars
      

      integer, intent(in) :: ncid
      logical, intent(out) :: ok
      character(len=*), intent(out) :: errmsg

      integer :: ierr
      integer :: lat_bnds_varid, lon_bnds_varid
      real(kind=rkind), allocatable :: tmp_lat_bnds(:,:)
      real(kind=rkind), allocatable :: tmp_lon_bnds(:,:)

      ok = .false.
      errmsg = "read_ncbounds: unknown error"

      ncfluxdata%has_bounds = .false.

      ierr = nf90_inq_varid(ncid, "lat_bnds", lat_bnds_varid)
      if (ierr /= nf90_noerr) then
        errmsg = "read_ncbounds: cannot find lat_bnds: " // trim(nf90_strerror(ierr))
        return
      end if

      ierr = nf90_inq_varid(ncid, "lon_bnds", lon_bnds_varid)
      if (ierr /= nf90_noerr) then
        errmsg = "read_ncbounds: cannot find lon_bnds: " // trim(nf90_strerror(ierr))
        return
      end if

      if (allocated(ncfluxdata%lat_bnds)) deallocate(ncfluxdata%lat_bnds)
      if (allocated(ncfluxdata%lon_bnds)) deallocate(ncfluxdata%lon_bnds)

      ! ncdump shows:
      !   lat_bnds(lat,bnds)
      !   lon_bnds(lon,bnds)
      !
      ! For NetCDF Fortran, read into reversed shape:
      !   tmp_lat_bnds(bnds,lat)
      !   tmp_lon_bnds(bnds,lon)
      allocate(tmp_lat_bnds(2, ncfluxdata%nlat))
      allocate(tmp_lon_bnds(2, ncfluxdata%nlon))

      ierr = nf90_get_var(ncid, lat_bnds_varid, tmp_lat_bnds)
      if (ierr /= nf90_noerr) then
        errmsg = "read_ncbounds: cannot read lat_bnds: " // trim(nf90_strerror(ierr))
        deallocate(tmp_lat_bnds, tmp_lon_bnds)
        return
      end if

      ierr = nf90_get_var(ncid, lon_bnds_varid, tmp_lon_bnds)
      if (ierr /= nf90_noerr) then
        errmsg = "read_ncbounds: cannot read lon_bnds: " // trim(nf90_strerror(ierr))
        deallocate(tmp_lat_bnds, tmp_lon_bnds)
        return
      end if

      ! Store internally as bounds(index,1:2)
      allocate(ncfluxdata%lat_bnds(ncfluxdata%nlat, 2))
      allocate(ncfluxdata%lon_bnds(ncfluxdata%nlon, 2))

      ncfluxdata%lat_bnds = transpose(tmp_lat_bnds)
      ncfluxdata%lon_bnds = transpose(tmp_lon_bnds)

      deallocate(tmp_lat_bnds, tmp_lon_bnds)

      ncfluxdata%has_bounds = .true.
      ok = .true.
      errmsg = "everything ok"
      
    end subroutine read_ncbounds
    
    
    subroutine readchannel()
      use typy
      use ncglobvars
      use readtools
      use debug_tools
      
      
      integer :: fileid, ierr
      integer(kind=ikind) :: counter, i
      real(kind=rkind), dimension(2) :: tmp
      
      open(newunit=fileid, file="drutes.conf/netcdf/channel.dat", iostat=ierr, status="old", action="read")
      
      if (ierr /= 0) then
        print *, "unable to open file drutes.conf/netcdf/channel.dat"
        ERROR STOP
      end if
      
      counter = 0
      
      do 
        call comment(fileid)
        read(fileid, fmt=*, iostat=ierr) tmp
        
        if (ierr == 0) then
          counter = counter + 1
        else
          EXIT
        end if
      end do
      
      if (counter < 2) then
        print *, "file drutes.conf/netcdf/channel.dat doesn't contain enough values, check the file! "
        ERROR STOP
      end if
      
      allocate(channel_nd%data(counter,2))
      channel_nd%kolik = counter
      allocate(channel_el%data(counter-1, 2))
      channel_el%kolik = counter - 1
      
      close(fileid)
      
      open(newunit=fileid, file="drutes.conf/netcdf/channel.dat", iostat=ierr, status="old", action="read")
      
      do i=1, counter
        call comment(fileid)
        read(fileid, fmt=*, iostat=ierr) channel_nd%data(i,:)
      end do
      
      do i=1, counter - 1
        channel_el%data(i,1) = i
        channel_el%data(i,2) = i+1
      end do

  
      
        
    
    end subroutine readchannel 






end module init_netcdf
