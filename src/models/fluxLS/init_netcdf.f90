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
      
      integer :: ierr
      integer(kind=ikind) :: i
      logical :: success
      real(kind=rkind) :: q
      character(len=1024) :: errmsg
      
      starttime%year = 2010
      starttime%month = 1
      starttime%day = 1
      starttime%hour = 0
      starttime%minute = 0
      starttime%second = 0
      
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
    
      call getmeshalt()
      
      ore_di_ini = difftime(ncstart, starttime, "hrs")
      
!      print *, ore_di_ini ; stop      
      do i=1, nodes%kolik
      
		call ncflux_get_xy(nodes%data(i,1), nodes%data(i,2), ore_di_ini, q, success, errmsg)
		
		print *, i, q
		

	  end do
      
      
      
      stop
       
      ! Get variable ID for Qrouted
!      ierr = nf90_inq_varid(netcdfID, "Qrouted", varid)
      
!        ! Get dimension IDs
!      ierr = nf90_inq_dimid(netcdfID, "time", dimid_time)
!      ierr = nf90_inq_dimid(netcdfID, "lat", dimid_lat)
!      ierr = nf90_inq_dimid(netcdfID, "lon", dimid_lon)
      
!        ! Get dimension lengths
!      ierr = nf90_inquire_dimension(netcdfID, dimid_time, len = time_len)
!      ierr = nf90_inquire_dimension(netcdfID, dimid_lat, len = lat_len)
!      ierr = nf90_inquire_dimension(netcdfID, dimid_lon, len = lon_len)
       
    
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
    
    







end module init_netcdf
