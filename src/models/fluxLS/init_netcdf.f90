module init_netcdf

  contains
  
    subroutine netcdf()
      use typy
      use netcdf
      use ncglobvars
      
      integer :: ierr
      
      ierr = nf90_open(path="drutes.conf/netcdf/flux.nc", mode=nf90_nowrite, netcdfID=netcdfID)
      
      if (ierr /= nf90_noerr) ERROR STOP "Error opening NetCDF file"
       
      ! Get variable ID for Qrouted
      ierr = nf90_inq_varid(netcdfID, "Qrouted", varid)
      
        ! Get dimension IDs
      ierr = nf90_inq_dimid(netcdfID, "time", dimid_time)
      ierr = nf90_inq_dimid(netcdfID, "lat", dimid_lat)
      ierr = nf90_inq_dimid(netcdfID, "lon", dimid_lon)
      
        ! Get dimension lengths
      ierr = nf90_inquire_dimension(netcdfID, dimid_time, len = time_len)
      ierr = nf90_inquire_dimension(netcdfID, dimid_lat, len = lat_len)
      ierr = nf90_inquire_dimension(netcdfID, dimid_lon, len = lon_len)
       
    
    end subroutine netcdf






end module init_netcdf
