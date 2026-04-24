module netcdfflux
  use typy
  use netcdf
  use ncglobvars
  use nctools
  

  public :: ncflux_init
  public :: ncflux_get_xy
  public :: ncflux_close

  contains

    subroutine ncflux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use globals 
      use global_objs
      use pde_objs
      use geom_tools
      use debug_tools
      use datetime
    
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
      
      real(kind=rkind), dimension(2) :: xy
      
      call getcoor(quadpnt, xy)
      
      
      
  
    end subroutine ncflux
    
    


	subroutine ncflux_init(ok, errmsg)
	  logical, intent(out) :: ok
	  character(len=*), intent(out) :: errmsg

	  integer :: ierr
	  integer :: dimid
	  integer :: n_default
	  integer, dimension(:), allocatable :: time_tmp

	  ok = .false.
	  errmsg = "ncflux_init: unknown error"

	  ierr = nf90_inq_varid(netcdfID, trim(ncfluxdata%lon_name), ncfluxdata%lon_varid)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot find longitude variable '"//trim(ncfluxdata%lon_name)//"': "//trim(nf90_strerror(ierr))
		return
	  end if

	  ierr = nf90_inq_varid(netcdfID, trim(ncfluxdata%lat_name), ncfluxdata%lat_varid)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot find latitude variable '"//trim(ncfluxdata%lat_name)//"': "//trim(nf90_strerror(ierr))
		return
	  end if

	  ierr = nf90_inq_varid(netcdfID, trim(ncfluxdata%time_name), ncfluxdata%time_varid)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot find time variable '"//trim(ncfluxdata%time_name)//"': "//trim(nf90_strerror(ierr))
		return
	  end if

	  ierr = nf90_inq_varid(netcdfID, trim(ncfluxdata%var_name), ncfluxdata%q_varid)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot find flux variable '"//trim(ncfluxdata%var_name)//"': "//trim(nf90_strerror(ierr))
		return
	  end if

	  ierr = nf90_inq_dimid(netcdfID, trim(ncfluxdata%lon_name), dimid)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot find lon dimension '"//trim(ncfluxdata%lon_name)//"': "//trim(nf90_strerror(ierr))
		return
	  end if

	  ierr = nf90_inquire_dimension(netcdfID, dimid, len=n_default)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot inquire lon dimension length: "//trim(nf90_strerror(ierr))
		return
	  end if
	  ncfluxdata%nlon = int(n_default, kind=ikind)

	  ierr = nf90_inq_dimid(netcdfID, trim(ncfluxdata%lat_name), dimid)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot find lat dimension '"//trim(ncfluxdata%lat_name)//"': "//trim(nf90_strerror(ierr))
		return
	  end if

	  ierr = nf90_inquire_dimension(netcdfID, dimid, len=n_default)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot inquire lat dimension length: "//trim(nf90_strerror(ierr))
		return
	  end if
	  ncfluxdata%nlat = int(n_default, kind=ikind)

	  ierr = nf90_inq_dimid(netcdfID, trim(ncfluxdata%time_name), dimid)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot find time dimension '"//trim(ncfluxdata%time_name)//"': "//trim(nf90_strerror(ierr))
		return
	  end if

	  ierr = nf90_inquire_dimension(netcdfID, dimid, len=n_default)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: cannot inquire time dimension length: "//trim(nf90_strerror(ierr))
		return
	  end if
	  ncfluxdata%ntime = int(n_default, kind=ikind)

	  if (allocated(ncfluxdata%lon)) deallocate(ncfluxdata%lon)
	  if (allocated(ncfluxdata%lat)) deallocate(ncfluxdata%lat)
	  if (allocated(ncfluxdata%time)) deallocate(ncfluxdata%time)
	  if (allocated(ncfluxdata%qslice)) deallocate(ncfluxdata%qslice)

	  allocate(ncfluxdata%lon(ncfluxdata%nlon))
	  allocate(ncfluxdata%lat(ncfluxdata%nlat))
	  allocate(ncfluxdata%time(ncfluxdata%ntime))
	  allocate(ncfluxdata%qslice(ncfluxdata%nlat, ncfluxdata%nlon))

	  ierr = nf90_get_var(netcdfID, ncfluxdata%lon_varid, ncfluxdata%lon)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: failed to read lon variable: "//trim(nf90_strerror(ierr))
		return
	  end if

	  ierr = nf90_get_var(netcdfID, ncfluxdata%lat_varid, ncfluxdata%lat)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: failed to read lat variable: "//trim(nf90_strerror(ierr))
		return
	  end if

	  allocate(time_tmp(ncfluxdata%ntime))
	  ierr = nf90_get_var(netcdfID, ncfluxdata%time_varid, time_tmp)
	  if (ierr /= nf90_noerr) then
		errmsg = "ncflux_init: failed to read time variable: "//trim(nf90_strerror(ierr))
		deallocate(time_tmp)
		return
	  end if

	  ncfluxdata%time = int(time_tmp, kind=ikind)
	  deallocate(time_tmp)

	  ncfluxdata%has_fill = .false.
	  ierr = nf90_get_att(netcdfID, ncfluxdata%q_varid, "_FillValue", ncfluxdata%fill_value)
	  if (ierr == nf90_noerr) then
		ncfluxdata%has_fill = .true.
	  else
		ierr = nf90_get_att(netcdfID, ncfluxdata%q_varid, "missing_value", ncfluxdata%fill_value)
		if (ierr == nf90_noerr) ncfluxdata%has_fill = .true.
	  end if

	  ncfluxdata%initialized = .true.
	  ncfluxdata%slice_loaded = .false.
	  ncfluxdata%current_time_index = -1_ikind

	  ok = .true.
	  write(errmsg, '(a, i0, a, i0, a, i0, a, es14.6, a, es14.6, a, es14.6, a, es14.6, a, i0, a, i0)') &
		"ncflux_init: everything ok; nlon=", ncfluxdata%nlon, &
		", nlat=", ncfluxdata%nlat, &
		", ntime=", ncfluxdata%ntime, &
		", lon range=", minval(ncfluxdata%lon), " to ", maxval(ncfluxdata%lon), &
		", lat range=", minval(ncfluxdata%lat), " to ", maxval(ncfluxdata%lat), &
		", time range=", minval(ncfluxdata%time), " to ", maxval(ncfluxdata%time)

	end subroutine ncflux_init


	 subroutine ncflux_get_xy(x, y, cur_hrs, qval, ok, errmsg)
	  real(kind=rkind), intent(in) :: x, y
	  integer(kind=ikind), intent(in) :: cur_hrs
	  real(kind=rkind), intent(out) :: qval
	  logical, intent(out) :: ok
	  character(len=*), intent(out) :: errmsg

	  real(kind=rkind) :: lat0, lon0
	  integer(kind=ikind) :: tidx

	  ok = .false.
	  qval = ncfluxdata%fill_value
	  errmsg = "unknown error"

	  if (.not. ncfluxdata%initialized) then
		errmsg = "ncfluxdata not initialized"
		return
	  end if

	  call find_time_index(cur_hrs, tidx, ok)

	  if (.not. ok) then
		errmsg = "time out of range or not present in NetCDF"
		return
	  end if

	  call load_qslice(tidx, ok, errmsg)

	  if (.not. ok) then
		return
	  end if

	  call utm2latlong(x, y, lat0, lon0)

	  call interpolate_flux_latlon(lat0, lon0, qval, ok, errmsg)

	  if (ok) then
		errmsg = "everything ok"
	  end if

	end subroutine ncflux_get_xy

	subroutine load_qslice(tidx, ok, errmsg)
	  integer(kind=ikind), intent(in) :: tidx
	  logical, intent(out) :: ok
	  character(len=*), intent(out) :: errmsg

	  integer :: ierr
	  integer, dimension(3) :: start, count
	  real(kind=rkind), dimension(:,:,:), allocatable :: tmp

	  ok = .false.
	  errmsg = "unknown error"

	  if (tidx < 1_ikind .or. tidx > ncfluxdata%ntime) then
		errmsg = "time index outside NetCDF range"
		return
	  end if

	  if (ncfluxdata%slice_loaded) then
		if (ncfluxdata%current_time_index == tidx) then
		  ok = .true.
		  errmsg = "everything ok"
		  return
		end if
	  end if

	  allocate(tmp(ncfluxdata%nlon, ncfluxdata%nlat, 1))

	  start = (/ 1, 1, int(tidx) /)
	  count = (/ int(ncfluxdata%nlon), int(ncfluxdata%nlat), 1 /)

	  ierr = nf90_get_var(netcdfID, ncfluxdata%q_varid, tmp, start=start, count=count)

	  if (ierr /= nf90_noerr) then
		errmsg = "NetCDF read error: "//trim(nf90_strerror(ierr))
		deallocate(tmp)
		return
	  end if

	  ncfluxdata%qslice = transpose(tmp(:,:,1))
	  deallocate(tmp)

	  ncfluxdata%current_time_index = tidx
	  ncfluxdata%slice_loaded = .true.

	  ok = .true.
	  errmsg = "everything ok"

	end subroutine load_qslice


	 subroutine interpolate_flux_latlon(lat0, lon0, qval, ok, errmsg)
	  real(kind=rkind), intent(in) :: lat0, lon0
	  real(kind=rkind), intent(out) :: qval
	  logical, intent(out) :: ok
	  character(len=*), intent(out) :: errmsg

	  integer(kind=ikind) :: i1, i2, j1, j2
	  logical :: ok_local
	  real(kind=rkind) :: x, y
	  real(kind=rkind) :: wx, wy
	  real(kind=rkind) :: q11, q12, q21, q22

	  ok = .false.
	  errmsg = "unknown error"
	  qval = ncfluxdata%fill_value

	  x = adjust_longitude_to_grid(ncfluxdata%lon, lon0)
	  y = lat0

	  call binary_search_bracket(ncfluxdata%lat, y, i1, i2, ok_local)
	  if (.not. ok_local) then
		errmsg = "latitude out of NetCDF range"
		return
	  end if

	  call binary_search_bracket(ncfluxdata%lon, x, j1, j2, ok_local)
	  if (.not. ok_local) then
		errmsg = "longitude out of NetCDF range"
		return
	  end if

	  wx = (x - ncfluxdata%lon(j1))/(ncfluxdata%lon(j2) - ncfluxdata%lon(j1))
	  wy = (y - ncfluxdata%lat(i1))/(ncfluxdata%lat(i2) - ncfluxdata%lat(i1))

	  q11 = ncfluxdata%qslice(i1,j1)
	  q12 = ncfluxdata%qslice(i1,j2)
	  q21 = ncfluxdata%qslice(i2,j1)
	  q22 = ncfluxdata%qslice(i2,j2)

	  if (ncfluxdata%has_fill) then
		if (is_fill(q11, ncfluxdata%fill_value) .or. &
			is_fill(q12, ncfluxdata%fill_value) .or. &
			is_fill(q21, ncfluxdata%fill_value) .or. &
			is_fill(q22, ncfluxdata%fill_value)) then

		  errmsg = "NetCDF returned fill value (no valid data)"
		  return
		end if
	  end if

	  qval = bilinear(q11, q12, q21, q22, wx, wy)

	  ok = .true.
	  errmsg = "everything ok"

	end subroutine interpolate_flux_latlon


	  subroutine find_time_index(cur_hrs, tidx, ok)
		integer(kind=ikind), intent(in) :: cur_hrs
		integer(kind=ikind), intent(out) :: tidx
		logical, intent(out) :: ok

		integer(kind=ikind) :: lo, hi, mid

		ok = .false.
		tidx = -1_ikind

		if (ncfluxdata%ntime < 1_ikind) return

		if (cur_hrs < minval(ncfluxdata%time) .or. cur_hrs > maxval(ncfluxdata%time)) then
		  return
		end if

		lo = 1_ikind
		hi = ncfluxdata%ntime

		do while (lo <= hi)
		  mid = (lo + hi)/2_ikind

		  if (ncfluxdata%time(mid) == cur_hrs) then
			tidx = mid
			ok = .true.
			return
		  else if (ncfluxdata%time(mid) < cur_hrs) then
			lo = mid + 1_ikind
		  else
			hi = mid - 1_ikind
		  end if
		end do

		print *, "find_time_index: cur_hrs is within NetCDF time range but exact time was not found"
		print *, "  requested cur_hrs = ", cur_hrs
		print *, "  time range        = ", minval(ncfluxdata%time), maxval(ncfluxdata%time)
	  end subroutine find_time_index


	  subroutine binary_search_bracket(arr, x, i1, i2, ok)
		real(kind=rkind), dimension(:), intent(in) :: arr
		real(kind=rkind), intent(in) :: x
		integer(kind=ikind), intent(out) :: i1, i2
		logical, intent(out) :: ok

		integer(kind=ikind) :: lo, hi, mid, n
		logical :: ascending

		ok = .false.
		i1 = -1_ikind
		i2 = -1_ikind

		n = size(arr, kind=ikind)
		if (n < 2_ikind) return

		ascending = arr(n) > arr(1)

		if (ascending) then
		  if (x < arr(1) .or. x > arr(n)) return
		else
		  if (x > arr(1) .or. x < arr(n)) return
		end if

		if (abs(x - arr(n)) < 100.0_rkind*epsilon(1.0_rkind)) then
		  i1 = n - 1_ikind
		  i2 = n
		  ok = .true.
		  return
		end if

		lo = 1_ikind
		hi = n

		do while (hi - lo > 1_ikind)
		  mid = (lo + hi)/2_ikind

		  if (ascending) then
			if (x >= arr(mid)) then
			  lo = mid
			else
			  hi = mid
			end if
		  else
			if (x <= arr(mid)) then
			  lo = mid
			else
			  hi = mid
			end if
		  end if
		end do

		i1 = lo
		i2 = hi
		ok = .true.
	  end subroutine binary_search_bracket


	  pure function adjust_longitude_to_grid(lon_arr, lon0) result(x)
		real(kind=rkind), dimension(:), intent(in) :: lon_arr
		real(kind=rkind), intent(in) :: lon0
		real(kind=rkind) :: x
		real(kind=rkind) :: lonmin, lonmax

		lonmin = minval(lon_arr)
		lonmax = maxval(lon_arr)
		x = lon0

		if (lonmin >= 0.0_rkind .and. x < 0.0_rkind) x = x + 360.0_rkind
		if (lonmax <= 180.0_rkind .and. x > 180.0_rkind) x = x - 360.0_rkind
	  end function adjust_longitude_to_grid


	  pure function bilinear(q11, q12, q21, q22, wx, wy) result(q)
		real(kind=rkind), intent(in) :: q11, q12, q21, q22
		real(kind=rkind), intent(in) :: wx, wy
		real(kind=rkind) :: q

		q = (1.0_rkind - wx)*(1.0_rkind - wy)*q11 + &
			 wx            *(1.0_rkind - wy)*q12 + &
			(1.0_rkind - wx)*wy            *q21 + &
			 wx            *wy            *q22
	  end function bilinear


	  pure logical function is_fill(x, fill_value)
		real(kind=rkind), intent(in) :: x, fill_value

		is_fill = abs(x - fill_value) < &
				  100.0_rkind*epsilon(1.0_rkind)*max(1.0_rkind, abs(fill_value))
	  end function is_fill


	  subroutine ncflux_close()
		if (allocated(ncfluxdata%lon)) deallocate(ncfluxdata%lon)
		if (allocated(ncfluxdata%lat)) deallocate(ncfluxdata%lat)
		if (allocated(ncfluxdata%time)) deallocate(ncfluxdata%time)
		if (allocated(ncfluxdata%qslice)) deallocate(ncfluxdata%qslice)

		ncfluxdata%initialized = .false.
		ncfluxdata%slice_loaded = .false.
		ncfluxdata%current_time_index = -1_ikind
	  end subroutine ncflux_close


end module netcdfflux
