

module ncdem
  use netcdf
  use typy


  private :: dem_cache
  private :: dem_open, dem_close, dem_get_altitude

  type, private :: dem_cache
    logical :: is_open = .false.

    ! NetCDF-facing identifiers
    integer :: ncid = -1
    integer :: lat_varid = -1
    integer :: lon_varid = -1
    integer :: dem_varid = -1

    integer(kind=ikind) :: nlat = 0_ikind
    integer(kind=ikind) :: nlon = 0_ikind
    integer(kind=ikind) :: dem_order = 0_ikind
    ! 1 => original DEM variable is (lat,lon)
    ! 2 => original DEM variable is (lon,lat)

    logical :: has_fill = .false.
    real(kind=rkind) :: fill_value = -9999.0_rkind

    real(kind=rkind), dimension(:), allocatable :: lat, lon
    real(kind=rkind), dimension(:,:), allocatable :: z
    ! stored internally always as z(lat,lon)
  end type dem_cache

contains

  subroutine dem_open(dem, filename, ok)
    type(dem_cache), intent(inout) :: dem
    character(len=*), intent(in)   :: filename
    logical, intent(out)           :: ok

    integer :: ncerr
    integer :: nvars, ngatts, unlimdimid
    integer :: ndims
    integer :: lat_ndims, lon_ndims
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
    integer, dimension(NF90_MAX_VAR_DIMS) :: lat_dimids, lon_dimids
    real(kind=rkind), dimension(:,:), allocatable :: tmp

    ok = .false.

    if (dem%is_open) then
      call dem_close(dem)
    end if

    ncerr = nf90_open(trim(filename), NF90_NOWRITE, dem%ncid)
    if (ncerr /= nf90_noerr) return

    ncerr = nf90_inquire(dem%ncid, nVariables=nvars, nAttributes=ngatts, unlimitedDimId=unlimdimid)
    if (ncerr /= nf90_noerr) then
      call close_ncid_only(dem)
      return
    end if

    call find_coordinate_variables(dem%ncid, nvars, dem%lat_varid, dem%lon_varid, ok)
    if (.not. ok) then
      call close_ncid_only(dem)
      return
    end if

    call find_dem_variable(dem%ncid, nvars, dem%lat_varid, dem%lon_varid, dem%dem_varid, ok)
    if (.not. ok) then
      call close_ncid_only(dem)
      return
    end if

    call read_axis_1d(dem%ncid, dem%lat_varid, dem%lat, dem%nlat, ok)
    if (.not. ok) then
      call close_ncid_only(dem)
      return
    end if

    call read_axis_1d(dem%ncid, dem%lon_varid, dem%lon, dem%nlon, ok)
    if (.not. ok) then
      call close_ncid_only(dem)
      return
    end if

    ncerr = nf90_inquire_variable(dem%ncid, dem%dem_varid, ndims=ndims, dimids=dimids)
    if (ncerr /= nf90_noerr) then
      call close_ncid_only(dem)
      return
    end if

    if (ndims /= 2) then
      call close_ncid_only(dem)
      return
    end if

    ncerr = nf90_inquire_variable(dem%ncid, dem%lat_varid, ndims=lat_ndims, dimids=lat_dimids)
    if (ncerr /= nf90_noerr) then
      call close_ncid_only(dem)
      return
    end if

    ncerr = nf90_inquire_variable(dem%ncid, dem%lon_varid, ndims=lon_ndims, dimids=lon_dimids)
    if (ncerr /= nf90_noerr) then
      call close_ncid_only(dem)
      return
    end if

    if (dimids(1) == lat_dimids(1) .and. dimids(2) == lon_dimids(1)) then
      dem%dem_order = 1_ikind
      allocate(dem%z(dem%nlat, dem%nlon))
      ncerr = nf90_get_var(dem%ncid, dem%dem_varid, dem%z)
      if (ncerr /= nf90_noerr) then
        call close_ncid_only(dem)
        if (allocated(dem%z)) deallocate(dem%z)
        return
      end if

    else if (dimids(1) == lon_dimids(1) .and. dimids(2) == lat_dimids(1)) then
      dem%dem_order = 2_ikind
      allocate(tmp(dem%nlon, dem%nlat))
      ncerr = nf90_get_var(dem%ncid, dem%dem_varid, tmp)
      if (ncerr /= nf90_noerr) then
        if (allocated(tmp)) deallocate(tmp)
        call close_ncid_only(dem)
        return
      end if

      allocate(dem%z(dem%nlat, dem%nlon))
      dem%z = transpose(tmp)
      deallocate(tmp)

    else
      call close_ncid_only(dem)
      return
    end if

    dem%has_fill = .false.
    ncerr = nf90_get_att(dem%ncid, dem%dem_varid, "_FillValue", dem%fill_value)
    if (ncerr == nf90_noerr) then
      dem%has_fill = .true.
    else
      ncerr = nf90_get_att(dem%ncid, dem%dem_varid, "missing_value", dem%fill_value)
      if (ncerr == nf90_noerr) dem%has_fill = .true.
    end if

    call close_ncid_only(dem)

    dem%is_open = .true.
    ok = .true.
  end subroutine dem_open


   subroutine dem_get_altitude(dem, lat0, lon0, altitude, ok)
    type(dem_cache), intent(in)      :: dem
    real(kind=rkind), intent(in)     :: lat0, lon0
    real(kind=rkind), intent(out)    :: altitude
    logical, intent(out)             :: ok

    integer(kind=ikind) :: i1, i2, j1, j2
    logical :: ok_local
    real(kind=rkind) :: x, y
    real(kind=rkind) :: x1, x2, y1, y2
    real(kind=rkind) :: wx, wy
    real(kind=rkind) :: z11, z12, z21, z22

    ok = .false.
    altitude = dem%fill_value

    if (.not. dem%is_open) then
      print *, "dem_get_altitude: DEM is not open."
      return
    end if

    x = adjust_longitude_to_grid(dem%lon, lon0)
    y = lat0

    call binary_search_bracket(dem%lat, y, i1, i2, ok_local)
    if (.not. ok_local) then
      print *, "dem_get_altitude: latitude is outside DEM range."
      print *, "  requested latitude = ", y
      print *, "  DEM latitude range = ", minval(dem%lat), maxval(dem%lat)
      print *, "  original input lat/lon = ", lat0, lon0
      return
    end if

    call binary_search_bracket(dem%lon, x, j1, j2, ok_local)
    if (.not. ok_local) then
      print *, "dem_get_altitude: longitude is outside DEM range."
      print *, "  requested longitude = ", x
      print *, "  DEM longitude range = ", minval(dem%lon), maxval(dem%lon)
      print *, "  original input lat/lon = ", lat0, lon0
      return
    end if

    y1 = dem%lat(i1)
    y2 = dem%lat(i2)
    x1 = dem%lon(j1)
    x2 = dem%lon(j2)

    if (abs(x2 - x1) < epsilon(1.0_rkind)) then
      wx = 0.0_rkind
    else
      wx = (x - x1)/(x2 - x1)
    end if

    if (abs(y2 - y1) < epsilon(1.0_rkind)) then
      wy = 0.0_rkind
    else
      wy = (y - y1)/(y2 - y1)
    end if

    wx = max(0.0_rkind, min(1.0_rkind, wx))
    wy = max(0.0_rkind, min(1.0_rkind, wy))

    z11 = dem%z(i1, j1)
    z12 = dem%z(i1, j2)
    z21 = dem%z(i2, j1)
    z22 = dem%z(i2, j2)

    if (dem%has_fill) then
      if (is_fill(z11, dem%fill_value) .or. is_fill(z12, dem%fill_value) .or. &
          is_fill(z21, dem%fill_value) .or. is_fill(z22, dem%fill_value)) then
!        print *, "dem_get_altitude: DEM fill value encountered."
!        print *, "  input lat/lon      = ", lat0, lon0
!        print *, "  adjusted longitude = ", x
!        print *, "  surrounding values = ", z11, z12, z21, z22
        altitude = dem%fill_value
        return
      end if
    end if

    altitude = bilinear(z11, z12, z21, z22, wx, wy)
    ok = .true.
  end subroutine dem_get_altitude


  subroutine dem_close(dem)
    type(dem_cache), intent(inout) :: dem

    if (allocated(dem%lat)) deallocate(dem%lat)
    if (allocated(dem%lon)) deallocate(dem%lon)
    if (allocated(dem%z))   deallocate(dem%z)

    call close_ncid_only(dem)

    dem%is_open    = .false.
    dem%lat_varid  = -1
    dem%lon_varid  = -1
    dem%dem_varid  = -1
    dem%nlat       = 0_ikind
    dem%nlon       = 0_ikind
    dem%dem_order  = 0_ikind
    dem%has_fill   = .false.
    dem%fill_value = -9999.0_rkind
  end subroutine dem_close


  subroutine close_ncid_only(dem)
    type(dem_cache), intent(inout) :: dem
    integer :: ncerr

    if (dem%ncid >= 0) then
      ncerr = nf90_close(dem%ncid)
      dem%ncid = -1
    end if
  end subroutine close_ncid_only


  pure function bilinear(z11, z12, z21, z22, wx, wy) result(z)
    real(kind=rkind), intent(in) :: z11, z12, z21, z22, wx, wy
    real(kind=rkind) :: z

    z = (1.0_rkind - wx)*(1.0_rkind - wy)*z11 + &
         wx            *(1.0_rkind - wy)*z12 + &
         (1.0_rkind - wx)*wy           *z21 + &
         wx            *wy            *z22
  end function bilinear


  subroutine binary_search_bracket(arr, x, i1, i2, ok)
    real(kind=rkind), dimension(:), intent(in) :: arr
    real(kind=rkind), intent(in)               :: x
    integer(kind=ikind), intent(out)           :: i1, i2
    logical, intent(out)                       :: ok

    integer(kind=ikind) :: lo, hi, mid, n
    logical :: ascending

    ok = .false.
    i1 = -1_ikind
    i2 = -1_ikind

    n = size(arr, kind=ikind)
    if (n < 2_ikind) return

    ascending = (arr(n) > arr(1))

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
    real(kind=rkind), intent(in)               :: lon0
    real(kind=rkind) :: x
    real(kind=rkind) :: lonmin, lonmax

    lonmin = minval(lon_arr)
    lonmax = maxval(lon_arr)
    x = lon0

    if (lonmin >= 0.0_rkind .and. x < 0.0_rkind) x = x + 360.0_rkind
    if (lonmax <= 180.0_rkind .and. x > 180.0_rkind) x = x - 360.0_rkind
  end function adjust_longitude_to_grid


  pure logical function is_fill(z, fill_value)
    real(kind=rkind), intent(in) :: z, fill_value

    is_fill = abs(z - fill_value) < 100.0_rkind*epsilon(1.0_rkind)*max(1.0_rkind, abs(fill_value))
  end function is_fill


  subroutine find_coordinate_variables(ncid, nvars, lat_varid, lon_varid, ok)
    integer, intent(in)  :: ncid, nvars
    integer, intent(out) :: lat_varid, lon_varid
    logical, intent(out) :: ok

    integer :: varid, ncerr, ndims, natts, xtype
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimids
    character(len=NF90_MAX_NAME) :: name
    character(len=256) :: standard_name, units, long_name
    logical :: found_lat, found_lon

    ok = .false.
    found_lat = .false.
    found_lon = .false.
    lat_varid = -1
    lon_varid = -1

    do varid = 1, nvars
      ncerr = nf90_inquire_variable(ncid, varid, name=name, xtype=xtype, ndims=ndims, dimids=dimids, natts=natts)
      if (ncerr /= nf90_noerr) return

      if (ndims /= 1) cycle

      standard_name = ""
      units = ""
      long_name = ""

      call try_get_text_att(ncid, varid, "standard_name", standard_name)
      call try_get_text_att(ncid, varid, "units", units)
      call try_get_text_att(ncid, varid, "long_name", long_name)

      if (.not. found_lat) then
        if (is_latitude_variable(name, standard_name, units, long_name)) then
          lat_varid = varid
          found_lat = .true.
        end if
      end if

      if (.not. found_lon) then
        if (is_longitude_variable(name, standard_name, units, long_name)) then
          lon_varid = varid
          found_lon = .true.
        end if
      end if
    end do

    ok = found_lat .and. found_lon
  end subroutine find_coordinate_variables


  subroutine find_dem_variable(ncid, nvars, lat_varid, lon_varid, dem_varid, ok)
    integer, intent(in)  :: ncid, nvars, lat_varid, lon_varid
    integer, intent(out) :: dem_varid
    logical, intent(out) :: ok

    integer :: ncerr
    integer :: lat_ndims, lon_ndims
    integer :: lat_dimid, lon_dimid
    integer :: varid, ndims, natts, xtype
    integer, dimension(NF90_MAX_VAR_DIMS) :: lat_dimids, lon_dimids, dimids
    character(len=NF90_MAX_NAME) :: name
    character(len=256) :: standard_name, long_name, units

    ok = .false.
    dem_varid = -1

    ncerr = nf90_inquire_variable(ncid, lat_varid, ndims=lat_ndims, dimids=lat_dimids)
    if (ncerr /= nf90_noerr) return

    ncerr = nf90_inquire_variable(ncid, lon_varid, ndims=lon_ndims, dimids=lon_dimids)
    if (ncerr /= nf90_noerr) return

    lat_dimid = lat_dimids(1)
    lon_dimid = lon_dimids(1)

    do varid = 1, nvars
      if (varid == lat_varid .or. varid == lon_varid) cycle

      ncerr = nf90_inquire_variable(ncid, varid, name=name, xtype=xtype, ndims=ndims, dimids=dimids, natts=natts)
      if (ncerr /= nf90_noerr) return

      if (ndims /= 2) cycle

      if (.not. ((dimids(1) == lat_dimid .and. dimids(2) == lon_dimid) .or. &
                 (dimids(1) == lon_dimid .and. dimids(2) == lat_dimid))) cycle

      standard_name = ""
      long_name = ""
      units = ""

      call try_get_text_att(ncid, varid, "standard_name", standard_name)
      call try_get_text_att(ncid, varid, "long_name", long_name)
      call try_get_text_att(ncid, varid, "units", units)

      if (looks_like_dem_variable(name, standard_name, long_name, units)) then
        dem_varid = varid
        ok = .true.
        return
      end if

      if (dem_varid < 0) dem_varid = varid
    end do

    if (dem_varid >= 0) ok = .true.
  end subroutine find_dem_variable


  subroutine read_axis_1d(ncid, varid, arr, n, ok)
    integer, intent(in) :: ncid, varid
    real(kind=rkind), dimension(:), allocatable, intent(out) :: arr
    integer(kind=ikind), intent(out) :: n
    logical, intent(out) :: ok

    integer :: ncerr, ndims, n_default
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimids

    ok = .false.
    n = 0_ikind

    ncerr = nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids)
    if (ncerr /= nf90_noerr) return

    if (ndims /= 1) return

    ncerr = nf90_inquire_dimension(ncid, dimids(1), len=n_default)
    if (ncerr /= nf90_noerr) return

    n = int(n_default, kind=ikind)

    allocate(arr(n))
    ncerr = nf90_get_var(ncid, varid, arr)
    if (ncerr /= nf90_noerr) then
      if (allocated(arr)) deallocate(arr)
      n = 0_ikind
      return
    end if

    ok = .true.
  end subroutine read_axis_1d


  subroutine try_get_text_att(ncid, varid, attname, value)
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in)  :: attname
    character(len=*), intent(out) :: value

    integer :: ncerr

    value = ""
    ncerr = nf90_get_att(ncid, varid, trim(attname), value)
  end subroutine try_get_text_att


  pure logical function is_latitude_variable(name, standard_name, units, long_name)
    character(len=*), intent(in) :: name, standard_name, units, long_name

    is_latitude_variable = .false.

    if (trim(lower(name)) == "lat") is_latitude_variable = .true.
    if (index(lower(name), "latitude") > 0) is_latitude_variable = .true.
    if (trim(lower(standard_name)) == "latitude") is_latitude_variable = .true.
    if (index(lower(units), "degrees_north") > 0) is_latitude_variable = .true.
    if (index(lower(long_name), "latitude") > 0) is_latitude_variable = .true.
  end function is_latitude_variable


  pure logical function is_longitude_variable(name, standard_name, units, long_name)
    character(len=*), intent(in) :: name, standard_name, units, long_name

    is_longitude_variable = .false.

    if (trim(lower(name)) == "lon") is_longitude_variable = .true.
    if (index(lower(name), "longitude") > 0) is_longitude_variable = .true.
    if (trim(lower(standard_name)) == "longitude") is_longitude_variable = .true.
    if (index(lower(units), "degrees_east") > 0) is_longitude_variable = .true.
    if (index(lower(long_name), "longitude") > 0) is_longitude_variable = .true.
  end function is_longitude_variable


  pure logical function looks_like_dem_variable(name, standard_name, long_name, units)
    character(len=*), intent(in) :: name, standard_name, long_name, units

    looks_like_dem_variable = .false.

    if (index(lower(name), "elev") > 0) looks_like_dem_variable = .true.
    if (index(lower(name), "height") > 0) looks_like_dem_variable = .true.
    if (index(lower(name), "dem") > 0) looks_like_dem_variable = .true.
    if (index(lower(name), "band") > 0) looks_like_dem_variable = .true.
    if (index(lower(standard_name), "height") > 0) looks_like_dem_variable = .true.
    if (index(lower(long_name), "elev") > 0) looks_like_dem_variable = .true.
    if (trim(lower(units)) == "m") looks_like_dem_variable = .true.
    if (index(lower(units), "meter") > 0) looks_like_dem_variable = .true.
  end function looks_like_dem_variable


  pure function lower(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out

    integer(kind=ikind) :: i, c

    do i = 1_ikind, len(s, kind=ikind)
      c = iachar(s(i:i), kind=ikind)
      if (c >= iachar('A', kind=ikind) .and. c <= iachar('Z', kind=ikind)) then
        out(i:i) = achar(c + 32)
      else
        out(i:i) = s(i:i)
      end if
    end do
  end function lower
  
  
  subroutine getmeshalt()
    use typy
    use globals
    use global_objs
    use ncglobvars
    use nctools
    use core_tools


    type(dem_cache) :: dem
    real(kind=rkind) :: z, lat, lon
    integer(kind=ikind) :: i
    logical :: success
    character(len=4096) :: msg

    call dem_open(dem, "drutes.conf/netcdf/dem.nc", success)
    
    if (.not. success) then
      print *, "Unable to open file drutes.conf/netcdf/dem.nc"
      ERROR STOP
    end if
    
    allocate(nodealt(nodes%kolik))
    
    do i=1, nodes%kolik
      call utm2latlong(nodes%data(i,1), nodes%data(i,2), lat, lon)
      call dem_get_altitude(dem, lat, lon, nodealt(i), success)
      if (.not. success) then
		write(msg, *) "W: failed to get dem altitude, check file drutes.conf/netcdf/dem.nc, for node:", i, &
						", node will be deactivated"
        call write_log(msg)
      end if
    end do
    

    call dem_close(dem)
  end subroutine getmeshalt

end module ncdem



