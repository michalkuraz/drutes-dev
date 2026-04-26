module ncfluxarea
  use typy
  use ncglobvars
  use nctools

  public :: ncflux_cell_area
  public :: ncflux_cell_area_xy

contains

  function ncflux_cell_area(ilat, ilon) result(area)
    integer(kind=ikind), intent(in) :: ilat, ilon
    real(kind=rkind) :: area

    real(kind=rkind), parameter :: radius_earth = 6371000.0_rkind
    real(kind=rkind) :: pi
    real(kind=rkind) :: lat1, lat2, lon1, lon2
    real(kind=rkind) :: phi1, phi2, lam1, lam2

    area = -1.0_rkind
    pi = 4.0_rkind*atan(1.0_rkind)

    if (.not. ncfluxdata%initialized) return
    if (.not. ncfluxdata%has_bounds) return

    if (ilat < 1_ikind .or. ilat > ncfluxdata%nlat) return
    if (ilon < 1_ikind .or. ilon > ncfluxdata%nlon) return

    lat1 = ncfluxdata%lat_bnds(ilat,1)
    lat2 = ncfluxdata%lat_bnds(ilat,2)
    lon1 = ncfluxdata%lon_bnds(ilon,1)
    lon2 = ncfluxdata%lon_bnds(ilon,2)

    phi1 = lat1*pi/180.0_rkind
    phi2 = lat2*pi/180.0_rkind
    lam1 = lon1*pi/180.0_rkind
    lam2 = lon2*pi/180.0_rkind

    area = radius_earth**2 * abs(lam2 - lam1) * abs(sin(phi2) - sin(phi1))
  end function ncflux_cell_area


  subroutine ncflux_cell_area_xy(x, y, area, ok, errmsg)
    real(kind=rkind), intent(in)  :: x, y
    real(kind=rkind), intent(out) :: area
    logical, intent(out)          :: ok
    character(len=*), intent(out) :: errmsg

    real(kind=rkind) :: lat0, lon0, lon_adj
    integer(kind=ikind) :: ilat, ilon
    logical :: ok_local

    ok = .false.
    area = -1.0_rkind
    errmsg = "ncflux_cell_area_xy: unknown error"

    if (.not. ncfluxdata%initialized) then
      errmsg = "ncflux_cell_area_xy: ncfluxdata not initialized"
      return
    end if

    if (.not. ncfluxdata%has_bounds) then
      errmsg = "ncflux_cell_area_xy: NetCDF file has no lat_bnds/lon_bnds"
      return
    end if

    call utm2latlong(x, y, lat0, lon0)

    lon_adj = adjust_longitude_to_grid(ncfluxdata%lon, lon0)

    call find_cell_from_bounds(ncfluxdata%lat_bnds, lat0, ilat, ok_local)
    if (.not. ok_local) then
      errmsg = "ncflux_cell_area_xy: latitude outside NetCDF bounds"
      return
    end if

    call find_cell_from_bounds(ncfluxdata%lon_bnds, lon_adj, ilon, ok_local)
    if (.not. ok_local) then
      errmsg = "ncflux_cell_area_xy: longitude outside NetCDF bounds"
      return
    end if

    area = ncflux_cell_area(ilat, ilon)

    if (area <= 0.0_rkind) then
      errmsg = "ncflux_cell_area_xy: failed to compute cell area"
      ok = .false.
      return
    end if

    ok = .true.
    errmsg = "everything ok"
  end subroutine ncflux_cell_area_xy


  subroutine find_cell_from_bounds(bounds, x, idx, ok)
    real(kind=rkind), dimension(:,:), intent(in) :: bounds
    real(kind=rkind), intent(in) :: x
    integer(kind=ikind), intent(out) :: idx
    logical, intent(out) :: ok

    integer(kind=ikind) :: i, n
    real(kind=rkind) :: b1, b2

    ok = .false.
    idx = -1_ikind

    n = size(bounds, 1, kind=ikind)

    do i = 1_ikind, n
      b1 = bounds(i,1)
      b2 = bounds(i,2)

      if (x >= min(b1,b2) .and. x <= max(b1,b2)) then
        idx = i
        ok = .true.
        return
      end if
    end do
  end subroutine find_cell_from_bounds


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

end module ncfluxarea



module ncfluxarea_fast
  use typy
  use ncglobvars
  use nctools
  implicit none
  private

  public :: ncfluxarea_fast_init
  public :: ncflux_cell_area
  public :: ncflux_cell_area_xy_fast
  public :: ncflux_cell_index_xy_fast

contains

  subroutine ncfluxarea_fast_init(ok, errmsg)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: errmsg

    ok = .false.
    errmsg = "ncfluxarea_fast_init: unknown error"

    if (.not. ncfluxdata%initialized) then
      errmsg = "ncfluxarea_fast_init: ncfluxdata not initialized"
      return
    end if

    if (.not. ncfluxdata%has_bounds) then
      errmsg = "ncfluxarea_fast_init: NetCDF file has no lat_bnds/lon_bnds"
      return
    end if

    if (ncfluxdata%nlat < 1_ikind .or. ncfluxdata%nlon < 1_ikind) then
      errmsg = "ncfluxarea_fast_init: invalid NetCDF grid size"
      return
    end if

    ncfluxdata%lat_b0 = ncfluxdata%lat_bnds(1,1)
    ncfluxdata%lon_b0 = ncfluxdata%lon_bnds(1,1)

    if (ncfluxdata%nlat > 1_ikind) then
      ncfluxdata%dlat_b = ncfluxdata%lat_bnds(2,1) - ncfluxdata%lat_bnds(1,1)
    else
      ncfluxdata%dlat_b = ncfluxdata%lat_bnds(1,2) - ncfluxdata%lat_bnds(1,1)
    end if

    if (ncfluxdata%nlon > 1_ikind) then
      ncfluxdata%dlon_b = ncfluxdata%lon_bnds(2,1) - ncfluxdata%lon_bnds(1,1)
    else
      ncfluxdata%dlon_b = ncfluxdata%lon_bnds(1,2) - ncfluxdata%lon_bnds(1,1)
    end if

    if (abs(ncfluxdata%dlat_b) <= epsilon(1.0_rkind)) then
      errmsg = "ncfluxarea_fast_init: zero latitude spacing"
      return
    end if

    if (abs(ncfluxdata%dlon_b) <= epsilon(1.0_rkind)) then
      errmsg = "ncfluxarea_fast_init: zero longitude spacing"
      return
    end if

    ncfluxdata%bounds_index_ready = .true.

    ok = .true.
    errmsg = "everything ok"
  end subroutine ncfluxarea_fast_init


  subroutine ncflux_cell_index_xy_fast(x, y, ilat, ilon, ok, errmsg)
    real(kind=rkind), intent(in) :: x, y
    integer(kind=ikind), intent(out) :: ilat, ilon
    logical, intent(out) :: ok
    character(len=*), intent(out) :: errmsg

    real(kind=rkind) :: lat0, lon0, lon_adj

    ok = .false.
    errmsg = "ncflux_cell_index_xy_fast: unknown error"
    ilat = -1_ikind
    ilon = -1_ikind

    if (.not. ncfluxdata%initialized) then
      errmsg = "ncflux_cell_index_xy_fast: ncfluxdata not initialized"
      return
    end if

    if (.not. ncfluxdata%has_bounds) then
      errmsg = "ncflux_cell_index_xy_fast: NetCDF file has no lat_bnds/lon_bnds"
      return
    end if

    if (.not. ncfluxdata%bounds_index_ready) then
      errmsg = "ncflux_cell_index_xy_fast: fast bounds index not initialized"
      return
    end if

    call utm2latlong(x, y, lat0, lon0)
    lon_adj = adjust_longitude_to_grid_fast(ncfluxdata%lon, lon0)

    call index_from_regular_bounds( &
      lat0, ncfluxdata%lat_bnds, ncfluxdata%lat_b0, ncfluxdata%dlat_b, ilat, ok)

    if (.not. ok) then
      errmsg = "ncflux_cell_index_xy_fast: latitude outside NetCDF bounds"
      return
    end if

    call index_from_regular_bounds( &
      lon_adj, ncfluxdata%lon_bnds, ncfluxdata%lon_b0, ncfluxdata%dlon_b, ilon, ok)

    if (.not. ok) then
      errmsg = "ncflux_cell_index_xy_fast: longitude outside NetCDF bounds"
      return
    end if

    ok = .true.
    errmsg = "everything ok"
  end subroutine ncflux_cell_index_xy_fast


  subroutine ncflux_cell_area_xy_fast(x, y, area, ok, errmsg)
    real(kind=rkind), intent(in) :: x, y
    real(kind=rkind), intent(out) :: area
    logical, intent(out) :: ok
    character(len=*), intent(out) :: errmsg

    integer(kind=ikind) :: ilat, ilon

    ok = .false.
    area = -1.0_rkind
    errmsg = "ncflux_cell_area_xy_fast: unknown error"

    call ncflux_cell_index_xy_fast(x, y, ilat, ilon, ok, errmsg)
    if (.not. ok) return

    area = ncflux_cell_area(ilat, ilon)

    if (area <= 0.0_rkind) then
      ok = .false.
      errmsg = "ncflux_cell_area_xy_fast: failed to compute cell area"
      return
    end if

    ok = .true.
    errmsg = "everything ok"
  end subroutine ncflux_cell_area_xy_fast


  function ncflux_cell_area(ilat, ilon) result(area)
    integer(kind=ikind), intent(in) :: ilat, ilon
    real(kind=rkind) :: area

    real(kind=rkind), parameter :: radius_earth = 6371000.0_rkind
    real(kind=rkind) :: pi
    real(kind=rkind) :: lat1, lat2, lon1, lon2
    real(kind=rkind) :: phi1, phi2, lam1, lam2

    area = -1.0_rkind
    pi = 4.0_rkind*atan(1.0_rkind)

    if (.not. ncfluxdata%initialized) return
    if (.not. ncfluxdata%has_bounds) return

    if (ilat < 1_ikind .or. ilat > ncfluxdata%nlat) return
    if (ilon < 1_ikind .or. ilon > ncfluxdata%nlon) return

    lat1 = ncfluxdata%lat_bnds(ilat,1)
    lat2 = ncfluxdata%lat_bnds(ilat,2)
    lon1 = ncfluxdata%lon_bnds(ilon,1)
    lon2 = ncfluxdata%lon_bnds(ilon,2)

    phi1 = lat1*pi/180.0_rkind
    phi2 = lat2*pi/180.0_rkind
    lam1 = lon1*pi/180.0_rkind
    lam2 = lon2*pi/180.0_rkind

    area = radius_earth**2 * abs(lam2 - lam1) * abs(sin(phi2) - sin(phi1))
  end function ncflux_cell_area


  subroutine index_from_regular_bounds(x, bounds, b0, db, idx, ok)
    real(kind=rkind), intent(in) :: x
    real(kind=rkind), dimension(:,:), intent(in) :: bounds
    real(kind=rkind), intent(in) :: b0, db
    integer(kind=ikind), intent(out) :: idx
    logical, intent(out) :: ok

    integer(kind=ikind) :: n
    integer(kind=ikind) :: guess

    ok = .false.
    idx = -1_ikind

    n = size(bounds, 1, kind=ikind)

    if (x < minval(bounds) .or. x > maxval(bounds)) return

    guess = int(floor((x - b0)/db), kind=ikind) + 1_ikind

    if (guess < 1_ikind) guess = 1_ikind
    if (guess > n) guess = n

    if (inside_bounds(x, bounds(guess,1), bounds(guess,2))) then
      idx = guess
      ok = .true.
      return
    end if

    if (guess > 1_ikind) then
      if (inside_bounds(x, bounds(guess-1_ikind,1), bounds(guess-1_ikind,2))) then
        idx = guess - 1_ikind
        ok = .true.
        return
      end if
    end if

    if (guess < n) then
      if (inside_bounds(x, bounds(guess+1_ikind,1), bounds(guess+1_ikind,2))) then
        idx = guess + 1_ikind
        ok = .true.
        return
      end if
    end if
  end subroutine index_from_regular_bounds


  pure logical function inside_bounds(x, b1, b2)
    real(kind=rkind), intent(in) :: x, b1, b2

    inside_bounds = (x >= min(b1,b2) .and. x <= max(b1,b2))
  end function inside_bounds


  pure function adjust_longitude_to_grid_fast(lon_arr, lon0) result(x)
    real(kind=rkind), dimension(:), intent(in) :: lon_arr
    real(kind=rkind), intent(in) :: lon0
    real(kind=rkind) :: x
    real(kind=rkind) :: lonmin, lonmax

    lonmin = minval(lon_arr)
    lonmax = maxval(lon_arr)
    x = lon0

    if (lonmin >= 0.0_rkind .and. x < 0.0_rkind) x = x + 360.0_rkind
    if (lonmax <= 180.0_rkind .and. x > 180.0_rkind) x = x - 360.0_rkind
  end function adjust_longitude_to_grid_fast

end module ncfluxarea_fast
