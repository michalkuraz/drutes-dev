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
