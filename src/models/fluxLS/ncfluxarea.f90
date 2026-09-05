module ncfluxarea
  use typy
  use ncglobvars
  use nctools

  public :: ncflux_cell_area
  public :: ncflux_cell_area_xy
  public :: ncflux_active_width

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


  ! Return the mapped hydrological cell width perpendicular to an element's flow direction.
  function ncflux_active_width(element_number) result(width)
    integer(kind=ikind), intent(in) :: element_number
    real(kind=rkind) :: width

    integer(kind=ikind) :: grid_element, i, node_number
    real(kind=rkind), dimension(2) :: direction, normal
    real(kind=rkind), dimension(4) :: projection
    real(kind=rkind) :: direction_length

    width = 0.0_rkind

    if (.not. allocated(el2ncgrid)) return
    if (.not. allocated(ncfluxdata%fluxvct)) return
    if (.not. allocated(ncelements%data)) return
    if (.not. allocated(ncnodes%data)) return

    if (element_number < 1_ikind .or. element_number > size(el2ncgrid, kind=ikind)) return
    if (element_number > size(ncfluxdata%fluxvct, 1, kind=ikind)) return

    grid_element = el2ncgrid(element_number)
    if (grid_element < 1_ikind .or. grid_element > ncelements%kolik) return

    direction = ncfluxdata%fluxvct(element_number,:)
    direction_length = norm2(direction)
    if (direction_length <= 10.0_rkind*epsilon(1.0_rkind)) return

    direction = direction/direction_length
    normal = (/ -direction(2), direction(1) /)

    do i = 1_ikind, 4_ikind
      node_number = ncelements%data(grid_element,i)
      if (node_number < 1_ikind .or. node_number > ncnodes%kolik) then
        width = 0.0_rkind
        return
      end if
      projection(i) = dot_product(ncnodes%data(node_number,1:2), normal)
    end do

    width = maxval(projection) - minval(projection)
  end function ncflux_active_width


  subroutine find_cell_from_bounds(bounds, value, idx, ok)
    use typy
    implicit none

    real(kind=rkind), intent(in) :: bounds(:,:)
    real(kind=rkind), intent(in) :: value
    integer(kind=ikind), intent(out) :: idx
    logical, intent(out) :: ok

    integer :: i
    real(kind=rkind) :: lo, hi
    real(kind=rkind), parameter :: eps = 1.0e-10_rkind

    ok = .false.
    idx = -1_ikind

    do i = 1, size(bounds,1)
      lo = min(bounds(i,1), bounds(i,2))
      hi = max(bounds(i,1), bounds(i,2))

      if (value >= lo - eps .and. value <= hi + eps) then
        idx = int(i, kind=ikind)
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

