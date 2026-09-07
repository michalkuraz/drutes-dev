module ncfluxarea
  use typy
  use ncglobvars
  use nctools
  use ncwidth_geometry, only: river_contact_width
  use netcdfflux, only: is_missing_flux
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  public :: ncflux_cell_area
  public :: ncflux_cell_area_xy
  public :: ncflux_active_width
  public :: ncflux_prepare_widths

  ! Geometry and the active channel mask are fixed during a simulation.
  real(kind=rkind), allocatable, private :: active_widths(:)

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


  ! Cache effective widths after mapping elements and assigning flow directions.
  ! Use the same initial NetCDF slice/Qmin as channel activation, excluding dry
  ! cells: a valid zero discharge is not an active river contact even if Qmin=0.
  subroutine ncflux_prepare_widths(ok, errmsg, zero_width_count)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: errmsg
    integer(kind=ikind), intent(out), optional :: zero_width_count
    logical, allocatable :: active_cells(:)
    real(kind=rkind), allocatable :: vertices(:,:,:)
    real(kind=rkind) :: neighbours(4,2,8), q
    integer(kind=ikind) :: cell, el, i, inode, ilat, ilon, row, col, adjacent, count

    ok = .false.
    if (present(zero_width_count)) zero_width_count = 0
    errmsg = "ncflux_prepare_widths: initial grid, discharge slice or directions are missing"
    if (allocated(active_widths)) deallocate(active_widths)
    if (.not. ncfluxdata%initialized) return
    if (.not. ncfluxdata%slice_loaded) return
    if (.not. allocated(ncfluxdata%qslice)) return
    if (.not. allocated(ncfluxdata%activeel)) return
    if (.not. allocated(ncfluxdata%fluxvct)) return
    if (.not. allocated(el2ncgrid)) return
    if (.not. allocated(ncelements%data)) return
    if (.not. allocated(ncnodes%data)) return
    if (ncfluxdata%nlat < 1 .or. ncfluxdata%nlon < 1) return
    errmsg = "ncflux_prepare_widths: incompatible hydrological grid or FE array dimensions"
    if (ncelements%kolik /= ncfluxdata%nlat*ncfluxdata%nlon) return
    if (size(ncelements%data,1,kind=ikind) /= ncelements%kolik) return
    if (size(ncelements%data,2) /= 4 .or. size(ncnodes%data,2) < 2) return
    if (size(ncfluxdata%qslice,1,kind=ikind) /= ncfluxdata%nlat) return
    if (size(ncfluxdata%qslice,2,kind=ikind) /= ncfluxdata%nlon) return
    if (size(ncfluxdata%activeel) /= size(el2ncgrid)) return
    if (size(ncfluxdata%fluxvct,1) /= size(el2ncgrid)) return
    if (size(ncfluxdata%fluxvct,2) /= 2) return

    allocate(vertices(4,2,ncelements%kolik), active_cells(ncelements%kolik))
    do cell = 1, ncelements%kolik
      ! getncmesh numbers cells with longitude varying fastest; it preserves
      ! the NetCDF axis order, including descending latitude/longitude axes.
      ilat = (cell-1)/ncfluxdata%nlon + 1
      ilon = mod(cell-1,ncfluxdata%nlon) + 1
      q = ncfluxdata%qslice(ilat,ilon)
      active_cells(cell) = .false.
      ! Fortran need not short-circuit logical expressions; do not compare NaN
      ! discharge values or pass them to the fill-value comparison.
      if (ieee_is_finite(q)) then
        if (.not. is_missing_flux(q)) active_cells(cell) = q > 0.0_rkind .and. q >= Qmin
      end if
      do i = 1, 4
        inode = ncelements%data(cell,i)
        errmsg = "ncflux_prepare_widths: invalid hydrological node index"
        if (inode < 1 .or. inode > size(ncnodes%data,1,kind=ikind)) return
        vertices(i,:,cell) = ncnodes%data(inode,1:2)
      end do
    end do

    allocate(active_widths(size(el2ncgrid)))
    active_widths = 0.0_rkind
    do el = 1, size(el2ncgrid,kind=ikind)
      if (.not. ncfluxdata%activeel(el)) cycle
      cell = el2ncgrid(el)
      if (cell < 1 .or. cell > ncelements%kolik) cycle
      if (.not. active_cells(cell)) cycle
      ilat = (cell-1)/ncfluxdata%nlon + 1
      ilon = mod(cell-1,ncfluxdata%nlon) + 1
      count = 0
      do row = max(1_ikind,ilat-1), min(ncfluxdata%nlat,ilat+1)
        do col = max(1_ikind,ilon-1), min(ncfluxdata%nlon,ilon+1)
          adjacent = (row-1)*ncfluxdata%nlon + col
          if (adjacent == cell) cycle
          if (.not. active_cells(adjacent)) cycle
          count = count + 1
          neighbours(:,:,count) = vertices(:,:,adjacent)
        end do
      end do
      active_widths(el) = river_contact_width(vertices(:,:,cell), &
        ncfluxdata%fluxvct(el,:), neighbours(:,:,1:count))
      if (present(zero_width_count)) then
        if (active_widths(el) <= 0.0_rkind) zero_width_count = zero_width_count + 1
      end if
    end do
    ok = .true.
    errmsg = "ncflux_prepare_widths: everything ok"
  end subroutine ncflux_prepare_widths


  ! Effective transverse width limited by the active inlet/outlet contacts.
  ! Call ncflux_prepare_widths again if the mask, mapping or directions change.
  function ncflux_active_width(element_number) result(width)
    integer(kind=ikind), intent(in) :: element_number
    real(kind=rkind) :: width

    width = 0.0_rkind
    if (.not. allocated(active_widths)) return
    if (.not. ncfluxdata%initialized) return
    if (element_number < 1 .or. element_number > size(active_widths,kind=ikind)) return
    width = active_widths(element_number)
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
