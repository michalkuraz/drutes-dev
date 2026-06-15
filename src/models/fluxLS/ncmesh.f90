module ncmesh
  use netcdf
  use typy
  use ncglobvars
  use global_objs
  use ncdem



  public :: getncmesh

contains


  subroutine getncmesh(ncid, ncnodes, ncelements, success, errmsg)
    integer, intent(in) :: ncid
    type(node), intent(inout) :: ncnodes
    type(element), intent(inout) :: ncelements
    logical, intent(out) :: success
    character(len=*), intent(out) :: errmsg

    type(dem_cache) :: dem

    integer :: ierr
    integer :: lat_varid, lon_varid
    integer :: lat_bnds_varid, lon_bnds_varid
    integer :: lat_dimid, lon_dimid
    integer :: nlat, nlon
    integer :: ilat, ilon
    integer :: jlat, jlon
    integer :: inode, ielem
    integer :: n1, n2, n3, n4

    integer, dimension(1) :: lat_dimids
    integer, dimension(1) :: lon_dimids

    integer(kind=ikind) :: nmissing_dem

    real(kind=rkind), dimension(:,:), allocatable :: lat_bnds
    real(kind=rkind), dimension(:,:), allocatable :: lon_bnds
    real(kind=rkind), dimension(:), allocatable :: lat_corner
    real(kind=rkind), dimension(:), allocatable :: lon_corner

    real(kind=rkind) :: lat0, lon0
    real(kind=rkind) :: x0, y0
    logical :: ok_dem

    success = .false.
    errmsg = "getncmesh: unknown error"

    call cleanup_mesh(ncnodes, ncelements)

    ierr = nf90_inq_varid(ncid, "lat", lat_varid)
    if (ierr /= nf90_noerr) then
      errmsg = "getncmesh: cannot find variable lat: " // trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_inq_varid(ncid, "lon", lon_varid)
    if (ierr /= nf90_noerr) then
      errmsg = "getncmesh: cannot find variable lon: " // trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_inquire_variable(ncid, lat_varid, dimids = lat_dimids)
    if (ierr /= nf90_noerr) then
      errmsg = "getncmesh: cannot inquire lat variable: " // trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_inquire_variable(ncid, lon_varid, dimids = lon_dimids)
    if (ierr /= nf90_noerr) then
      errmsg = "getncmesh: cannot inquire lon variable: " // trim(nf90_strerror(ierr))
      return
    end if

    lat_dimid = lat_dimids(1)
    lon_dimid = lon_dimids(1)

    ierr = nf90_inquire_dimension(ncid, lat_dimid, len = nlat)
    if (ierr /= nf90_noerr) then
      errmsg = "getncmesh: cannot inquire lat dimension: " // trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_inquire_dimension(ncid, lon_dimid, len = nlon)
    if (ierr /= nf90_noerr) then
      errmsg = "getncmesh: cannot inquire lon dimension: " // trim(nf90_strerror(ierr))
      return
    end if

    if (nlat <= 0 .or. nlon <= 0) then
      errmsg = "getncmesh: invalid nlat/nlon"
      return
    end if

    ierr = nf90_inq_varid(ncid, "lat_bnds", lat_bnds_varid)
    if (ierr /= nf90_noerr) then
      errmsg = "getncmesh: cannot find variable lat_bnds: " // trim(nf90_strerror(ierr))
      return
    end if

    ierr = nf90_inq_varid(ncid, "lon_bnds", lon_bnds_varid)
    if (ierr /= nf90_noerr) then
      errmsg = "getncmesh: cannot find variable lon_bnds: " // trim(nf90_strerror(ierr))
      return
    end if

    call read_bounds_2d(ncid, lat_bnds_varid, nlat, lat_bnds, success, errmsg)
    if (.not. success) then
      return
    end if

    call read_bounds_2d(ncid, lon_bnds_varid, nlon, lon_bnds, success, errmsg)
    if (.not. success) then
      call cleanup_bounds(lat_bnds, lon_bnds, lat_corner, lon_corner)
      return
    end if

    call bounds_to_corners(lat_bnds, lat_corner, success, errmsg)
    if (.not. success) then
      call cleanup_bounds(lat_bnds, lon_bnds, lat_corner, lon_corner)
      return
    end if

    call bounds_to_corners(lon_bnds, lon_corner, success, errmsg)
    if (.not. success) then
      call cleanup_bounds(lat_bnds, lon_bnds, lat_corner, lon_corner)
      return
    end if

    if (size(lat_corner) /= nlat + 1) then
      errmsg = "getncmesh: invalid number of latitude corners"
      call cleanup_bounds(lat_bnds, lon_bnds, lat_corner, lon_corner)
      return
    end if

    if (size(lon_corner) /= nlon + 1) then
      errmsg = "getncmesh: invalid number of longitude corners"
      call cleanup_bounds(lat_bnds, lon_bnds, lat_corner, lon_corner)
      return
    end if

    call dem_open(dem, ok_dem)

    if (.not. ok_dem) then
      errmsg = "getncmesh: unable to open DEM file"
      call cleanup_bounds(lat_bnds, lon_bnds, lat_corner, lon_corner)
      return
    end if

    ncnodes%kolik = int((nlat + 1)*(nlon + 1), kind=ikind)
    ncelements%kolik = int(nlat*nlon, kind=ikind)

    allocate(ncnodes%data(ncnodes%kolik, 3))
    allocate(ncelements%data(ncelements%kolik, 4))

    nmissing_dem = 0_ikind

    inode = 0

    do jlat = 1, nlat + 1
      do jlon = 1, nlon + 1

        inode = inode + 1

        lat0 = lat_corner(jlat)
        lon0 = lon_corner(jlon)

        call latlong2utm(lat0, lon0, x0, y0)

        ncnodes%data(inode, 1) = x0
        ncnodes%data(inode, 2) = y0

        call dem_get_altitude(dem, lat0, lon0, ncnodes%data(inode, 3), ok_dem)

        if (.not. ok_dem) then
          ncnodes%data(inode, 3) = -9999.0_rkind
          nmissing_dem = nmissing_dem + 1_ikind
        end if

      end do
    end do

    ielem = 0

    do ilat = 1, nlat
      do ilon = 1, nlon

        ielem = ielem + 1

        n1 = nc_node_id(ilat    , ilon    , nlon)
        n2 = nc_node_id(ilat    , ilon + 1, nlon)
        n3 = nc_node_id(ilat + 1, ilon + 1, nlon)
        n4 = nc_node_id(ilat + 1, ilon    , nlon)

        ncelements%data(ielem, :) = [n1, n2, n3, n4]

      end do
    end do

    call dem_close(dem)
    call cleanup_bounds(lat_bnds, lon_bnds, lat_corner, lon_corner)

    if (nmissing_dem > 0_ikind) then
      write(*,*) "getncmesh: DEM altitude missing for ", nmissing_dem, &
                 " NetCDF mesh nodes; filled with -9999."
    end if

    success = .true.
    errmsg = "everything ok"

  end subroutine getncmesh


  subroutine read_bounds_2d(ncid, varid, ncell, bounds, success, errmsg)
    integer, intent(in) :: ncid
    integer, intent(in) :: varid
    integer, intent(in) :: ncell
    real(kind=rkind), dimension(:,:), allocatable, intent(out) :: bounds
    logical, intent(out) :: success
    character(len=*), intent(out) :: errmsg

    integer :: ierr
    real(kind=rkind), dimension(:,:), allocatable :: tmp

    success = .false.
    errmsg = "read_bounds_2d: unknown error"

    if (allocated(bounds)) then
      deallocate(bounds)
    end if

    if (ncell <= 0) then
      errmsg = "read_bounds_2d: invalid ncell"
      return
    end if

    allocate(tmp(2, ncell))

    ierr = nf90_get_var(ncid, varid, tmp)

    if (ierr /= nf90_noerr) then
      errmsg = "read_bounds_2d: cannot read bounds variable: " // trim(nf90_strerror(ierr))
      if (allocated(tmp)) then
        deallocate(tmp)
      end if
      return
    end if

    allocate(bounds(ncell, 2))

    bounds = transpose(tmp)

    deallocate(tmp)

    success = .true.
    errmsg = "everything ok"

  end subroutine read_bounds_2d


  subroutine bounds_to_corners(bounds, corners, success, errmsg)
    real(kind=rkind), dimension(:,:), intent(in) :: bounds
    real(kind=rkind), dimension(:), allocatable, intent(out) :: corners
    logical, intent(out) :: success
    character(len=*), intent(out) :: errmsg

    integer :: ncell
    integer :: i

    success = .false.
    errmsg = "bounds_to_corners: unknown error"

    ncell = size(bounds, 1)

    if (ncell <= 0) then
      errmsg = "bounds_to_corners: invalid number of cells"
      return
    end if

    if (size(bounds, 2) /= 2) then
      errmsg = "bounds_to_corners: bounds second dimension must be 2"
      return
    end if

    if (allocated(corners)) then
      deallocate(corners)
    end if

    allocate(corners(ncell + 1))

    do i = 1, ncell
      corners(i) = min(bounds(i, 1), bounds(i, 2))
    end do

    corners(ncell + 1) = maxval(bounds)

    call sort_real_array(corners)

    success = .true.
    errmsg = "everything ok"

  end subroutine bounds_to_corners


  subroutine sort_real_array(arr)
    real(kind=rkind), dimension(:), intent(inout) :: arr

    integer :: i, j, n
    real(kind=rkind) :: tmp

    n = size(arr)

    do i = 1, n - 1
      do j = i + 1, n
        if (arr(j) < arr(i)) then
          tmp = arr(i)
          arr(i) = arr(j)
          arr(j) = tmp
        end if
      end do
    end do

  end subroutine sort_real_array


  integer function nc_node_id(jlat, jlon, nlon) result(id)
    integer, intent(in) :: jlat
    integer, intent(in) :: jlon
    integer, intent(in) :: nlon

    id = (jlat - 1)*(nlon + 1) + jlon

  end function nc_node_id


  subroutine latlong2utm(latit, longit, x, y)
    real(kind=rkind), intent(in)  :: latit, longit
    real(kind=rkind), intent(out) :: x, y

    real(kind=rkind) :: a, f, e2, ep2, k0
    real(kind=rkind) :: pi
    real(kind=rkind) :: lat_rad, lon_rad, lambda0
    real(kind=rkind) :: n, t, c, aa, meridian_arc
    logical :: is_south

    pi = 4.0_rkind * atan(1.0_rkind)

    a  = 6378137.0_rkind
    f  = 1.0_rkind / 298.257223563_rkind
    e2 = f * (2.0_rkind - f)
    ep2 = e2 / (1.0_rkind - e2)
    k0 = 0.9996_rkind

    lat_rad = latit  * pi / 180.0_rkind
    lon_rad = longit * pi / 180.0_rkind

    lambda0 = (geograzone * 6.0_rkind - 183.0_rkind) * pi / 180.0_rkind

    n  = a / sqrt(1.0_rkind - e2 * sin(lat_rad)**2)
    t  = tan(lat_rad)**2
    c  = ep2 * cos(lat_rad)**2
    aa = cos(lat_rad) * (lon_rad - lambda0)

    meridian_arc = a * ( &
        (1.0_rkind - e2/4.0_rkind - 3.0_rkind*e2**2/64.0_rkind - 5.0_rkind*e2**3/256.0_rkind) * lat_rad &
      - (3.0_rkind*e2/8.0_rkind + 3.0_rkind*e2**2/32.0_rkind + 45.0_rkind*e2**3/1024.0_rkind) * sin(2.0_rkind*lat_rad) &
      + (15.0_rkind*e2**2/256.0_rkind + 45.0_rkind*e2**3/1024.0_rkind) * sin(4.0_rkind*lat_rad) &
      - (35.0_rkind*e2**3/3072.0_rkind) * sin(6.0_rkind*lat_rad) )

    x = k0 * n * ( &
          aa &
        + (1.0_rkind - t + c) * aa**3 / 6.0_rkind &
        + (5.0_rkind - 18.0_rkind*t + t**2 + 72.0_rkind*c - 58.0_rkind*ep2) * aa**5 / 120.0_rkind )

    x = x + 500000.0_rkind

    y = k0 * ( &
          meridian_arc &
        + n * tan(lat_rad) * ( &
            aa**2 / 2.0_rkind &
          + (5.0_rkind - t + 9.0_rkind*c + 4.0_rkind*c**2) * aa**4 / 24.0_rkind &
          + (61.0_rkind - 58.0_rkind*t + t**2 + 600.0_rkind*c - 330.0_rkind*ep2) * aa**6 / 720.0_rkind ) )

    is_south = latit < 0.0_rkind

    if (is_south) then
      y = y + 10000000.0_rkind
    end if

  end subroutine latlong2utm


  subroutine cleanup_bounds(lat_bnds, lon_bnds, lat_corner, lon_corner)
    real(kind=rkind), dimension(:,:), allocatable, intent(inout) :: lat_bnds
    real(kind=rkind), dimension(:,:), allocatable, intent(inout) :: lon_bnds
    real(kind=rkind), dimension(:), allocatable, intent(inout) :: lat_corner
    real(kind=rkind), dimension(:), allocatable, intent(inout) :: lon_corner

    if (allocated(lat_bnds)) then
      deallocate(lat_bnds)
    end if

    if (allocated(lon_bnds)) then
      deallocate(lon_bnds)
    end if

    if (allocated(lat_corner)) then
      deallocate(lat_corner)
    end if

    if (allocated(lon_corner)) then
      deallocate(lon_corner)
    end if

  end subroutine cleanup_bounds


  subroutine cleanup_mesh(ncnodes, ncelements)
    type(node), intent(inout) :: ncnodes
    type(element), intent(inout) :: ncelements

    if (allocated(ncnodes%data)) then
      deallocate(ncnodes%data)
    end if

    if (allocated(ncelements%data)) then
      deallocate(ncelements%data)
    end if

    ncnodes%kolik = 0_ikind
    ncelements%kolik = 0_ikind

  end subroutine cleanup_mesh
  
  subroutine point_in_quad(xp, yp, xq, yq, inside)
    real(kind=rkind), intent(in) :: xp, yp
    real(kind=rkind), dimension(4), intent(in) :: xq, yq
    logical, intent(out) :: inside

    logical :: inside1
    logical :: inside2

    call point_in_triangle(xp, yp, xq(1), yq(1), xq(2), yq(2), xq(3), yq(3), inside1)
    call point_in_triangle(xp, yp, xq(1), yq(1), xq(3), yq(3), xq(4), yq(4), inside2)

    inside = inside1 .or. inside2

  end subroutine point_in_quad


  subroutine point_in_triangle(xp, yp, x1, y1, x2, y2, x3, y3, inside)
    real(kind=rkind), intent(in) :: xp, yp
    real(kind=rkind), intent(in) :: x1, y1
    real(kind=rkind), intent(in) :: x2, y2
    real(kind=rkind), intent(in) :: x3, y3
    logical, intent(out) :: inside

    real(kind=rkind) :: d1, d2, d3
    logical :: has_neg
    logical :: has_pos
    real(kind=rkind), parameter :: eps = 1.0e-10_rkind

    d1 = triangle_sign(xp, yp, x1, y1, x2, y2)
    d2 = triangle_sign(xp, yp, x2, y2, x3, y3)
    d3 = triangle_sign(xp, yp, x3, y3, x1, y1)

    has_neg = (d1 < -eps) .or. (d2 < -eps) .or. (d3 < -eps)
    has_pos = (d1 >  eps) .or. (d2 >  eps) .or. (d3 >  eps)

    inside = .not. (has_neg .and. has_pos)

  end subroutine point_in_triangle


  pure function triangle_sign(xp, yp, x1, y1, x2, y2) result(s)
    real(kind=rkind), intent(in) :: xp, yp
    real(kind=rkind), intent(in) :: x1, y1
    real(kind=rkind), intent(in) :: x2, y2
    real(kind=rkind) :: s

    s = (xp - x2)*(y1 - y2) - (x1 - x2)*(yp - y2)

  end function triangle_sign

end module ncmesh
