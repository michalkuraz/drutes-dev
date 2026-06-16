module nctools
  public :: utm2latlong
  contains
  
  subroutine utm2latlong(x, y, latit, longit)
    use typy
    use ncglobvars
    
    real(kind=rkind), intent(in)  :: x, y
    real(kind=rkind), intent(out) :: latit, longit

    real(kind=rkind) :: a, f, e2, ep2, k0
    real(kind=rkind) :: x_rel, y_rel
    real(kind=rkind) :: mu, e1, phi1
    real(kind=rkind) :: n1, r1, t1, c1, d
    real(kind=rkind) :: lambda0
    real(kind=rkind) :: pi
    logical          :: is_south

    pi = 4.0_rkind*atan(1.0_rkind)

    ! WGS-84 ellipsoid
    a  = 6378137.0_rkind
    f  = 1.0_rkind/298.257223563_rkind
    e2 = f*(2.0_rkind - f)
    ep2 = e2/(1.0_rkind - e2)
    k0 = 0.9996_rkind

    ! Central meridian of UTM zone (radians)
    lambda0 = (geograzone*6.0_rkind - 183.0_rkind) * pi/180.0_rkind

    ! Remove false easting
    x_rel = x - 500000.0_rkind

    ! Automatic hemisphere detection
    is_south = (y >= 10000000.0_rkind)

    if (is_south) then
        y_rel = (y - 10000000.0_rkind) / k0
    else
        y_rel = y / k0
    endif

    ! Footpoint latitude
    e1 = (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2))

    mu = y_rel / (a*(1 - e2/4.0_rkind - 3*e2**2/64.0_rkind - 5*e2**3/256.0_rkind))

    phi1 = mu + (3*e1/2.0_rkind - 27*e1**3/32.0_rkind)*sin(2*mu) &
              + (21*e1**2/16.0_rkind - 55*e1**4/32.0_rkind)*sin(4*mu) &
              + (151*e1**3/96.0_rkind)*sin(6*mu) &
              + (1097*e1**4/512.0_rkind)*sin(8*mu)

    ! Intermediate values
    n1 = a / sqrt(1 - e2*sin(phi1)**2)
    r1 = a*(1-e2) / (1 - e2*sin(phi1)**2)**1.5_rkind
    t1 = tan(phi1)**2
    c1 = ep2 * cos(phi1)**2
    d  = x_rel / (n1*k0)

    ! Latitude (radians)
    latit = phi1 - (n1*tan(phi1)/r1) * ( &
          d**2/2 - &
          (5 + 3*t1 + 10*c1 - 4*c1**2 - 9*ep2)*d**4/24 + &
          (61 + 90*t1 + 298*c1 + 45*t1**2 - 252*ep2 - 3*c1**2)*d**6/720.0_rkind )

    ! Longitude (radians)
    longit = lambda0 + ( &
          d - (1 + 2*t1 + c1)*d**3/6.0_rkind + &
         (5 - 2*c1 + 28*t1 - 3*c1**2 + 8*ep2 + 24*t1**2)*d**5/120.0_rkind ) / cos(phi1)

    ! Convert to degrees
    latit = latit * 180.0_rkind / pi
    longit = longit * 180.0_rkind / pi

  end subroutine utm2latlong
  
  
  subroutine terrain_slopes()
    use typy
    use globals
    use global_objs
    use ncglobvars
    use core_tools

    
    integer(kind=ikind) :: el, nd
    real(kind=rkind), dimension(3,3) :: pts
    logical :: elfine
    
    allocate(elslopes(elements%kolik, 2))

    elslopes = 0.0_rkind    
    do el=1, elements%kolik
      elfine = .true.
      do nd = 1,3
        pts(nd,1:2) = nodes%data(elements%data(el,nd),:)
        pts(nd,3) = nodealt(elements%data(el,nd))
        if (int(pts(nd,3)) == missing) then
          nodes%edge(elements%data(el,:)) = addedbc
          elfine = .false.
          EXIT
        end if
          
      end do

      if (elfine)  call plane_derivative(pts(1,:), pts(2,:), pts(3,:), elslopes(el,1), elslopes(el,2))
    end do
    
	end subroutine terrain_slopes
  
  subroutine ncelslope()
    use typy
    use global_objs
    use ncglobvars
    use globals


    integer(kind=ikind) :: ielem
    integer(kind=ikind) :: i
    integer(kind=ikind) :: inode

    integer(kind=ikind), dimension(4) :: nds

    real(kind=rkind), dimension(4) :: x
    real(kind=rkind), dimension(4) :: y
    real(kind=rkind), dimension(4) :: z

    real(kind=rkind) :: xbar, ybar, zbar
    real(kind=rkind) :: xx, yy, xy
    real(kind=rkind) :: xz, yz
    real(kind=rkind) :: det
    real(kind=rkind) :: dzdx, dzdy
    logical :: valid_element

    if (allocated(ncelements%ders)) then
      deallocate(ncelements%ders)
    end if

    allocate(ncelements%ders(ncelements%kolik, drutes_config%dimen, 1))

    ncelements%ders = missing

    do ielem = 1, ncelements%kolik

      nds(:) = ncelements%data(ielem, 1:4)

      valid_element = .true.

      do i = 1, 4

        inode = nds(i)

        if (inode <= 0_ikind .or. inode > ncnodes%kolik) then
          valid_element = .false.
          exit
        end if

        x(i) = ncnodes%data(inode, 1)
        y(i) = ncnodes%data(inode, 2)
        z(i) = ncnodes%data(inode, 3)

        if (abs(z(i) - missing) < 1.0e-8_rkind) then
          valid_element = .false.
          exit
        end if

      end do

      if (.not. valid_element) then
        ncelements%ders(ielem, 1, 1) = missing
        ncelements%ders(ielem, 2, 1) = missing
        cycle
      end if

      xbar = sum(x) / 4.0_rkind
      ybar = sum(y) / 4.0_rkind
      zbar = sum(z) / 4.0_rkind

      xx = 0.0_rkind
      yy = 0.0_rkind
      xy = 0.0_rkind
      xz = 0.0_rkind
      yz = 0.0_rkind

      do i = 1_ikind, 4_ikind
        xx = xx + (x(i) - xbar) * (x(i) - xbar)
        yy = yy + (y(i) - ybar) * (y(i) - ybar)
        xy = xy + (x(i) - xbar) * (y(i) - ybar)

        xz = xz + (x(i) - xbar) * (z(i) - zbar)
        yz = yz + (y(i) - ybar) * (z(i) - zbar)
      end do

      det = xx*yy - xy*xy

      if (abs(det) < epsilon(det)) then
        ncelements%ders(ielem, 1, 1) = missing
        ncelements%ders(ielem, 2, 1) = missing
        cycle
      end if

      dzdx = (xz*yy - yz*xy) / det
      dzdy = (yz*xx - xz*xy) / det

      ncelements%ders(ielem, 1, 1) = dzdx
      ncelements%ders(ielem, 2, 1) = dzdy

    end do

  end subroutine ncelslope


end module nctools



