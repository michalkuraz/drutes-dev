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
  
  
  subroutine getncalt(xy, alt) 
    use typy
    use netcdf
    use ncglobvars
    
    real(kind=rkind), dimension(2), intent(in) :: xy
    real(kind=rkind), intent(out) :: alt
    
    real(kind=rkind) :: latit, longit
    
    call utm2latlong(xy(1), xy(2), latit, longit)
    
    print *, xy
    print *, latit, longit
    
  end subroutine getncalt
    
    
    
    

    


end module nctools
