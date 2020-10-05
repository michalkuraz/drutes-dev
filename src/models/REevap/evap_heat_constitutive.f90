module evap_heat_constitutive

  public :: heatcap_TT, heatcap_hT

  private :: water_cap, vapour_cap, latentheat
  
  contains
  
    function heatcap_TT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      use re_globals
      use heat_globals
      use re_constitutive
      use evap_RE_constitutive
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: ths, theta, theta_v
      
      if (.not. present(quadpnt)) then
        print *, "runtime error evap_heat_constitutive::heatcap_TT"
        ERROR STOP
      end if
      
      ths = vgset(layer)%ths
      theta = vangen(pde(re_ord), layer, quadpnt)
      theta_v = thetav(pde(re_ord), layer, quadpnt) 
      
      val = heatpar(layer)%C*(1-ths) + water_cap(quadpnt)*theta + vapour_cap(quadpnt)*theta_v + &
            latentheat(quadpnt)*dthetav_dtemp(pde(re_ord), layer, quadpnt)
      
    end function heatcap_TT
    
    function heatcap_hT(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use evapglob
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      if (.not. present(quadpnt)) then
        print *, "runtime error evap_heat_constitutive::heatcap_hT"
        ERROR STOP
      end if
      
    end function heatcap_hT
    
    function water_cap(quadpnt) result(val)
      use typy
      use global_objs
      use evapglob
      use evap_RE_constitutive
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      val = pde(heat_ord)%getval(quadpnt) * dens_liquid(quadpnt)
    
    end function water_cap
    
    
    function vapour_cap(quadpnt) result(val)
      use typy
      use global_objs
      use evapglob
      use evap_RE_constitutive
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      val = C_vap * relhumid(quadpnt) * dens_satvap(quadpnt)
    
    end function vapour_cap
    
    
    function latentheat(quadpnt) result(val)
      use typy
      use global_objs
      use evapglob
      use pde_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: val
      
      val = 2.501e6 - 2369.2*pde(heat_ord)%getval(quadpnt)
    
    end function latentheat


end module evap_heat_constitutive
