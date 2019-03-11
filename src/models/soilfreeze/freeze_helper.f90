module freeze_helper
  use pde_objs
  use typy
  use freeze_globs
  use debug_tools
  use RE_constitutive

  public :: iceswitch, rho_icewat, Q_reduction, surf_tens_deriv, Kliquid_temp, hl
      
      
  
  contains

      !> switch for freezing condition based on Clapeyron equation 
    function iceswitch(quadpnt) result(sw)
      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      logical :: sw
      
      real(kind=rkind) :: Tf
      type(integpnt_str) :: quadpnt_loc
      
      quadpnt_loc = quadpnt
      quadpnt_loc%column = 1
      
      Tf = Tref*exp(pde(1)%getval(quadpnt_loc)*grav/Lf)
      Tf = Tf - 273.15_rkind

      if (pde(2)%getval(quadpnt_loc) > Tf) then
      !> melting
        sw = .FALSE.
      else
      !> freezing
        sw = .TRUE.
      end if
          
    end function iceswitch

    function rho_icewat(quadpnt) result(rho)

      use typy
      use global_objs
      
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: rho
      
      integer(kind=ikind) :: layer, el
      real(kind=rkind) :: thl, thall, thice
      
      
      if (quadpnt%type_pnt == "ndpt" ) then
        el = nodes%element(quadpnt%order)%data(1)
      else
        el = quadpnt%element
      end if
      
      layer = elements%material(el)
      thl = vangen(pde(1), layer, x=(/hl(quadpnt)/))
      thall = vangen(pde(1), layer, quadpnt)
      thice = thall - thl
      rho = (thl * rho_wat + thice * rho_ice)/thall
       
    end function rho_icewat
    
    function Q_reduction(layer, quadpnt, x) result(val)

      use typy
      use global_objs
      use re_globals
      
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      integer(kind=ikind) :: el
      real(kind=rkind) :: thl, thall, thice, val
      if(present(quadpnt)) then
        thall = vangen(pde(1), layer, quadpnt)
        thl = vangen(pde(1), layer, x=(/hl(quadpnt)/))
     ! else if (present(x)) then
      !  thall = vangen(pde(1), layer, x = x)
      !  thl = vangen(pde(1), layer, x = x)
      end if
      
      thice = thall - thl
      val = thice/(thall- vgset(layer)%Thr)
       
    end function Q_reduction
    
    subroutine Kliquid_temp(pde_loc, layer, quadpnt, x, T, tensor, scalar) 
      use typy
      use global_objs
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind),intent(in), optional    :: T
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      real(kind=rkind), dimension(3,3) :: Klh, Klt
      integer(kind=ikind) :: D

      D = drutes_config%dimen

      if (present(tensor)) then
        call mualem(pde_loc, layer, quadpnt, tensor = Klt(1:D, 1:D))

        if (present(quadpnt)) then
          tensor = Klt(1:D, 1:D)*gwt*pde(1)%getval(quadpnt)*surf_tens_deriv(pde_loc, layer, quadpnt)
        !else if (present(T)) then
        !  if (present (x)) then
        !    tensor = Klt(1:D, 1:D)*gwt*x(1)*surf_tens_deriv(pde_loc, layer, T = T)
        !  end if
        else
          print *, "runtime error"
          print *, "exited from Kliquid_temp::freeze_fnc"
          ERROR STOP
        end if
      else
         print *, "ERROR! output tensor undefined, exited from Kliquid_temp::freeze_fnc"
      end if    

      
    end subroutine Kliquid_temp
    
    
    function surf_tens_deriv(pde_loc, layer, quadpnt, T) result(val)
      use typy
      use global_objs
      use freeze_globs
      use pde_objs
    
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), intent(in), optional    ::  T
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
    
      real(kind=rkind) :: temp
    
      temp = pde(2)%getval(quadpnt)
      if (present(T)) then
        temp = T
      end if
      
      val = -0.1425-2.38e-4*temp
      
    end function surf_tens_deriv
    
    function hl(quadpnt) result(val)
      use typy
      use global_objs
      use freeze_globs
     
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: val, T_melt
      
      real(kind=rkind) :: hw, temp
      
      hw = pde(1)%getval(quadpnt)
      
      temp = pde(2)%getval(quadpnt)
      T_melt = Tref*exp(hw*grav/Lf)
      
      if(iceswitch(quadpnt)) then
        val = hw+Lf/grav*log((temp+273.15_rkind)/T_melt) !units
      else
        val = hw
      end if
      
    end function hl
    
    function thetai(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: layer
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      real(kind=rkind) :: val
      
      real(kind=rkind) :: thl, thall
      
      thl = vangen(pde(1), layer, x=(/hl(quadpnt)/))
      thall = vangen(pde(1), layer, quadpnt)
      
      val = (thall * rho_icewat(quadpnt) - thl * rho_wat)/rho_ice
    end function thetai
end module freeze_helper
