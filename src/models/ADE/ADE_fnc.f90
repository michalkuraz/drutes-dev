module ADE_fnc
  public :: dispersion
  public :: convection
  public :: reaction
  
  contains
    subroutine dispersion(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use globals
      use ADE_globals
      use re_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return tensor
      real(kind=rkind), dimension(:,:), intent(out), optional :: tensor
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional                 :: scalar
      
      real(kind=rkind), dimension(3,3) :: identity
      real(kind=rkind), dimension(3) :: q_w
      real(kind=rkind) :: theta, q_abs, tortuo
      integer(kind=ikind) :: D
      
      
      D = drutes_config%dimen
      identity = 0.0
      do i=1, D
	identity(i,i) = 1.0
      end do
      
      if (drutes_config%name=="ADEstd") then
	tensor = (ADEpar(layer)%diff*abs(ADEpar(layer)%convection) + &
	  identity(1:D, 1:D)*adepar(layer)%difmol)* &
	  adepar(layer)%water_cont
	  RETURN
      end if
      
      theta = pde(1)%mass(layer, quadpnt)
      
      call pde(1)%pde_fnc(1)%convection(layer, quadpnt, vector_out = q_w(1:D))
      
      q_abs = 0.0
      do i=1, D
	q_abs =q_abs + q_w(i)*q_w(i)
      end do
      q_abs = sqrt(q_abs)
      tortuo = theta**(10.0/3.0)/(vgset(layer)%ths*vgset(layer)%ths)
      
      tensor = theta * (adepar(layer)%diff*q_abs + adepar(layer)%diffmol*identity(1:D, 1:D))
      
      
    
    end subroutine dispersion
    

    
    function reaction(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
    
    end function reaction
    
    

end module ADE_fnc