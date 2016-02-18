module ADE_fnc
  public :: ADEdispersion
  public :: reaction
  public :: ADE_std_convection
  public :: ADE_tder_coef
  public :: ADE_mass
  public :: ADE_reaction, ADE_zerorder
  
  contains
    subroutine ADEdispersion(pde_loc, layer, quadpnt, x, tensor, scalar)
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
      integer(kind=ikind) :: D, i
      
      
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
      
      call pde(1)%pde_fnc(1)%convection(pde(1), layer, quadpnt, vector_out = q_w(1:D))
      
      q_abs = 0.0
      do i=1, D
	q_abs =q_abs + q_w(i)*q_w(i)
      end do
      q_abs = sqrt(q_abs)
      tortuo = theta**(10.0/3.0)/(vgset(layer)%ths*vgset(layer)%ths)
      
      tensor = theta * (adepar(layer)%diff*q_abs + adepar(layer)%difmol*identity(1:D, 1:D))
      
      
    
    end subroutine ADEdispersion
    

    
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
    
    subroutine ADE_std_convection(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use global_objs
      use pde_objs
      use ADE_globals
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> input vector
      real(kind=rkind), dimension(:), intent(in), optional  :: vector_in
      !> output vector
      real(kind=rkind), dimension(:), intent(out), optional :: vector_out
      !> relative scalar value of the nonlinear function 
      real(kind=rkind), intent(out), optional               :: scalar
      
      if (present(vector_out)) then
	vector_out = adepar(layer)%convection
      end if
      
      
      if (present(scalar)) then
	scalar = abs(adepar(layer)%convection)
      end if
      
      
    end subroutine ADE_std_convection
    
    
    function ADE_tder_coef(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use ADE_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val
      
      real(kind=rkind) :: theta, n, ka, kd, csmax, cl
      
      if (drutes_config%name=="ADE_wr") then
	theta = pde(1)%mass(layer, quadpnt)
      else
	theta = adepar(layer)%water_cont
      end if
      
      if (.not. adepar(layer)%sorption%kinetic) then
	ka = adepar(layer)%sorption%adsorb
	kd = adepar(layer)%sorption%desorb
	if (ka > 10*epsilon(ka) .and. kd > 10*epsilon(kd)) then 
	  select case(adepar(layer)%sorption%name)
	    case("freund")
	      n = adepar(layer)%sorption%third
	      if (abs(n-1.0_rkind)>10*epsilon(n)) then
		cl = pde_loc%getval(quadpnt)
		val = theta+(1-theta)*ka/kd*adepar(layer)%bd*cl**(n-1)
	      else
		val = theta+(1-theta)*ka/kd*adepar(layer)%bd
	      end if
	    
	    case("langmu")
	      cl = pde_loc%getval(quadpnt)
	      csmax = adepar(layer)%sorption%third
	      val = theta + (1-theta)*adepar(layer)%bd*(ka*csmax)/(kd+ka*cl)
	    case default
	      print *, "unsupported sorption type, runtime error, called from ADE_fnc::ADE_tder_coef"
	      ERROR STOP
	    
	  end select
	else
	  val = theta
	end if
	
      else
	val = theta
      end if
      
      
      
      
    
    end function ADE_tder_coef
    
    
    function ADE_mass(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use ADE_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      
      val = pde_loc%getval(quadpnt)
      
    end function ADE_mass
    
    function ADE_reaction(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use ADE_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      
      integer(kind=ikind) :: n, i
      real(kind=rkind) :: theta, cl
      
     
      if (drutes_config%name=="ADE_wr") then
	theta = pde(1)%mass(layer, quadpnt)
      else
	theta = adepar(layer)%water_cont
      end if
      
      val = 0.0_rkind
      
      do i=1, ubound(adepar(layer)%orders,1)
	n = 10
	if (abs(adepar(layer)%orders(i) - 1.0_rkind) < 100*epsilon(1.0_rkind)) n = 1
	if (abs(adepar(layer)%orders(i)) < 100*epsilon(1.0_rkind)) n = 0
	select case(n)
	  case(0)
	    CONTINUE
	  case(1)
	    val = val + theta*adepar(layer)%lambda(i)
	  case default
	    cl = pde_loc%getval(quadpnt)
	    val = theta*adepar(layer)%lambda(i)*cl**(adepar(layer)%orders(i)-1)
	end select
      end do
	
	  
      
      
    end function ADE_reaction
    
    function ADE_zerorder(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use ADE_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      
      integer(kind=ikind) :: n, i
      real(kind=rkind) :: theta
      
      if (drutes_config%name=="ADE_wr") then
	theta = pde(1)%mass(layer, quadpnt)
      else
	theta = adepar(layer)%water_cont
      end if
      
       do i=1, ubound(adepar(layer)%orders,1)
	n = 10
	if (abs(adepar(layer)%orders(i)) < 100*epsilon(1.0_rkind)) n = 0
	select case(n)
	  case(0)
	    val =  val + theta*adepar(layer)%lambda(i)
	  case default
	    CONTINUE
	end select
      end do
      
      
    end function ADE_zerorder
    

end module ADE_fnc