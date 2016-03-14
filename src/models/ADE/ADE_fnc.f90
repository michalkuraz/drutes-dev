module ADE_fnc
  public :: ADEdispersion
  public :: ADE_convection
  public :: ADE_tder_coef, ADE_tder_cscl, ADE_tder_cscs
  public :: ADE_mass
  public :: ADE_reaction, ADE_zerorder, ADE_flux, ADE_icond, ADE_csbc
  public :: ADE_cscl_react, ADE_cscs_react
  public :: ADEcs_icond, ADEcs_mass
  public :: ADE_dirichlet, ADE_neumann, ADE_null_bc
  
  contains
    subroutine ADEdispersion(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use globals
      use ADE_globals
      use re_globals
      use debug_tools
      
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
      real(kind=rkind) :: theta, q_abs, tortuo, ths
      integer(kind=ikind) :: D, i
      
     
      D = drutes_config%dimen
      identity = 0.0
      do i=1, D
	identity(i,i) = 1.0
      end do
      
      select case(drutes_config%name)
	case("ADEstd", "ADEstd_kinsorb")
	  q_w = ADEpar(layer)%convection
	  theta = adepar(layer)%water_cont
	  ths = adepar(layer)%water_cont
	case default
	  theta = pde(1)%mass(layer, quadpnt)
	  call pde(1)%flux(layer, quadpnt, vector_out = q_w(1:D))
	  ths = vgset(layer)%ths
      end select
      

      q_abs = 0.0
      do i=1, D
	q_abs =q_abs + q_w(i)*q_w(i)
      end do
      q_abs = sqrt(q_abs)
      tortuo = theta**(10.0/3.0)/(ths*ths)
      
      
      if (present(tensor)) then
	tensor = theta * (adepar(layer)%diff*q_abs + adepar(layer)%difmol*identity(1:D, 1:D))	
      end if
      

      
      
    
    end subroutine ADEdispersion
    
    
    subroutine ADE_convection(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
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
      
      
      if (pde_loc%order == 1) then
	if (present(vector_out)) then
	  vector_out = adepar(layer)%convection
	end if
	
	
	if (present(scalar)) then
	  scalar = abs(adepar(layer)%convection)
	end if
      else
	if (present(vector_out)) then
	  call pde(1)%flux(layer, quadpnt, vector_out=vector_out)
	end if
	
        if (present(scalar)) then
	  call pde(1)%flux(layer, quadpnt, scalar=scalar)
	end if
	
      end if
      
      
    end subroutine ADE_convection
    
    
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
      
      
      if (pde_loc%order == 2) then
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
    
    
    function ADE_tder_cscl(pde_loc, layer, quadpnt, x) result(val)
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
      
      real(kind=rkind) :: theta
      
      if (pde_loc%order == 2) then
	theta = pde(1)%mass(layer, quadpnt)
      else
	theta = adepar(layer)%water_cont
      end if
      
      val = 1.0_rkind-theta
      
    end function ADE_tder_cscl
    
    
    function ADE_tder_cscs(pde_loc, layer, quadpnt, x) result(val)
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
      
      val = 1.0_rkind
    
    end function ADE_tder_cscs
    
    
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
      real(kind=rkind)                :: theta
      
      if (pde_loc%order == 1) then
	theta = adepar(layer)%water_cont
      else
	theta = pde(1)%mass(layer, quadpnt)
      end if
      
      if (present(quadpnt)) then
	val = theta*pde_loc%getval(quadpnt)
      else
        val = theta * x(1)
      end if
      
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
      
     
      if (pde_loc%order == 2) then
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
      
      if (pde_loc%order == 2) then
	theta = pde(1)%mass(layer, quadpnt)
      else
	theta = adepar(layer)%water_cont
      end if
      
       val = 0.0_rkind
       do i=1, ubound(adepar(layer)%orders,1)
	if (abs(adepar(layer)%orders(i)) < 100*epsilon(1.0_rkind)) then
	  val =  val + theta*adepar(layer)%lambda(i)
	end if
      end do
      
      
    end function ADE_zerorder
    
    subroutine ADE_dirichlet(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      
      integer(kind=ikind) :: edge_id, i, j, proc
      real(kind=rkind) :: tempval
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      
      if (present(value)) then
	if (pde_loc%bc(edge_id)%file) then
	  do i=1, ubound(pde_loc%bc(edge_id)%series,1)
	    if (pde_loc%bc(edge_id)%series(i,1) > time) then
	      if (i > 1) then
		j = i-1
	      else
		j = i
	      end if
	      tempval = pde_loc%bc(edge_id)%series(j,2)
	      EXIT
	    end if
	  end do
	else
	  tempval =  pde_loc%bc(edge_id)%value
	end if
	value = tempval 
      end if


      
      if (present(code)) then
	code = 1
      end if
      

    end subroutine ADE_dirichlet
    
    
    subroutine ADE_neumann(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      
      integer(kind=ikind) :: edge_id, i, j, proc
      real(kind=rkind) :: tempval
      
      edge_id = nodes%edge(elements%data(el_id, node_order))
      
      if (present(value)) then
	if (pde_loc%bc(edge_id)%file) then
	  do i=1, ubound(pde_loc%bc(edge_id)%series,1)
	    if (pde_loc%bc(edge_id)%series(i,1) > time) then
	      if (i > 1) then
		j = i-1
	      else
		j = i
	      end if
	      tempval = pde_loc%bc(edge_id)%series(j,2)
	      EXIT
	    end if
	  end do
	else
	  tempval =  pde_loc%bc(edge_id)%value
	end if
	value = tempval 
      end if


      
      if (present(code)) then
	code = 2
      end if
      

    end subroutine ADE_neumann
    
    
     subroutine ADE_null_bc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code

      if (present(value)) then
	value = 0.0_rkind
      end if

      if (present(code)) then
	code = 2
      end if
	
    end subroutine ADE_null_bc
    
    
    subroutine ADE_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use debug_tools
      use ADE_globals
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
    
      real(kind=rkind), dimension(:,:), allocatable, save  :: Dhm
      real(kind=rkind), dimension(:), allocatable, save :: q_w, gradC
      real(kind=rkind) :: c, cmax
      
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
	print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
	print *, "exited from ADE_fnc::ADE_flux"
	ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from ADE_fnc::ADE_flux"
	ERROR stop
      end if
      
      if (.not. allocated(q_w)) allocate(q_w(drutes_config%dimen))

      if (.not. allocated(gradC)) allocate(gradC(drutes_config%dimen))
      
      if (.not. allocated(Dhm)) allocate(Dhm(drutes_config%dimen, drutes_config%dimen))
      
      if (present(quadpnt)) then
	c = pde_loc%getval(quadpnt)
	call pde_loc%getgrad(quadpnt, gradC)
      else
        if (ubound(x,1) /=1) then
	  print *, "ERROR: van Genuchten function is a function of a single variable h"
	  print *, "       your input data has:", ubound(x,1), "variables"
	  print *, "exited from ADE_fnc::ADE_flux"
	  ERROR STOP
	end if
	c = x(1)
	gradC = grad
      end if
      
      
      call pde_loc%pde_fnc(1)%dispersion(pde_loc, layer, quadpnt, tensor=Dhm)
      
      select case(pde_loc%order)
	case(1)
	  q_w = adepar(layer)%convection
	case(2)
	  call pde(1)%flux(layer, quadpnt, vector_out=q_w)
      end select
      
      select case(adepar(layer)%icondtype)
	case("ca")
	  cmax = 1.0_rkind
	 case("cr")
	  cmax = adepar(layer)%cmax
      end select
      
      
      if (present(flux)) then
	flux = cmax*matmul(Dhm, gradC) + cmax*q_w*c
      end if
    
    end subroutine ADE_flux
    
    subroutine ADE_icond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use ADE_globals

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
   
      D = drutes_config%dimen
      do i=1, elements%kolik
	layer = elements%material(i,1)
	do j=1, ubound(elements%data,2)
	  k = elements%data(i,j)
	  l = nodes%edge(k)
	  m = pde_loc%permut(k)
	  if (m == 0) then
	    call pde_loc%bc(l)%value_fnc(pde_loc, i, j, value)
	    pde_loc%solution(k) =  value 
	  else
	    select case (adepar(layer)%icondtype)
	      case("ca")
		pde_loc%solution(k) = adepar(layer)%cinit
	      case("cr")
		pde_loc%solution(k) = adepar(layer)%cinit * adepar(layer)%cmax
	    end select
	  end if
	end do   
      end do

    
    
    end subroutine ADE_icond
    
    function ADE_cscs_react(pde_loc, layer, quadpnt, x) result(val)
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
      
      integer(kind=ikind) :: proc_cl
      real(kind=rkind) :: cs, cl
      
      proc_cl = pde_loc%order - 1
      
      
      select case(adepar(layer)%sorption%name)
	case("langmu")
	  cs = pde_loc%getval(quadpnt) 
          cl = pde(proc_cl)%getval(quadpnt)
          val = -adepar(layer)%sorption%adsorb*cl - adepar(layer)%sorption%desorb
	case("freund")
	  val = -adepar(layer)%sorption%desorb
      end select
      
      
  end function ADE_cscs_react
  
   function ADE_cscl_react(pde_loc, layer, quadpnt, x) result(val)
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
      
      integer(kind=ikind) :: proc_cl
      real(kind=rkind) :: cl 
      

      
      select case(adepar(layer)%sorption%name)
	case("langmu")
	  val = adepar(layer)%sorption%adsorb*adepar(layer)%sorption%third
	case("freund")
	  if (abs(adepar(layer)%sorption%third - 1.0_rkind ) > 10*epsilon(1.0_rkind)) then
	    proc_cl = pde_loc%order - 1
	    cl = pde(proc_cl)%getval(quadpnt)
	    val = adepar(layer)%sorption%adsorb*cl**(1-adepar(layer)%sorption%third)
	  else
	    val = adepar(layer)%sorption%adsorb
	  end if
      end select
      
   end function ADE_cscl_react
  
  
    subroutine ADE_csbc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use debug_tools
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      
      integer(kind=ikind) :: edge_id, i, j, proc
      real(kind=rkind) :: tempval

     
      
      if (present(value)) then
	value = 0.0
      end if


      
      if (present(code)) then
	code = 0
      end if
      

    end subroutine ADE_csbc
    
    subroutine ADEcs_icond(pde_loc) 
      use typy
      use pde_objs

      
      class(pde_str), intent(in out) :: pde_loc
      
      pde_loc%solution = 0.0_rkind
    
    end subroutine ADEcs_icond
    
    
    function ADEcs_mass(pde_loc, layer, quadpnt, x) result(val)
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
      real(kind=rkind)                :: theta
      
      if (pde_loc%order == 2) then
	theta = adepar(layer)%water_cont
      else 
	theta = pde(1)%mass(layer, quadpnt)
      end if
      
      if (present(quadpnt)) then
	val = pde_loc%getval(quadpnt)*(1-theta)
      else
	print *, "exited from ADE_fnc::ADEcs_mass"
	ERROR STOP
      end if
    
    end function ADEcs_mass
    

end module ADE_fnc