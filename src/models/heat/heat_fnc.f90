module heat_fnc
  public :: heat_conduct
  public :: heat_convect
  public :: heat_source
  public :: heat_dirichlet, heat_neumann
  public :: heat_flux
  public :: heat_icond
  public :: heat_elast
  
  contains
    subroutine heat_conduct(pde_loc, layer, quadpnt, x, tensor, scalar)
      use typy
      use global_objs
      use pde_objs
      use globals
      use heat_globals
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

      integer(kind=ikind) :: D, i
      
     
      D = drutes_config%dimen

      
      if (present(tensor)) then
	tensor =  heatpar(layer)%lambda
      end if
      
    
    end subroutine heat_conduct
    
    
    subroutine heat_convect(pde_loc, layer, quadpnt, x, vector_in, vector_out, scalar)
      use typy
      use global_objs
      use pde_objs
      use heat_globals
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
	  vector_out = heatpar(layer)%C_w*heatpar(layer)%convection
	end if
	
	
	if (present(scalar)) then
	  scalar = norm2(heatpar(layer)%C_w*heatpar(layer)%convection)
	end if
   
      
    end subroutine heat_convect
    
    function heat_elast(pde_loc,layer, quadpnt, x) result(E)
      use typy
      use heat_globals
      use pde_objs
      use core_tools

      class(pde_str), intent(in) :: pde_loc 
      integer(kind=ikind), intent(in) :: layer
      !> pressure head
      real(kind=rkind), intent(in), dimension(:),  optional :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      real(kind=rkind) :: h
      !> resulting system elasticity
      real(kind=rkind) :: E

       
      
      if (present(quadpnt) .and. present(x)) then
        print *, "ERROR: the function can be called either with integ point or x value definition, not both of them"
        print *, "exited from re_constitutive::vangen_elast"
        ERROR stop
      else if (.not. present(quadpnt) .and. .not. present(x)) then
        print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from re_constitutive::vangen_elast"
        ERROR stop
      end if


      E = heatpar(layer)%C
      

    end function heat_elast 

    function heat_source(pde_loc, layer, quadpnt, x) result(val)
      use typy
      use global_objs
      use pde_objs
      use heat_globals
      
      class(pde_str), intent(in) :: pde_loc
      !> value of the nonlinear function
      real(kind=rkind), dimension(:), intent(in), optional    :: x
      !> Gauss quadrature point structure (element number and rank of Gauss quadrature point)
      type(integpnt_str), intent(in), optional :: quadpnt
      !> material ID
      integer(kind=ikind), intent(in) :: layer
      !> return value
      real(kind=rkind)                :: val 
      
      
      val = heatpar(layer)%source
      
      
    end function heat_source
    
    subroutine heat_dirichlet(pde_loc, el_id, node_order, value, code) 
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
      

    end subroutine heat_dirichlet
    
    
    subroutine heat_neumann(pde_loc, el_id, node_order, value, code) 
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
      

    end subroutine heat_neumann
    
    
    
    
    subroutine heat_flux(pde_loc, layer, quadpnt, x, grad,  flux, flux_length)
      use typy
      use pde_objs
      use global_objs
      use debug_tools
      use heat_globals
       
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)                          :: layer
      type(integpnt_str), intent(in), optional :: quadpnt    
      real(kind=rkind), intent(in), dimension(:), optional                   :: x
      !> this value is optional, because it is required by the vector_fnc procedure pointer global definition
      real(kind=rkind), dimension(:), intent(in), optional     :: grad
      real(kind=rkind), dimension(:), intent(out), optional    :: flux
      real(kind=rkind), intent(out), optional                  :: flux_length
    

      real(kind=rkind), dimension(:), allocatable, save :: gradT

      
      
      if (present(quadpnt) .and. (present(grad) .or. present(x))) then
	print *, "ERROR: the function can be called either with integ point or x value definition and gradient, not both of them"
	print *, "exited from heat_fnc::heat_flux"
	ERROR stop
      else if ((.not. present(grad) .or. .not. present(x)) .and. .not. present(quadpnt)) then
	print *, "ERROR: you have not specified either integ point or x value"
        print *, "exited from heat_fnc::heat_flux"
	ERROR stop
      end if   

      if (.not. allocated(gradT)) allocate(gradT(drutes_config%dimen))

      if (present(quadpnt)) then
	call pde_loc%getgrad(quadpnt, gradT)
      else
	gradT = grad
      end if
      
      
      if (present(flux)) then
	flux = matmul(heatpar(layer)%lambda, gradT) 
      end if
      
      if (present(flux_length)) then
        flux_length = norm2(matmul(heatpar(layer)%lambda, gradT))
      end if
    
    end subroutine heat_flux
    
    subroutine heat_icond(pde_loc) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use heat_globals

      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i, j, k,l, m, layer, D
      real(kind=rkind) :: value
      
   
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
	    pde_loc%solution(k) = heatpar(layer)%Tinit
	  end if
	end do   
      end do

    
    
    end subroutine heat_icond
    
    

end module heat_fnc