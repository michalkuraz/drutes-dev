module ADE_pointers
  public :: ADE
  public :: ADEkinsorb
  
  contains
  
    subroutine ADE(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use ADE_fnc
      use ADE_reader
      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i
      
      call ADE_read(pde_loc)
	    
      pde_loc%pde_fnc(pde_loc%order)%dispersion => ADEdispersion
      
      pde_loc%pde_fnc(pde_loc%order)%convection => ADE_convection

      pde_loc%pde_fnc(pde_loc%order)%elasticity => ADE_tder_coef

      pde_loc%mass => ADE_mass

      pde_loc%pde_fnc(pde_loc%order)%reaction => ADE_reaction
            
      pde_loc%pde_fnc(pde_loc%order)%zerord => ADE_zerorder
	  
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
	select case(pde_loc%bc(i)%code)
	  case(1)
	    pde_loc%bc(i)%value_fnc => ADE_dirichlet
	  case(2)
	    pde_loc%bc(i)%value_fnc => ADE_neumann
	end select
      end do    
	
      pde_loc%flux => ADE_flux
      
      pde_loc%initcond => ADE_icond  
      
    
    end subroutine ADE
    
    subroutine ADEkinsorb(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use ADE_fnc
      use ADE_reader
      
      class(pde_str), intent(in out) :: pde_loc  
      integer(kind=ikind) :: i
      
      call ADEcs_read(pde_loc)
      
      pde_loc%pde_fnc(pde_loc%order)%elasticity => ADE_tder_cscs
      
      pde(pde_loc%order-1)%pde_fnc(pde_loc%order)%elasticity => ADE_tder_cscl
      
      pde_loc%pde_fnc(pde_loc%order-1)%reaction => ADE_cscl_react
      
      pde_loc%pde_fnc(pde_loc%order)%reaction => ADE_cscs_react
      
      allocate(pde_loc%bc(lbound(pde(pde_loc%order-1)%bc,1) : (ubound(pde(pde_loc%order-1)%bc,1) )  ))
      
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
	pde_loc%bc(i)%code = 2
	pde_loc%bc(i)%value_fnc => ADE_null_bc
      end do 
      
      pde_loc%initcond => ADEcs_icond
      
    
    end subroutine ADEkinsorb
    
      
  

end module ADE_pointers