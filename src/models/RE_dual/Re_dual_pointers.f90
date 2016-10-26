module Re_dual_pointers
  public :: RE_fracture
  public :: RE_matrix
  
  contains

    subroutine RE_matrix(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_por
      use Re_dual_reader
      use RE_constitutive
      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i
      
      call Re_dual_readm(pde_loc)
	  call Re_dual_var() 

      pde_loc%pde_fnc(pde_loc%order)%dispersion => dual_mualemm
      
      pde_loc%pde_fnc(pde_loc%order)%reaction => dual_coupling

      pde_loc%pde_fnc(pde_loc%order)%elasticity => dual_ret_capm
      
      pde_loc%mass => vangen_d_m
      
      pde_loc%flux => darcy_law_d!implement darcy law for fluxes
      
      ! boundary condition defined as different type boundary_vals
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
	   select case(pde_loc%bc(i)%code)
	     case(1)
	       pde_loc%bc(i)%value_fnc => re_dirichlet_bc
	     case(2)
	       pde_loc%bc(i)%value_fnc => re_neumann_bc
	   end select
      end do  
    pde_loc%initcond => dual_inicond_m
   
    end subroutine RE_matrix
    
    subroutine RE_fracture(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_por
      use Re_dual_reader
      use RE_constitutive
      
      class(pde_str), intent(in out) :: pde_loc  
      integer(kind=ikind) :: i
      
      call Re_dual_readf(pde_loc)
      
      pde_loc%pde_fnc(pde_loc%order)%dispersion => dual_mualemf
      
      pde_loc%pde_fnc(pde_loc%order)%reaction => dual_coupling_f

      pde_loc%pde_fnc(pde_loc%order)%elasticity => dual_ret_capf
      
      pde_loc%mass => vangen_d_f
      
      pde_loc%flux => darcy_law_d
      
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
	   select case(pde_loc%bc(i)%code)
	     case(1)
	       pde_loc%bc(i)%value_fnc => re_dirichlet_bc
	     case(2)
	       pde_loc%bc(i)%value_fnc => re_neumann_bc
	   end select
     end do  
     
    pde_loc%initcond => dual_inicond_f
    end subroutine RE_fracture
      

end module Re_dual_pointers