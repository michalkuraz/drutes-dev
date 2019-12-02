module vapour_pointers
  public :: vapour_processes
  public :: vapour_linker
  
  contains
    subroutine vapour_processes(processes)
      use typy
      
      integer(kind=ikind), intent(out) :: processes
      
      processes = 2
      
    end subroutine vapour_processes
    
    
    subroutine vapour_linker()
      use typy
      use debug_tools
      use RE_pointers
      use pde_objs
      use evap_reader
      use evap_fnc
      use evap_globals
      use heat_pointers
      use heat_fnc
      use re_constitutive
      
!       call RE_std(pde(re_order))
      call allREpointers(pde(re_order))
      
      call set_evapbc(pde(re_order))
      
      pde(re_order)%flux => liquid_flux
      
      pde(re_order)%initcond => re_initcond
      
      call evap_var()
      
      !> Richards modified equation
      
      pde(re_order)%pde_fnc(re_order)%dispersion => difussion_hh
      
      pde(re_order)%pde_fnc(heat_order)%dispersion => difussion_hT
      
      pde(re_order)%pde_fnc(re_order)%zerord  =>  source_h
      
      call heatlinker(pde(heat_order))
      
      !> Heat modified equation
      
      pde(heat_order)%pde_fnc(heat_order)%elasticity => capacity_T
      
      pde(heat_order)%pde_fnc(heat_order)%dispersion => difussion_TT
      
      pde(heat_order)%pde_fnc(re_order)%dispersion => difussion_Th
      
      pde(heat_order)%pde_fnc(heat_order)%convection => convection_T
      
      pde(heat_order)%pde_fnc(heat_order)%zerord  =>  source_T
      
      pde(heat_order)%initcond => heat_icond  
      
      
      pde(heat_order)%flux => heatmod_flux
    
    end subroutine vapour_linker
    
    
    
    subroutine set_evapbc(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use evap_bc
      
      class(pde_str), intent(in out) :: pde_loc
      
      integer(kind=ikind) :: i
      
      
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
        select case(pde_loc%bc(i)%code)
          case(-1)
            pde_loc%bc(i)%value_fnc => re_dirichlet_height_bc
          case(0)
            pde_loc%bc(i)%value_fnc => re_null_bc
          case(1)
            pde_loc%bc(i)%value_fnc => re_dirichlet_bc
          case(2)
            pde_loc%bc(i)%value_fnc => re_neumann_bc
          case(3)
            pde_loc%bc(i)%value_fnc => re_null_bc
          case(4)
            print *, "seepage face is not created yet for pressure head RE"
            ERROR STOP
          case(5)
            pde_loc%bc(i)%value_fnc => water_evap
          case default
            print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
            print *, "the incorrect boundary code specified is:", pde_loc%bc(i)%code
            ERROR stop
        end select
      end do
	

    end subroutine set_evapbc
    

end module vapour_pointers
