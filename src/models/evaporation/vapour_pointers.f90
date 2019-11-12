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
      
      call RE_std(pde(re_order))

      call evap_var()
      
      pde(re_order)%pde_fnc(re_order)%dispersion => difussion_hh
      
      pde(re_order)%pde_fnc(heat_order)%dispersion => difussion_hT
      
      !> function dtheta_vapordt doesn't contain root water uptake term, must be reformulated
      pde(re_order)%pde_fnc(re_order)%zerord  =>  dtheta_vapordt
      
      call heatlinker(pde(heat_order))
      
      pde(heat_order)%pde_fnc(heat_order)%elasticity => capacity_T
      
      pde(heat_order)%pde_fnc(heat_order)%dispersion => difussion_TT
      
      pde(heat_order)%pde_fnc(re_order)%dispersion => difussion_Th
      
      pde(heat_order)%pde_fnc(heat_order)%convection => convection_T
    
    end subroutine vapour_linker

end module vapour_pointers
