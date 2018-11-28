module freeze_pointers
  public :: freeze_processes, frz_pointers
  
  contains
  
    subroutine freeze_processes(processes)
      use typy
      use globals
      
      integer(kind=ikind), intent(out) :: processes
      
      
      select case (drutes_config%name)
        case ("freeze")
          processes = 2
          
        case ("LTNE")
          processes = 4
          
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::freeze_processes"
          error stop
      end select
        
    
    end subroutine freeze_processes
  
  
  
    subroutine frz_pointers()
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use freeze_globs
      use freeze_fnc
      use heat_pointers
      
      
      call heat(pde(1:2))
      
      pde(1)%pde_fnc(1)%elasticity => capacityhh
      
      pde(1)%pde_fnc(2)%elasticity => capacityhT
      
      pde(1)%pde_fnc(1)%dispersion => diffhh
      
      pde(1)%pde_fnc(2)%dispersion => diffhT
      
      pde(2)%pde_fnc(1)%elasticity => capacityTh
      
      pde(2)%pde_fnc(2)%elasticity => capacityTT
      
      pde(2)%pde_fnc(2)%dispersion => diffTT
      
      pde(2)%pde_fnc(2)%convection => convectTT
      
      if (drutes_config%fnc_method == 0) then
        Kliquid => mualem
        rwcap => vangen_elast
        theta => vangen
      else
        Kliquid  => mualem_tab		
        rwcap => vangen_elast_tab
        theta => vangen_tab
      end if
      
      pde(1)%problem_name(1) = "RE_freeze_thaw"
      pde(1)%problem_name(2) = "Richards' equation with freezing andd thawing processes"

      pde(1)%solution_name(1) = "press_head" !nazev vystupnich souboru
      pde(1)%solution_name(2) = "h  [L]" !popisek grafu

      pde(1)%flux_name(1) = "flux"  
      pde(1)%flux_name(2) = "Darcian flow [L.T^{-1}]"

      pde(1)%mass_name(1) = "theta"
      pde(1)%mass_name(2) = "theta [-]"
      
      
    
    end subroutine frz_pointers
  
  



end module freeze_pointers
