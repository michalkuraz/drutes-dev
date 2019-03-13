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
      use re_total
      use re_reader
      use re_pointers
      use freeze_globs
      use freeze_fnc
      use heat_pointers
      use heat_reader
      use heat_fnc
      
      integer(kind=ikind) :: i
      
      call heat_read(pde(2))
      call res_read(pde(1))
      !call RE_totheadbc(pde(1))
      call RE_pressheadbc(pde(1))
    ! pointers for water flow model
      pde(1)%pde_fnc(1)%elasticity => capacityhh
      
      pde(1)%pde_fnc(2)%elasticity => capacityhT
      
      pde(1)%pde_fnc(1)%dispersion => diffhh
      
      !pde(1)%pde_fnc(1)%convection => convz

      pde(1)%pde_fnc(2)%dispersion => diffhT
      
      pde(1)%flux => all_fluxes
      !pde(1)%flux => darcy4totH

      pde(1)%initcond => retot_initcond
      
      if (drutes_config%fnc_method == 0) then
        rwcap => vangen_elast
      else
        rwcap => vangen_elast_tab
      end if
      
      pde(1)%problem_name(1) = "RE_freeze_thaw"
      pde(1)%problem_name(2) = "Richards' equation with freezing and thawing processes"

      pde(1)%solution_name(1) = "press_head" 
      pde(1)%solution_name(2) = "h [L]" 

      pde(1)%flux_name(1) = "flux"  
      pde(1)%flux_name(2) = "Darcian flow [L.T^{-1}]"

      deallocate(pde(1)%mass_name)
      deallocate(pde(1)%mass)
      
      allocate(pde(1)%mass_name(2,2))
      allocate(pde(1)%mass(2))

      pde(1)%mass_name(1,1) = "theta_l"
      pde(1)%mass_name(1,2) = "theta_l [-]"
      
      pde(1)%mass_name(2,1) = "theta_i"
      pde(1)%mass_name(2,2) = "theta_i [-]"
            
      
      pde(1)%mass(1)%val => vangen
      
      pde(1)%mass(2)%val => thetai
      ! allocate porosity as mass(3)?

    ! pointers for heat flow model
      pde(2)%pde_fnc(1)%elasticity => capacityTh
      
      pde(2)%pde_fnc(2)%elasticity => capacityTT
      
      pde(2)%pde_fnc(2)%dispersion => diffTT
      
      pde(2)%pde_fnc(2)%convection => convectTT
      
      pde(2)%flux => heat_flux
      
      pde(2)%initcond => heat_icond  
      
      do i=lbound(pde(2)%bc,1), ubound(pde(2)%bc,1)
        select case(pde(2)%bc(i)%code)
          case(1)
            pde(2)%bc(i)%value_fnc => heat_dirichlet
          case(2)
            pde(2)%bc(i)%value_fnc => heat_neumann
          case(0)
            pde(2)%bc(i)%value_fnc => re_null_bc
        end select
      end do  
    
    end subroutine frz_pointers
  
  



end module freeze_pointers
