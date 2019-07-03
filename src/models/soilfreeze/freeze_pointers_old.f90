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
      use freeze_helper
      use freeze_read
      use heat_pointers
      use heat_reader
      use heat_fnc
      
      integer(kind=ikind) :: i
      
      select case (drutes_config%name)
        case ("freeze")
          call freeze_reader(pde(wat))
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::frz_pointers"
          error stop
      end select
   
      
      call RE_totheadbc(pde(wat))
      pde(wat)%getval => getval_retotfr
      
    ! pointers for water flow model
      pde(wat)%pde_fnc(wat)%elasticity => capacityhh
      
      pde(wat)%pde_fnc(heat_proc)%elasticity => capacityhT
      
      pde(wat)%pde_fnc(wat)%dispersion => diffhh
      
      pde(wat)%pde_fnc(heat_proc)%dispersion => diffhT
      
      pde(wat)%flux => all_fluxes
      do i=lbound(pde(wat)%bc,1), ubound(pde(wat)%bc,1)
        select case(pde(wat)%bc(i)%code)
          case(6)
            pde(wat)%bc(i)%value_fnc => Dirichlet_mass_bc
          case(7)
            pde(wat)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
        end select
      end do
      pde(wat)%initcond => wat_initcond
      
      rwcap => vangen_elast_fr
      
      pde(wat)%problem_name(1) = "RE_freeze_thaw"
      pde(wat)%problem_name(2) = "Richards' equation with freezing and thawing"

      pde(wat)%solution_name(1) = "press_head" 
      pde(wat)%solution_name(2) = "h [L]" 

      pde(wat)%flux_name(1) = "flux"  
      pde(wat)%flux_name(2) = "Darcian flow [L.T^{-1}]"
      
      pde(wat)%print_mass = .true.
      deallocate(pde(wat)%mass)
      
      allocate(pde(wat)%mass_name(4,2))
      allocate(pde(wat)%mass(4))

      pde(wat)%mass_name(1,1) = "theta_tot"
      pde(wat)%mass_name(1,2) = "theta_tot [-]"
      
      pde(wat)%mass_name(2,1) = "theta_i"
      pde(wat)%mass_name(2,2) = "theta_i [-]"
      
      pde(wat)%mass_name(3,1) = "theta_l"
      pde(wat)%mass_name(3,2) = "theta_l [-]"
      
      pde(wat)%mass_name(4,1) = "h_l"
      pde(wat)%mass_name(4,2) = "h_l [L]"
      

      pde(wat)%mass(1)%val => vangen_fr
      
      pde(wat)%mass(2)%val => thetai
      
      pde(wat)%mass(3)%val => thetal
      
      pde(wat)%mass(4)%val => hl
      ! allocate porosity as mass(3)?

    ! pointers for heat flow model
      pde(heat_proc)%problem_name(1) = "heat"
      pde(heat_proc)%problem_name(2) = "Heat conduction equation with convection"

      pde(heat_proc)%solution_name(1) = "temperature" 
      pde(heat_proc)%solution_name(2) = "T " 

      pde(heat_proc)%flux_name(1) = "heat_flux"  
      pde(heat_proc)%flux_name(2) = "heat flux [W.L-2]"
      
      allocate(pde(heat_proc)%mass_name(0,2))
      pde(heat_proc)%print_mass = .false.
      
      pde(heat_proc)%pde_fnc(wat)%elasticity => capacityTh
      
      pde(heat_proc)%pde_fnc(heat_proc)%elasticity => capacityTT
      
      pde(heat_proc)%pde_fnc(heat_proc)%dispersion => diffTT
      
      pde(heat_proc)%pde_fnc(heat_proc)%convection => convectTT
      
      pde(heat_proc)%flux => heat_flux_freeze
      
      pde(heat_proc)%initcond => temp_initcond 
      
      do i=lbound(pde(heat_proc)%bc,1), ubound(pde(heat_proc)%bc,1)
        select case(pde(heat_proc)%bc(i)%code)
          case(1)
            pde(heat_proc)%bc(i)%value_fnc => heat_dirichlet
          case(2)
            pde(heat_proc)%bc(i)%value_fnc => heat_neumann
          case(0)
            pde(heat_proc)%bc(i)%value_fnc => re_null_bc
          case(3)
            pde(heat_proc)%bc(i)%value_fnc => freeze_coolant_bc
           case(4)
            pde(heat_proc)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
        end select
      end do  
        
    end subroutine frz_pointers
  
end module freeze_pointers
