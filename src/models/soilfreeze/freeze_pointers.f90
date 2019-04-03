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
          processes = 3
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
          call freeze_reader(pde(1))
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::frz_pointers"
          error stop
      end select
      
      call RE_totheadbc(pde(1))
      pde(1)%getval => getval_retotfr
      
    ! pointers for water flow model
      pde(1)%pde_fnc(1)%elasticity => capacityhh
      
      pde(1)%pde_fnc(2)%elasticity => capacityhT
      
      pde(1)%pde_fnc(1)%dispersion => diffhh
      
      pde(1)%pde_fnc(2)%dispersion => diffhT
      
      pde(1)%flux => all_fluxes
      do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
        select case(pde(1)%bc(i)%code)
          case(6)
            pde(1)%bc(i)%value_fnc => Dirichlet_mass_bc
          case(7)
            pde(1)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
        end select
      end do
      pde(1)%initcond => wat_initcond
      
      rwcap => vangen_elast_fr
      
      pde(1)%problem_name(1) = "RE_freeze_thaw"
      pde(1)%problem_name(2) = "Richards' equation with freezing and thawing"

      pde(1)%solution_name(1) = "press_head" 
      pde(1)%solution_name(2) = "h [L]" 

      pde(1)%flux_name(1) = "flux"  
      pde(1)%flux_name(2) = "Darcian flow [L.T^{-1}]"
      
      pde(1)%print_mass = .true.
      deallocate(pde(1)%mass)
      
      allocate(pde(1)%mass_name(4,2))
      allocate(pde(1)%mass(4))

      pde(1)%mass_name(1,1) = "theta_tot"
      pde(1)%mass_name(1,2) = "theta_tot [-]"
      
      pde(1)%mass_name(2,1) = "theta_i"
      pde(1)%mass_name(2,2) = "theta_i [-]"
      
      pde(1)%mass_name(3,1) = "theta_l"
      pde(1)%mass_name(3,2) = "theta_l [-]"
      
      pde(1)%mass_name(4,1) = "h_l"
      pde(1)%mass_name(4,2) = "h_l [L]"
      

      pde(1)%mass(1)%val => vangen_fr
      
      pde(1)%mass(2)%val => thetai
      
      pde(1)%mass(3)%val => thetal
      
      pde(1)%mass(4)%val => hl
      ! allocate porosity as mass(3)?

    ! pointers for heat flow model
      pde(2)%problem_name(1) = "heat"
      pde(2)%problem_name(2) = "Heat conduction equation with convection"

      pde(2)%solution_name(1) = "temperature" 
      pde(2)%solution_name(2) = "T " 

      pde(2)%flux_name(1) = "heat_flux"  
      pde(2)%flux_name(2) = "heat flux [W.L-2]"
      
      allocate(pde(2)%mass_name(0,2))
      pde(2)%print_mass = .false.
      
      pde(2)%pde_fnc(1)%elasticity => capacityTh
      
      pde(2)%pde_fnc(2)%elasticity => capacityTT
      
      pde(2)%pde_fnc(2)%dispersion => diffTT
      
      pde(2)%pde_fnc(2)%convection => convectTT
      
      pde(2)%flux => heat_flux_freeze
      
      pde(2)%initcond => temp_initcond 
      
      do i=lbound(pde(2)%bc,1), ubound(pde(2)%bc,1)
        select case(pde(2)%bc(i)%code)
          case(1)
            pde(2)%bc(i)%value_fnc => heat_dirichlet
          case(2)
            pde(2)%bc(i)%value_fnc => heat_neumann
          case(0)
            pde(2)%bc(i)%value_fnc => re_null_bc
          case(3)
            pde(2)%bc(i)%value_fnc => freeze_coolant_bc
           case(4)
            pde(2)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
        end select
      end do  
        
    end subroutine frz_pointers
  
end module freeze_pointers
