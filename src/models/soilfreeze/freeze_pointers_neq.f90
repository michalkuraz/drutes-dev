module freeze_pointers_neq
  public :: freeze_processes_neq, frz_pointers_neq
  
  contains
  
    subroutine freeze_processes_neq(processes)
      use typy
      use globals
      use freeze_globs
      use freeze_read_neq
      
      integer(kind=ikind), intent(out) :: processes
      call read_frrate_neq()
      select case (drutes_config%name)
        case ("freezeneq")
            processes = 3
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::freeze_processes"
          error stop
      end select
        
    
    end subroutine freeze_processes_neq
  
  
  
    subroutine frz_pointers_neq()
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
      use freeze_read_neq
      use heat_pointers
      use heat_reader
      use heat_fnc
      use debug_tools
      
      integer(kind=ikind) :: i
      
      select case (drutes_config%name)
        case ("freezeneq")
          call freeze_reader_neq(pde(wat))
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::frz_pointers"
          error stop 
      end select


      call RE_totheadbc(pde(wat))
      pde(wat)%getval => getval_retotfr
      
    ! pointers for water flow model
      pde(wat)%pde_fnc(wat)%elasticity => capacityhh_neq

      pde(wat)%pde_fnc(ice)%elasticity => capacityhice
      
      pde(wat)%pde_fnc(wat)%dispersion => diffhh
      
      pde(wat)%pde_fnc(heat_proc)%dispersion => diffhT_neq
      
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
            
      pde(wat)%problem_name(1) = "RE_freeze_thaw"
      pde(wat)%problem_name(2) = "Richards' equation with freezing and thawing"

      pde(wat)%solution_name(1) = "press_head" 
      pde(wat)%solution_name(2) = "h [L]" 

      pde(wat)%flux_name(1) = "flux"  
      pde(wat)%flux_name(2) = "Darcian flow [L.T^{-1}]"
      pde(wat)%print_mass = .true.
      
      pde(wat)%print_mass = .false.
      allocate(pde(wat)%mass_name(0,2))


    ! pointers for heat flow model
      pde(heat_proc)%problem_name(1) = "heat"
      pde(heat_proc)%problem_name(2) = "Heat conduction equation with convection"

      pde(heat_proc)%solution_name(1) = "temperature" 
      pde(heat_proc)%solution_name(2) = "T" 

      pde(heat_proc)%flux_name(1) = "flux"  
      pde(heat_proc)%flux_name(2) = "heat flux [W.L-2]"

      allocate(pde(heat_proc)%mass_name(0,2))
      pde(heat_proc)%print_mass = .false.

      pde(heat_proc)%pde_fnc(wat)%elasticity => capacityTh_neq
      
      pde(heat_proc)%pde_fnc(heat_proc)%elasticity => capacityTT_neq

      pde(heat_proc)%pde_fnc(ice)%elasticity => capacityTice
      
      pde(heat_proc)%pde_fnc(heat_proc)%dispersion => diffTT
      
      pde(heat_proc)%pde_fnc(wat)%dispersion => diffTh
      
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
          case(42)
            pde(heat_proc)%bc(i)%value_fnc => freeze_coolant_bc_bot
        end select
      end do  


        allocate(pde(ice)%mass_name(0,2))
        pde(ice)%print_mass = .false.
        pde(ice)%problem_name(1) = "ice"
        pde(ice)%problem_name(2) = "Ice with freezing rate"
            
        pde(ice)%solution_name(1) = "ice_content" 
        pde(ice)%solution_name(2) = "theta_i " 

        pde(ice)%pde_fnc(ice)%elasticity => capacity_ice
        pde(ice)%pde_fnc(ice)%zerord => ice_rate
        

        do i=lbound(pde(ice)%bc,1), ubound(pde(ice)%bc,1)
          select case(pde(ice)%bc(i)%code)
            case(1)
                pde(ice)%bc(i)%value_fnc => heat_dirichlet
            case(2)
                pde(ice)%bc(i)%value_fnc => heat_neumann
          end select
        end do

        pde(ice)%initcond => ice_initcond

    end subroutine frz_pointers_neq

end module freeze_pointers_neq
