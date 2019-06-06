module LTNE_pointers
  public :: LTNE_processes, LTNE_points
  contains
  
    subroutine LTNE_processes(processes)
      use typy
      use globals
      
      integer(kind=ikind), intent(out) :: processes
      
      select case (drutes_config%name)
        case ("LTNE")
          processes = 3
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from LTNE_pointers::LTNE_processes"
          error stop
      end select
        
    
    end subroutine LTNE_processes
  
  
  
    subroutine LTNE_points()
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_constitutive
      use re_total
      use re_pointers
      use LTNE_globs
      use LTNE_fnc
      use LTNE_helper
      use LTNE_read
      use heat_pointers
      use heat_fnc
      
      
      integer(kind=ikind) :: i

      
      select case (drutes_config%name)
        case ("LTNE")
          call LTNE_reader(pde(wat))
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from LTNE_pointers::frz_pointers"
          error stop
      end select
      
      call RE_totheadbc(pde(wat))

      
      pde(wat)%getval => getval_retotltne
      
    ! pointers for water flow model
      pde(wat)%pde_fnc(wat)%elasticity => capacityhh
            
      pde(wat)%pde_fnc(wat)%dispersion => diffhh
      pde(wat)%pde_fnc(heat_proc)%dispersion => diffhT

      pde(wat)%flux => all_fluxes_LTNE
      do i=lbound(pde(wat)%bc,1), ubound(pde(wat)%bc,1)
        select case(pde(wat)%bc(i)%code)
          case(6)
            pde(wat)%bc(i)%value_fnc => Dirichlet_mass_bc
          case(7)
            pde(wat)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
        end select
      end do

      pde(wat)%initcond => wat_init
      
      rwcap => vangen_elast_LTNE
      
      
      pde(wat)%problem_name(1) = "RE_LTNE_thaw"
      pde(wat)%problem_name(2) = "Richards' equation with freezing and thawing"

      pde(wat)%solution_name(1) = "press_head" 
      pde(wat)%solution_name(2) = "h [L]" 

      pde(wat)%flux_name(1) = "flux"  
      pde(wat)%flux_name(2) = "Darcian flow [L.T^{-1}]"
      
      pde(wat)%print_mass = .true.
      deallocate(pde(wat)%mass)
      
      allocate(pde(wat)%mass_name(5,2))
      allocate(pde(wat)%mass(5))

      pde(wat)%mass_name(1,1) = "theta_tot"
      pde(wat)%mass_name(1,2) = "theta_tot [-]"
      
      pde(wat)%mass_name(2,1) = "theta_i"
      pde(wat)%mass_name(2,2) = "theta_i [-]"
      
      pde(wat)%mass_name(3,1) = "theta_l"
      pde(wat)%mass_name(3,2) = "theta_l [-]"
      
      pde(wat)%mass_name(4,1) = "h_l"
      pde(wat)%mass_name(4,2) = "h_l [L]"
      
      pde(wat)%mass_name(5,1) = "T_m"
      pde(wat)%mass_name(5,2) = "T_m [deg C]"

      pde(wat)%mass(1)%val => vangen_ltne
      
      pde(wat)%mass(2)%val => thetai
      
      pde(wat)%mass(3)%val => thetal
      
      pde(wat)%mass(4)%val => hl
      
      pde(wat)%mass(5)%val => T_m

      ! allocate porosity as mass(3)?

    ! pointers for heat flow model
      pde(heat_proc)%problem_name(1) = "heat_liquid"
      pde(heat_proc)%problem_name(2) = "Heat conduction equation with convection"

      pde(heat_proc)%solution_name(1) = "T_liquid" 
      pde(heat_proc)%solution_name(2) = "T_l " 

      pde(heat_proc)%flux_name(1) = "heat_l_flux"  
      pde(heat_proc)%flux_name(2) = "heat liquid flux [W.L-2]"
      
      allocate(pde(heat_proc)%mass_name(0,2))
      pde(heat_proc)%print_mass = .false.
      
      pde(heat_proc)%pde_fnc(wat)%elasticity => capacityTlh
      
      pde(heat_proc)%pde_fnc(heat_proc)%elasticity => capacityTlTl
      
      pde(heat_proc)%pde_fnc(heat_proc)%dispersion => diffTlTl
      
      pde(heat_proc)%pde_fnc(heat_proc)%convection => convectTlTl
      
      pde(heat_proc)%pde_fnc(heat_solid)%reaction => qsl_pos
      pde(heat_proc)%pde_fnc(heat_proc)%reaction => qsl_neg

      pde(heat_proc)%flux => heat_flux_l_LTNE
      
      pde(heat_proc)%initcond => temp_l_initcond 
      
      do i=lbound(pde(heat_proc)%bc,1), ubound(pde(heat_proc)%bc,1)
        select case(pde(heat_proc)%bc(i)%code)
          case(1)
            pde(heat_proc)%bc(i)%value_fnc => heat_dirichlet
          case(2)
            pde(heat_proc)%bc(i)%value_fnc => heat_neumann
          case(0)
            pde(heat_proc)%bc(i)%value_fnc => re_null_bc
          case(3)
            pde(heat_proc)%bc(i)%value_fnc => ltne_coolant_bc
          case(4)
            pde(heat_proc)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
        end select
      end do  
        
        
      !pointers for solid heat flow model
      pde(heat_solid)%problem_name(1) = "heat_solid"
      pde(heat_solid)%problem_name(2) = "Heat conduction equation with convection"

      pde(heat_solid)%solution_name(1) = "T_solid" 
      pde(heat_solid)%solution_name(2) = "T_s " 

      pde(heat_solid)%flux_name(1) = "heat_s_flux"  
      pde(heat_solid)%flux_name(2) = "heat solid flux [W.L-2]"
      
      allocate(pde(heat_solid)%mass_name(0,2))
      pde(heat_solid)%print_mass = .false.
            
      pde(heat_solid)%pde_fnc(heat_solid)%elasticity => capacityTsTs
      
      pde(heat_solid)%pde_fnc(heat_solid)%dispersion => diffTsTs
            
      pde(heat_solid)%flux => heat_flux_s_LTNE
      
      pde(heat_solid)%initcond => temp_s_initcond 
      
      pde(heat_solid)%pde_fnc(heat_solid)%reaction => qsl_neg
      pde(heat_solid)%pde_fnc(heat_proc)%reaction => qsl_pos
      
      do i=lbound(pde(heat_solid)%bc,1), ubound(pde(heat_proc)%bc,1)
        select case(pde(heat_solid)%bc(i)%code)
          case(1)
            pde(heat_solid)%bc(i)%value_fnc => heat_dirichlet
          case(2)
            pde(heat_solid)%bc(i)%value_fnc => heat_neumann
          case(0)
            pde(heat_solid)%bc(i)%value_fnc => re_null_bc
          case(3)
            pde(heat_solid)%bc(i)%value_fnc => ltne_coolant_bc
          case(4)
            pde(heat_solid)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
        end select
      end do  
           
        
    end subroutine LTNE_points
  
  

end module LTNE_pointers
