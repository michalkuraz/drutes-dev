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
          call LTNE_reader(pde(1))
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from LTNE_pointers::frz_pointers"
          error stop
      end select
      
      call RE_totheadbc(pde(1))

      
      pde(1)%getval => getval_retotltne
      
    ! pointers for water flow model
      pde(1)%pde_fnc(1)%elasticity => capacityhh
            
      pde(1)%pde_fnc(1)%dispersion => diffhh
            
      pde(1)%flux => all_fluxes_LTNE
      do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
        select case(pde(1)%bc(i)%code)
          case(6)
            pde(1)%bc(i)%value_fnc => Dirichlet_mass_bc
        end select
      end do

      pde(1)%initcond => wat_init
      
      rwcap => vangen_elast_LTNE
      
      pde(1)%problem_name(1) = "RE_LTNE_thaw"
      pde(1)%problem_name(2) = "Richards' equation with freezing and thawing"

      pde(1)%solution_name(1) = "press_head" 
      pde(1)%solution_name(2) = "h [L]" 

      pde(1)%flux_name(1) = "flux"  
      pde(1)%flux_name(2) = "Darcian flow [L.T^{-1}]"
      
      pde(1)%print_mass = .true.
      deallocate(pde(1)%mass)
      
      allocate(pde(1)%mass_name(5,2))
      allocate(pde(1)%mass(5))

      pde(1)%mass_name(1,1) = "theta_tot"
      pde(1)%mass_name(1,2) = "theta_tot [-]"
      
      pde(1)%mass_name(2,1) = "theta_i"
      pde(1)%mass_name(2,2) = "theta_i [-]"
      
      pde(1)%mass_name(3,1) = "theta_l"
      pde(1)%mass_name(3,2) = "theta_l [-]"
      
      pde(1)%mass_name(4,1) = "h_l"
      pde(1)%mass_name(4,2) = "h_l [L]"
      
      pde(1)%mass_name(5,1) = "T_m"
      pde(1)%mass_name(5,2) = "T_m [deg C]"

      pde(1)%mass(1)%val => vangen_ltne
      
      pde(1)%mass(2)%val => thetai
      
      pde(1)%mass(3)%val => thetal
      
      pde(1)%mass(4)%val => hl
      
      pde(1)%mass(5)%val => T_m

      ! allocate porosity as mass(3)?

    ! pointers for heat flow model
      pde(2)%problem_name(1) = "heat_liquid"
      pde(2)%problem_name(2) = "Heat conduction equation with convection"

      pde(2)%solution_name(1) = "T_liquid" 
      pde(2)%solution_name(2) = "T_l " 

      pde(2)%flux_name(1) = "heat_l_flux"  
      pde(2)%flux_name(2) = "heat liquid flux [W.L-2]"
      
      allocate(pde(2)%mass_name(0,2))
      pde(2)%print_mass = .false.
      
      pde(2)%pde_fnc(1)%elasticity => capacityTlh
      
      pde(2)%pde_fnc(2)%elasticity => capacityTlTl
      
      pde(2)%pde_fnc(2)%dispersion => diffTlTl
      
      pde(2)%pde_fnc(2)%convection => convectTlTl
      
      pde(2)%pde_fnc(3)%reaction => qsl_pos
      pde(2)%pde_fnc(2)%reaction => qsl_neg

      pde(2)%flux => heat_flux_l_LTNE
      
      pde(2)%initcond => temp_l_initcond 
      
      do i=lbound(pde(2)%bc,1), ubound(pde(2)%bc,1)
        select case(pde(2)%bc(i)%code)
          case(1)
            pde(2)%bc(i)%value_fnc => heat_dirichlet
          case(2)
            pde(2)%bc(i)%value_fnc => heat_neumann
          case(0)
            pde(2)%bc(i)%value_fnc => re_null_bc
          case(3)
            pde(2)%bc(i)%value_fnc => LTNE_coolant_bc
        end select
      end do  
        
        
      !pointers for solid heat flow model
      pde(3)%problem_name(1) = "heat_solid"
      pde(3)%problem_name(2) = "Heat conduction equation with convection"

      pde(3)%solution_name(1) = "T_solid" 
      pde(3)%solution_name(2) = "T_s " 

      pde(3)%flux_name(1) = "heat_s_flux"  
      pde(3)%flux_name(2) = "heat solid flux [W.L-2]"
      
      allocate(pde(3)%mass_name(0,2))
      pde(3)%print_mass = .false.
            
      pde(3)%pde_fnc(3)%elasticity => capacityTsTs
      
      pde(3)%pde_fnc(3)%dispersion => diffTsTs
            
      pde(3)%flux => heat_flux_s_LTNE
      
      pde(3)%initcond => temp_s_initcond 
      
      pde(3)%pde_fnc(3)%reaction => qsl_neg
      pde(3)%pde_fnc(2)%reaction => qsl_pos
      
      do i=lbound(pde(3)%bc,1), ubound(pde(2)%bc,1)
        select case(pde(3)%bc(i)%code)
          case(1)
            pde(3)%bc(i)%value_fnc => heat_dirichlet
          case(2)
            pde(3)%bc(i)%value_fnc => heat_neumann
          case(0)
            pde(3)%bc(i)%value_fnc => re_null_bc
          case(3)
            pde(3)%bc(i)%value_fnc => LTNE_coolant_bc
        end select
      end do  
           
        
    end subroutine LTNE_points
  
  

end module LTNE_pointers
