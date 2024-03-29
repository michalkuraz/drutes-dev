module freeze_pointers
  public :: freeze_processes, frz_pointers
  
  contains
  
    subroutine freeze_processes(processes)
      use typy
      use globals
      use freeze_globs
      use freeze_read
      
      integer(kind=ikind), intent(out) :: processes
      call read_frrate()
      select case (drutes_config%name)
        case ("freeze")
            processes = 2
         case ("LTNE")
            processes = 3
         case ("ICENE")
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
      use debug_tools
      
      integer(kind=ikind) :: i
      
      select case (drutes_config%name)
        case ("freeze")
          call freeze_reader(pde(wat))
        case ("LTNE")
          call freeze_reader(pde(wat))
        case ("ICENE")
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
      select case (drutes_config%name)
        case ("freeze")
          pde(wat)%pde_fnc(heat_proc)%elasticity => capacityhT
        case ("LTNE")
          pde(wat)%pde_fnc(heat_proc)%elasticity => capacityhT
        case ("ICENE")
          !pde(wat)%pde_fnc(ice)%elasticity => capacityhice
           pde(wat)%pde_fnc(wat)%zerord => zerordhice

        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::frz_pointers"
          error stop 
      end select

      pde(wat)%pde_fnc(wat)%dispersion => diffhh
      pde(wat)%pde_fnc(heat_proc)%dispersion => diffhT
      pde(wat)%flux => all_fluxes
      do i=lbound(pde(wat)%bc,1), ubound(pde(wat)%bc,1)
        select case(pde(wat)%bc(i)%code)
          case(7)
            pde(wat)%bc(i)%value_fnc => Dirichlet_mass_bc
          case(8)
            pde(wat)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
          case(9)
            pde(wat)%bc(i)%value_fnc => Infiltration_bc
          case(10)
            pde(wat)%bc(i)%value_fnc => Infiltration_bc2
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
      deallocate(pde(wat)%mass)

      select case (drutes_config%name)
        case ("freeze")
            allocate(pde(wat)%mass_name(5,2))
            allocate(pde(wat)%mass(5))
            
            pde(wat)%mass_name(4,1) = "theta_i"
            pde(wat)%mass_name(4,2) = "theta_i [-]"
            
            pde(wat)%mass_name(5,1) = "theta_i_eqwat"
            pde(wat)%mass_name(5,2) = "thi eq. wat [-]"
            
            pde(wat)%mass(4)%val => thetai
            pde(wat)%mass(5)%val => thetai_wat_eq
            pde(wat)%mass_name(1,1) = "theta_tot"
            pde(wat)%mass_name(1,2) = "theta_tot [-]"
      
            pde(wat)%mass_name(2,1) = "theta_l"
            pde(wat)%mass_name(2,2) = "theta_l [-]"
      
            pde(wat)%mass_name(3,1) = "h_l"
            pde(wat)%mass_name(3,2) = "h_l [L]"
      
            pde(wat)%mass(1)%val => vangen_fr
            
            pde(wat)%mass(2)%val => thetal
      
            pde(wat)%mass(3)%val => hl
      
        case ("LTNE")
            allocate(pde(wat)%mass_name(5,2))
            allocate(pde(wat)%mass(5))

            pde(wat)%mass_name(4,1) = "T_m"
            pde(wat)%mass_name(4,2) = "T_m [deg C]"
            
            pde(wat)%mass_name(5,1) = "theta_i"
            pde(wat)%mass_name(5,2) = "theta_i [-]"
            
            pde(wat)%mass(4)%val => T_m        

            pde(wat)%mass(5)%val => thetai
            pde(wat)%mass_name(1,1) = "theta_tot"
            pde(wat)%mass_name(1,2) = "theta_tot [-]"
      
            pde(wat)%mass_name(2,1) = "theta_l"
            pde(wat)%mass_name(2,2) = "theta_l [-]"
      
            pde(wat)%mass_name(3,1) = "h_l"
            pde(wat)%mass_name(3,2) = "h_l [L]"
      
            pde(wat)%mass(1)%val => vangen_fr
            
            pde(wat)%mass(2)%val => thetal
      
           pde(wat)%mass(3)%val => hl
        case ("ICENE")
            allocate(pde(wat)%mass_name(6,2))
            allocate(pde(wat)%mass(6))
            
                  
            pde(wat)%mass_name(1,1) = "theta_l"
            pde(wat)%mass_name(1,2) = "theta_l [-]"

            pde(wat)%mass_name(2,1) = "theta_cl"
            pde(wat)%mass_name(2,2) = "theta_cl [-]"
      
            pde(wat)%mass_name(3,1) = "theta_i"
            pde(wat)%mass_name(3,2) = "theta_i [-]"
            
            pde(wat)%mass_name(4,1) = "theta_i_eqwat"
            pde(wat)%mass_name(4,2) = "thi eq. wat [-]"

            pde(wat)%mass_name(5,1) = "icerate"
            pde(wat)%mass_name(5,2) = "vf [1/T]"
            
            pde(wat)%mass_name(6,1) = "theta_i_eq"
            pde(wat)%mass_name(6,2) = "ice_eq"
            
            pde(wat)%mass(1)%val => vangen_fr
                  
            pde(wat)%mass(2)%val => theta_cl
            
            pde(wat)%mass(3)%val => thetai
            pde(wat)%mass(4)%val => thetai_wat_eq
            
            pde(wat)%mass(5)%val => icerate
            
            pde(wat)%mass(6)%val => thetai_eq
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::frz_pointers"
          error stop
      end select
      


      
    ! pointers for heat flow model
      pde(heat_proc)%problem_name(1) = "heat"
      pde(heat_proc)%problem_name(2) = "Heat conduction equation with convection"

      pde(heat_proc)%solution_name(1) = "temperature" 
      pde(heat_proc)%solution_name(2) = "T" 

      pde(heat_proc)%flux_name(1) = "flux"  
      pde(heat_proc)%flux_name(2) = "heat flux [W.L-2]"
      
      allocate(pde(heat_proc)%mass_name(0,2))
      pde(heat_proc)%print_mass = .false.
      
      pde(heat_proc)%pde_fnc(wat)%elasticity => capacityTh
      
      pde(heat_proc)%pde_fnc(heat_proc)%elasticity => capacityTT
      
      pde(heat_proc)%pde_fnc(heat_proc)%dispersion => diffTT

      pde(heat_proc)%pde_fnc(wat)%dispersion => diffTh

      pde(heat_proc)%pde_fnc(heat_proc)%convection => convectTT
      
      pde(heat_proc)%flux => heat_flux_freeze
      
      pde(heat_proc)%initcond => temp_initcond 
      
      select case (drutes_config%name)
        case ("freeze")
          CONTINUE
        case ("LTNE")
          CONTINUE
        case ("ICENE")
          !pde(heat_proc)%pde_fnc(ice)%elasticity => capacityTice
          pde(heat_proc)%pde_fnc(heat_proc)%zerord => zerordTice
          pde(heat_proc)%pde_fnc(heat_proc)%reaction => reactTice

        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::frz_pointers"
          error stop 
      end select
      
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
          !pointers for solid heat flow model

      select case (drutes_config%name)
        case ("freeze", "ICENE")
         CONTINUE
        case ("LTNE")
          pde(heat_solid)%problem_name(1) = "heat_solid"
          pde(heat_solid)%problem_name(2) = "Heat conduction equation with convection"

          pde(heat_solid)%solution_name(1) = "Temperature" 
          pde(heat_solid)%solution_name(2) = "T_s " 

          pde(heat_solid)%flux_name(1) = "flux"  
          pde(heat_solid)%flux_name(2) = "heat solid flux [W.L-2]"
        
          allocate(pde(heat_solid)%mass_name(0,2))
          pde(heat_solid)%print_mass = .false.
                
          pde(heat_solid)%pde_fnc(heat_solid)%elasticity => capacityTsTs
        
          pde(heat_solid)%pde_fnc(heat_solid)%dispersion => diffTsTs
                
          pde(heat_solid)%flux => heat_flux_s_LTNE
        
          pde(heat_solid)%initcond => temp_s_initcond 
                    
          pde(heat_solid)%pde_fnc(heat_solid)%reaction => qsl_neg
          pde(heat_solid)%pde_fnc(heat_proc)%reaction => qsl_pos
          
          pde(heat_proc)%pde_fnc(heat_solid)%reaction => qsl_pos
          pde(heat_proc)%pde_fnc(heat_proc)%reaction => qsl_neg

        
          do i=lbound(pde(heat_solid)%bc,1), ubound(pde(heat_proc)%bc,1)
            select case(pde(heat_solid)%bc(i)%code)
              case(1)
                pde(heat_solid)%bc(i)%value_fnc => heat_dirichlet
              case(2)
                pde(heat_solid)%bc(i)%value_fnc => heat_neumann
              case(0)
                pde(heat_solid)%bc(i)%value_fnc => re_null_bc
              case(3)
                pde(heat_solid)%bc(i)%value_fnc => freeze_coolant_bc
              case(4)
                pde(heat_solid)%bc(i)%value_fnc => Dirichlet_Neumann_switch_bc
              case(42)
                pde(heat_solid)%bc(i)%value_fnc => freeze_coolant_bc_bot
            end select
        end do 
        
        
        
        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::frz_pointers"
          error stop
      end select
      
      select case (drutes_config%name)
        case ("freeze")
          CONTINUE
        case ("LTNE")
          CONTINUE
        case ("ICENE")
          pde(ice)%problem_name(1) = "ice"
          pde(ice)%problem_name(2) = "vol. ice"

          pde(ice)%solution_name(1) = "content" 
          pde(ice)%solution_name(2) = "i" 

          pde(ice)%flux_name(1) = "flux"  
          pde(ice)%flux_name(2) = "ice flux [L.T^{-1}]"

          allocate(pde(ice)%mass_name(0,2))
          pde(ice)%print_mass = .false.
          pde(ice)%pde_fnc(ice)%elasticity => capacityiceice
          pde(ice)%pde_fnc(ice)%zerord => icerate
          do i=lbound(pde(ice)%bc,1), ubound(pde(ice)%bc,1)
            select case(pde(ice)%bc(i)%code)
              case(1)
                pde(ice)%bc(i)%value_fnc => heat_dirichlet
              case(2)
                pde(ice)%bc(i)%value_fnc => heat_neumann
              case(0)
                pde(ice)%bc(i)%value_fnc => re_null_bc
            end select
        end do 
        pde(ice)%initcond => ice_initcond 
        pde(ice)%print_mass = .false.

        case default
          print *, "procedure called when unexpected problem name"
          print *, "exited from freeze_pointers::frz_pointers"
          error stop 
      end select

    end subroutine frz_pointers
  
end module freeze_pointers
