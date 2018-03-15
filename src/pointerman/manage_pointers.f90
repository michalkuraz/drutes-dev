module manage_pointers
  public :: set_pointers

  contains


    !> set pointers for the entire problem, except the file read pointers
    subroutine set_pointers()
      use typy
      use global_objs
      use pde_objs
      use globals
      use core_tools
      use capmat
      use fem_tools
      use feminittools
      use fem
      use femmat
      use schwarz_dd
      use schwarz_dd2subcyc
      use solver_interfaces
      use debug_tools
      use bousspointers
      use re_pointers
      use ADE_pointers
      use Re_dual_pointers
      use heat_pointers
      use drutes_init

      integer(kind=ikind) :: i, processes
      
      

      select case(cut(drutes_config%name))
        case("REstd")
          write(unit=drutes_config%fullname, fmt=*) "Richards' equation, in", drutes_config%dimen, &
            "D, in pressure head form."
        
        
          call RE_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)
          
          call RE_std(pde(1))
        case("RE")
          write(unit=drutes_config%fullname, fmt=*) "Richards' equation, in", drutes_config%dimen, &
            "D, in total hydraulic head form."
          call RE_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)

          call REstdH(pde(1))
        
              
        case("boussi")   
           write(unit=drutes_config%fullname, fmt=*) " Boussinesq equation for hillslope runoff", &
                 "(1D approximation of groundwater flow)."
            
           call bouss_processes(pde_common%processes)
           call pde_constructor(pde_common%processes)      
           call boussi(pde(1))
                 
              
        case("ADE") 
        
          call ade_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)
          write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation in", &
                 drutes_config%dimen, "D"
          call ade()

        case("Re_dual")
            write(unit=drutes_config%fullname, fmt=*) " Richards equation ", &
              "in total hydraulic head form for dual (fracture and matrix) medium"	
             pde_common%processes = 2
             call pde_constructor(pde_common%processes)
             call RE_matrix()
             call RE_fracture()  
      
    
        case("REtest")
          write(unit=drutes_config%fullname, fmt=*) "DRUtES debugs itself"
          pde_common%processes = 3
          call pde_constructor(3_ikind)
          do i=1, 3
            call RE_std(pde(i))
          end do
	  
        case("heat")
        
          call heat_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)
          write(unit=drutes_config%fullname, fmt=*) "DRUtES solves heat conduction with convection"
          call heat(pde(:))
	  
        case default
          print *, "your new model: ", trim(drutes_config%name), " requires pointer linking"
          print *, "exited from manage_pointers::set_pointers"
          ERROR stop

	  
      end select

      select case(drutes_config%dimen)
        case(1)
            solve_matrix => LDU_face
!        	    solve_matrix => CG_normal_face
        case(2)
      !           solve_matrix => pcg
      ! 	    solve_matrix => LDU_face
            solve_matrix => CG_normal_face
      ! 	    solve_matrix => sparse_gem_pig_AtA
      ! 	    solve_matrix => jacobi_face
      end select
      
      select case (drutes_config%it_method)
        case(0) 
          pde_common%treat_pde => solve_picard
        case(1)
          pde_common%treat_pde => schwarz_picard
        case(2)
          pde_common%treat_pde => schwarz_subcyc
      end select
          
      select case(pde_common%timeint_method)
        case(0)
          pde_common%time_integ => steady_state_int
        case(1)
          pde_common%time_integ => impl_euler_np_diag
        case(2)
          pde_common%time_integ => impl_euler_np_nondiag
      end select


  end subroutine set_pointers


  

end module manage_pointers
