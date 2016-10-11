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
      use Re_dual_pointers ! added J

      integer(kind=ikind) :: i
      select case(adjustl(trim(drutes_config%name)))
	case("RE_std")
	  write(unit=drutes_config%fullname, fmt=*) "Richards' equation, in", drutes_config%dimen, &
	  "D, in pressure head form."
	  call RE_std(pde(1))
	case("REstdH")
	  write(unit=drutes_config%fullname, fmt=*) "Richards' equation, in", drutes_config%dimen, &
	  "D, in total hydraulic head form."
	  call REstdH(pde(1))
	case("RE_rot")
	  write(unit=drutes_config%fullname, fmt=*) "Richards' equation, in", drutes_config%dimen, &
	  "D, in cylindric coordinates and axysimmetric flow, in pressure head form."
	  call RE_rot(pde(1))
	case("RErotH")
	  write(unit=drutes_config%fullname, fmt=*) "Richards' equation, in", drutes_config%dimen, &
	  "D, in cylindric coordinates and axysimmetric flow, in total hydraulic head form."
	  call RErotH(pde(1))
	
	      
	case("boussi")   
	   write(unit=drutes_config%fullname, fmt=*) " Boussinesq equation for hillslope runoff", &
           "(1D approximation of groundwater flow)."
           call boussi(pde(1))
           
        
	case("ADEstd") 
	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation in", &
           drutes_config%dimen, "D, convection is specified in input files, equilibrium sorption"
           call ade(pde(1))

      case("ADEstd_kinsorb")
      
      	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation in", &
           drutes_config%dimen, "D, convection is specified in input files, kinetic sorption"
           call ade(pde(1))
           call adekinsorb(pde(2))
           
      case("ADE_RE_std")
     	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation", &
     	  "  Richards equation in pressure head form in", &
           drutes_config%dimen, "D, convection is computed, equilibrium sorption"	
           
           call RE_std(pde(1))
           call ade(pde(2))
	
      case("ADE_REstdH")
     	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation", &
     	  " and Richards equation in total hydraulic head form in", &
           drutes_config%dimen, "D, convection is computed, equilibrium sorption"	
           
           call REstdH(pde(1))
           call ade(pde(2))


      case("ADE_RE_rot")
     	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation", &
     	  " and Richards equation in pressure head form", &
           "flow is axisymmetric, convection is computed, equilibrium sorption"	
           
           call RE_rot(pde(1))
           call ade(pde(2))
           
      case("ADE_RErotH")
     	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation", &
     	  " and Richards equation in total hydraulic head form", &
           "flow is axisymmetric, convection is computed, equilibrium sorption"	
           
           call RE_rot(pde(1))
           call ade(pde(2))           
	  
	  ! added J 13/6/16
	  case("Re_dual_totH")
	  		write(unit=drutes_config%fullname, fmt=*) " Richards equation ", &
     	  "in total hydraulic head form for dual (fracture and matrix) medium"	
         
           call RE_matrix(pde(1))
           call RE_fracture(pde(2))  


      case("ADE_RE_std_kinsorb")
     	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation", &
     	  " and Richards equation in pressure head form in", &
           drutes_config%dimen, "D, convection is computed, kinetic sorption"	
           
           call RE_std(pde(1))
           call ade(pde(2))
           call adekinsorb(pde(3))
	
      case("ADE_REstdH_kinsorb")
     	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation", &
     	  " and Richards equation in total hydraulic head form in", &
           drutes_config%dimen, "D, convection is computed, kinetic sorption"	
           
           call REstdH(pde(1))
           call ade(pde(2))
           call adekinsorb(pde(3))

      case("ADE_RE_rot_kinsorb")
     	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation", &
     	  " and Richards equation in pressure head form", &
           " flow is axisymmetric, convection is computed, kinetic sorption"	
           
           call RE_rot(pde(1))
           call ade(pde(2))
           call adekinsorb(pde(3))
           
      case("ADE_RErotH_kinsorb")
     	  write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation", &
     	  " and Richards equation in total hydraulic head form", &
           " flow is axisymmetric, convection is computed, kinetic sorption"	
           
           call RE_rot(pde(1))
           call ade(pde(2))
           call adekinsorb(pde(3))
		
	
       case("REtest")
	  write(unit=drutes_config%fullname, fmt=*) "DRUtES debugs itself"
	  do i=1, 3
	    call RE_std(pde(i))
	  end do

	  
    end select

      select case(drutes_config%dimen)
	case(1)
 	    solve_matrix => LDU_face
	case(2)
!  solve_matrix => pcg
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