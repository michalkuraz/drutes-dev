module manage_pointers
  public :: set_pointers
  public :: set_readers

  contains

    subroutine set_readers()
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_reader
      use modRE_reader
      use boussread
      use ADE_reader
      
      integer :: i


      select case(drutes_config%name)
        case("RE_std", "RE_rot", "REstdH", "RErotH")
          pde(1)%read_parameters => res_read
        case("RE_mod")
	  pde(1)%read_parameters => modre_read
	  pde(2)%read_parameters => modheat_read
	  pde(3)%read_parameters => modsolute_read
	  
	case("REtest")
	  do i=1, ubound(pde,1)
	    pde(i)%read_parameters => res_read
	  end do
	  
	case("boussi")
	  pde(1)%read_parameters => boussreader
	  
	case("ADEstd")
	  pde(1)%read_parameters => ADE_read
	  
	case("ADE_wr")
	  pde(1)%read_parameters => res_read
	  pde(2)%read_parameters => ADE_read
	  
        
        case default
          print *, "unsupported problem type, terminated from manage_pointers::set_readers"
      end select
        


    end subroutine set_readers


    !> set pointers for the entire problem, except the file read pointers
    subroutine set_pointers()
      use typy
      use global_objs
      use pde_objs
      use globals
      use re_constitutive
      use re_total
      use core_tools
      use capmat
      use fem_tools
      use feminittools
      use fem
      use femmat
      use schwarz_dd
      use schwarz_dd2subcyc
      use solver_interfaces
      use modRE_constitutive
      use modRE_junctions
      use debug_tools
      use boussfnc

      

      integer(kind=ikind) :: i,j
      type(integpnt_str) :: quadpnt
      
      do i=1, ubound(pde,1)
	pde(i)%getval => getvalp1
      end do


      select case(drutes_config%name)
       !general setting for Richards equation in all modes
	case("RE_std", "RE_rot", "RErotH", "REstdH")
          call domainswitch("m")
	  pde_common%nonlinear = .true.
	  if (drutes_config%fnc_method == 0) then
	    pde(1)%pde_fnc(1)%dispersion => mualem
	    pde(1)%pde_fnc(1)%convection => dmualem_dh
	    pde(1)%pde_fnc(1)%elasticity => vangen_elast
	    pde(1)%mass => vangen
	  else
	    call tabvalues(Kfnc=mualem, dKdhfnc = dmualem_dh, Cfnc=vangen_elast, thetafnc=vangen)
	    pde(1)%pde_fnc(1)%dispersion  => mualem_tab		
	    pde(1)%pde_fnc(1)%convection => dmualem_dh_tab
	    pde(1)%pde_fnc(1)%elasticity => vangen_elast_tab
	    pde(1)%mass => vangen_tab
	  end if
	  pde(1)%pde_fnc(1)%reaction => dummy_scalar
	  pde(1)%pde_fnc(1)%der_convect => dummy_vector
          pde(1)%pde_fnc(1)%zerord => dummy_scalar
          

	  do i=1, ubound(vgmatrix,1)
	    if (vgmatrix(i)%rcza_set%use) then
	      call init_zones(vgmatrix)
              pde(1)%dt_check => rcza_check
	      EXIT
	    else
	      pde(1)%dt_check => time_check_ok
	    end if
	  end do

	  ! adjustments for each specific type
	  select case(drutes_config%name)
	    case("RE_std", "RE_rot")
	  
	      do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
		select case(pde(1)%bc(i)%code)
		  case(-1)
		      pde(1)%bc(i)%value_fnc => re_dirichlet_height_bc
		  case(0)
			pde(1)%bc(i)%value_fnc => re_null_bc
		  case(1)
			pde(1)%bc(i)%value_fnc => re_dirichlet_bc
		  case(2)
			pde(1)%bc(i)%value_fnc => re_neumann_bc
		  case(3)
			pde(1)%bc(i)%value_fnc => re_null_bc
		  case default
			print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
			print *, "the incorrect boundary code specified is:", pde(1)%bc(i)%code
			ERROR stop
		end select
	      end do
	      
	      if (drutes_config%name == "RE_rot") then
		pde(1)%pde_fnc(1)%convection => convection_rerot
	      end if
	      pde(1)%flux => darcy_law
	      pde(1)%initcond => re_initcond

	      
	    case("RErotH", "REstdH")
	      do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
		select case(pde(1)%bc(i)%code)
		  case(-1)
		      pde(1)%bc(i)%value_fnc => retot_dirichlet_height_bc
		  case(0)
			pde(1)%bc(i)%value_fnc => re_null_bc
		  case(1)
			pde(1)%bc(i)%value_fnc => retot_dirichlet_bc
		  case(2)
			pde(1)%bc(i)%value_fnc => retot_neumann_bc
		  case(3)
			pde(1)%bc(i)%value_fnc => retot_freedrainage
		  case default
			print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
			print *, "the incorrect boundary code specified is:", pde(1)%bc(i)%code
			ERROR stop
		end select
	      end do
	      
	      if (drutes_config%name == "REstdH") then
		pde(1)%pde_fnc(1)%convection => dummy_vector
	      else if (drutes_config%name == "RErotH") then
		pde(1)%pde_fnc(1)%convection => convection_rerot
	      end if
	      
	      pde(1)%flux => darcy4totH
	      pde(1)%initcond => retot_initcond
! 	      pde(1)%initcond => iconddebouss
	      pde(1)%getval => getval_retot
	      
	    end select
	      
	      
	case("boussi")   

	  pde(1)%pde_fnc(1)%dispersion => bouss_cond
	  pde(1)%pde_fnc(1)%convection => bouss_adv
! pde(1)%pde_fnc(1)%convection =>dummy_vector
	  pde(1)%pde_fnc(1)%elasticity => bouss_elast
	  pde(1)%mass => dummy_scalar

	  pde(1)%pde_fnc(1)%reaction => dummy_scalar
	  pde(1)%pde_fnc(1)%der_convect => dummy_vector
          pde(1)%pde_fnc(1)%zerord => boussreact
	      
	  do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
	    pde(1)%bc(i)%value_fnc => bouss_bc
	  end do    
	   
	  pde(1)%flux => darcy4bouss
	  pde(1)%initcond => boussicond    
 
	  pde(1)%dt_check => time_check_ok

	  
	case("RE_mod")
	  !pde(1) Richards eq.
	  !pde(2) Heat flux
	  !pde(3) Transport of solutes
	  
	  pde_common%nonlinear = .true.
	  pde(1)%mass => mass_R
	  pde(1)%pde_fnc(1)%elasticity => elas_R_to_potential
	  pde(1)%pde_fnc(1)%dispersion => disp_R_to_potential
	  pde(1)%pde_fnc(1)%convection => conv_R_to_gravity
	  pde(1)%pde_fnc(1)%der_convect => dummy_vector
	  pde(1)%pde_fnc(1)%reaction => dummy_scalar
	  
	  pde(1)%pde_fnc(2)%elasticity => dummy_scalar  
	  pde(1)%pde_fnc(2)%dispersion => disp_R_to_heat
! 	  pde(1)%pde_fnc(2)%dispersion => dummy_tensor
	  pde(1)%pde_fnc(2)%convection => dummy_vector
	  pde(1)%pde_fnc(2)%der_convect => dummy_vector
	  pde(1)%pde_fnc(2)%reaction => dummy_scalar
	  
	  pde(1)%pde_fnc(3)%elasticity => dummy_scalar  
	  pde(1)%pde_fnc(3)%dispersion => disp_R_to_solute
! 	  pde(1)%pde_fnc(3)%dispersion => dummy_tensor
	  pde(1)%pde_fnc(3)%convection => dummy_vector
	  pde(1)%pde_fnc(3)%der_convect => dummy_vector
	  pde(1)%pde_fnc(3)%reaction => dummy_scalar
	  
	  
	  pde(2)%mass => dummy_scalar  
	  pde(2)%pde_fnc(1)%elasticity => dummy_scalar
	  pde(2)%pde_fnc(1)%dispersion => disp_H_to_potential
! 	  pde(2)%pde_fnc(1)%dispersion => dummy_tensor
	  pde(2)%pde_fnc(1)%convection => dummy_vector
	  pde(2)%pde_fnc(1)%der_convect => dummy_vector
	  pde(2)%pde_fnc(1)%reaction => dummy_scalar
	  
	  pde(2)%pde_fnc(2)%elasticity => elas_H_to_heat 
	  pde(2)%pde_fnc(2)%dispersion => disp_H_to_heat
	  pde(2)%pde_fnc(2)%convection => conv_H_to_heat
	  pde(2)%pde_fnc(2)%der_convect => dummy_vector
	  pde(2)%pde_fnc(2)%reaction => dummy_scalar
	  
	  pde(2)%pde_fnc(3)%elasticity => dummy_scalar  
	  pde(2)%pde_fnc(3)%dispersion => disp_H_to_solute
! 	  pde(2)%pde_fnc(3)%dispersion => dummy_tensor
	  pde(2)%pde_fnc(3)%convection => dummy_vector
	  pde(2)%pde_fnc(3)%der_convect => dummy_vector
	  pde(2)%pde_fnc(3)%reaction => dummy_scalar
	  
	  
	  pde(3)%mass => dummy_scalar  
	  pde(3)%pde_fnc(1)%elasticity => dummy_scalar
	  pde(3)%pde_fnc(1)%dispersion => dummy_tensor
	  pde(3)%pde_fnc(1)%convection => dummy_vector
	  pde(3)%pde_fnc(1)%der_convect => dummy_vector
	  pde(3)%pde_fnc(1)%reaction => dummy_scalar
	  
	  pde(3)%pde_fnc(2)%elasticity => dummy_scalar  
	  pde(3)%pde_fnc(2)%dispersion => dummy_tensor
	  pde(3)%pde_fnc(2)%convection => dummy_vector
	  pde(3)%pde_fnc(2)%der_convect => dummy_vector
	  pde(3)%pde_fnc(2)%reaction => dummy_scalar
	  
	  pde(3)%pde_fnc(3)%elasticity => elas_S_to_solute
	  pde(3)%pde_fnc(3)%dispersion => disp_S_to_solute
	  pde(3)%pde_fnc(3)%convection => conv_S_to_solute
	  pde(3)%pde_fnc(3)%der_convect => dummy_vector
	  pde(3)%pde_fnc(3)%reaction => dummy_scalar
	  
! 	  ! ! ! ! Pocatecni podminky ! ! ! ! 
	  pde(1)%initcond => modre_initcond
	  pde(2)%initcond => modheat_initcond
	  pde(3)%initcond => modsolute_initcond
	  ! ! ! ! 
	  

	  do i=1, ubound(pde,1)  
	    pde(i)%dt_check => time_check_ok
	  end do
	  
	  ! okrajove podminky transportu vody
	  do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
	    select case(pde(1)%bc(i)%code)
	      case(-1)
		  pde(1)%bc(i)%value_fnc => re_dirichlet_height_bc
	      case(0)
		    pde(1)%bc(i)%value_fnc => re_null_bc
	      case(1)
		    pde(1)%bc(i)%value_fnc => re_dirichlet_bc
	      case(2)
		    pde(1)%bc(i)%value_fnc => re_neumann_bc
	      case(3)
		    pde(1)%bc(i)%value_fnc => re_null_bc
	      case default
		    print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		    print *, "the incorrect boundary code specified is:", pde(1)%bc(i)%code
		    ERROR stop
	    end select
	  end do

	  ! okrajove podminky transportu tepla
	  do i=lbound(pde(2)%bc,1), ubound(pde(2)%bc,1)
	    select case(pde(2)%bc(i)%code)
	      case(0)
		    pde(2)%bc(i)%value_fnc => re_null_bc
	      case(1)
		    pde(2)%bc(i)%value_fnc => T_dirichlet_bc
	      case(2)
		    pde(2)%bc(i)%value_fnc => T_neumann_bc
	      case(3)
		    pde(2)%bc(i)%value_fnc => re_null_bc
	      case default
		    print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		    print *, "the incorrect boundary code specified is:", pde(1)%bc(i)%code
		    ERROR stop
	    end select
	  end do
	
	  ! okrajove podminky transportu kontaminantu
	  do i=lbound(pde(3)%bc,1), ubound(pde(3)%bc,1)
	    select case(pde(1)%bc(i)%code)
	      case(0)
		    pde(3)%bc(i)%value_fnc => re_null_bc
	      case(1)
		    pde(3)%bc(i)%value_fnc => C_dirichlet_bc
	      case(2)
		    pde(3)%bc(i)%value_fnc => C_neumann_bc
	      case(3)
		    pde(3)%bc(i)%value_fnc => re_null_bc
	      case default
		    print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		    print *, "the incorrect boundary code specified is:", pde(1)%bc(i)%code
		    ERROR stop
	    end select
	  end do
	  	    
		
	
	case("REtest")
          !pde(1) Richards eq.
	  !pde(2) Richards eq.
	  !pde(3) Richards eq.
	  pde_common%nonlinear = .true.
	  do i=1, ubound(pde,1)
	    pde(i)%pde_fnc(i)%dispersion => mualem
	    pde(i)%pde_fnc(i)%convection => dmualem_dh
	    pde(i)%pde_fnc(i)%elasticity => vangen_elast
	    pde(i)%mass => vangen
	    pde(i)%pde_fnc(i)%reaction => dummy_scalar
	    pde(i)%pde_fnc(i)%der_convect => dummy_vector
	    pde(i)%flux => darcy_law
	    pde(i)%initcond => re_initcond
	  end do
	  
	  pde(1)%pde_fnc(2)%elasticity => dummy_scalar  
	  pde(1)%pde_fnc(2)%dispersion => dummy_tensor
	  pde(1)%pde_fnc(2)%convection => dummy_vector
	  pde(1)%pde_fnc(2)%der_convect => dummy_vector
	  pde(1)%pde_fnc(2)%reaction => dummy_scalar
	  
	  pde(1)%pde_fnc(3)%elasticity => dummy_scalar  
	  pde(1)%pde_fnc(3)%dispersion => dummy_tensor
	  pde(1)%pde_fnc(3)%convection => dummy_vector
	  pde(1)%pde_fnc(3)%der_convect => dummy_vector
	  pde(1)%pde_fnc(3)%reaction => dummy_scalar
 	  
	  pde(2)%mass => dummy_scalar  
	  pde(2)%pde_fnc(1)%elasticity => dummy_scalar
	  pde(2)%pde_fnc(1)%dispersion => dummy_tensor
	  pde(2)%pde_fnc(1)%convection => dummy_vector
	  pde(2)%pde_fnc(1)%der_convect => dummy_vector
	  pde(2)%pde_fnc(1)%reaction => dummy_scalar
	  
	  pde(2)%pde_fnc(3)%elasticity => dummy_scalar  
	  pde(2)%pde_fnc(3)%dispersion => dummy_tensor 
	  pde(2)%pde_fnc(3)%convection => dummy_vector
	  pde(2)%pde_fnc(3)%der_convect => dummy_vector
	  pde(2)%pde_fnc(3)%reaction => dummy_scalar
	  
	  
	  pde(3)%mass => dummy_scalar  
	  pde(3)%pde_fnc(1)%elasticity => dummy_scalar
	  pde(3)%pde_fnc(1)%dispersion => dummy_tensor
	  pde(3)%pde_fnc(1)%convection => dummy_vector
	  pde(3)%pde_fnc(1)%der_convect => dummy_vector
	  pde(3)%pde_fnc(1)%reaction => dummy_scalar
	  
	  pde(3)%pde_fnc(2)%elasticity => dummy_scalar  
	  pde(3)%pde_fnc(2)%dispersion => dummy_tensor
	  pde(3)%pde_fnc(2)%convection => dummy_vector
	  pde(3)%pde_fnc(2)%der_convect => dummy_vector
	  pde(3)%pde_fnc(2)%reaction => dummy_scalar


	  do i=1, ubound(pde,1)  
	    pde(i)%dt_check => time_check_ok
	  end do
	    
	

	
	  do j=1, ubound(pde,1)
	  ! okrajove podminky transportu vody
	    do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
	      select case(pde(1)%bc(i)%code)
		case(-1)
		    pde(j)%bc(i)%value_fnc => re_dirichlet_height_bc
		case(0)
		      pde(j)%bc(i)%value_fnc => re_null_bc
		case(1)
		      pde(j)%bc(i)%value_fnc => re_dirichlet_bc
		case(2)
		      pde(j)%bc(i)%value_fnc => re_neumann_bc
		case(3)
		      pde(j)%bc(i)%value_fnc => re_null_bc
		case(5)
		      pde(j)%bc(i)%value_fnc => retot_atmospheric
		case default
		      print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		      print *, "the incorrect boundary code specified is:", pde(1)%bc(i)%code
		      ERROR stop
	      end select
	    end do	
	  end do
	  
	  vgset =>vgmatrix

      end select

      
      select case(drutes_config%dimen)
	case(1)
! 		    solve_matrix => CG_face
 	    solve_matrix => LDU_face
! 	    solve_matrix => CG_normal_face
	case(2)
!  solve_matrix => pcg
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