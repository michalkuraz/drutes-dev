module Re_dual_pointers
  public :: RE_fracture
  public :: RE_matrix
  
  contains
! change pointers to take tabular version
! tabularize coupling function with coupling pars
    subroutine RE_matrix(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_por
      use Re_dual_reader
      use RE_constitutive
      use debug_tools
      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i

           
      pde_loc%getval => getval_retot_dual
      call Re_dual_readm(pde_loc)
      call Re_dual_var() 
      

      
      pde_loc%flux => darcy_law_d!implement darcy law for fluxes
      if (drutes_config%fnc_method == 0) then
	    pde_loc%pde_fnc(pde_loc%order)%dispersion => dual_mualemm
	    pde_loc%pde_fnc(pde_loc%order)%reaction => dual_coupling
	    pde_loc%pde_fnc(pde_loc%order)%elasticity => dual_ret_capm
	    pde_loc%mass => vangen_d_m
      else
	    call dual_tabvalues(pde_loc, Kfnc=dual_mualemm, Cfnc=dual_ret_capm,&
	     thetafnc=vangen_d_m,Kfnc_f=dual_mualemf, Cfnc_f=dual_ret_capf, &
	     thetafnc_f=vangen_d_f,ex_K_fnc=dual_coupling_K)
	    pde_loc%pde_fnc(pde_loc%order)%dispersion  => dual_mualem_m_tab		
	    pde_loc%pde_fnc(pde_loc%order)%reaction => dual_coupling_tab
	    pde_loc%pde_fnc(pde_loc%order)%elasticity => dual_ret_capm_tab
	    pde_loc%mass => vangen_d_m_tab
      end if
      
      ! boundary condition defined as different type boundary_vals
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
	   select case(pde_loc%bc(i)%code)
	   	 case(0)
	       pde_loc%bc(i)%value_fnc => re_null_bc
	     case(1)
	       pde_loc%bc(i)%value_fnc => re_dirichlet_bc
	     case(2)
	       pde_loc%bc(i)%value_fnc => re_neumann_bc
	   end select
      end do  
    pde_loc%initcond => dual_inicond_m
   
    end subroutine RE_matrix
    
    subroutine RE_fracture(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_por
      use Re_dual_reader
      use RE_constitutive
      use debug_tools
      
      class(pde_str), intent(in out) :: pde_loc  
      integer(kind=ikind) :: i
      
      pde_loc%getval => getval_retot_dual
      
      call Re_dual_readf(pde_loc)
      
     if (drutes_config%fnc_method == 0) then
	    pde_loc%pde_fnc(pde_loc%order)%dispersion => dual_mualemf
	    pde_loc%pde_fnc(pde_loc%order)%reaction => dual_coupling_f
	    pde_loc%pde_fnc(pde_loc%order)%elasticity => dual_ret_capf
	    pde_loc%mass => vangen_d_m
      else
	    pde_loc%pde_fnc(pde_loc%order)%dispersion  => dual_mualem_f_tab		
	    pde_loc%pde_fnc(pde_loc%order)%reaction => dual_coupling_f_tab
	    pde_loc%pde_fnc(pde_loc%order)%elasticity => dual_ret_capf_tab
	    pde_loc%mass => vangen_d_f_tab
      end if
      
      pde_loc%flux => darcy_law_d
      
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
	   select case(pde_loc%bc(i)%code)
	   	 case(0)
	       pde_loc%bc(i)%value_fnc => re_null_bc
	     case(1)
	       pde_loc%bc(i)%value_fnc => re_dirichlet_bc
	     case(2)
	       pde_loc%bc(i)%value_fnc => re_neumann_bc
	   end select
     end do  
     
    pde_loc%initcond => dual_inicond_f
    end subroutine RE_fracture
      

end module Re_dual_pointers