module Re_dual_pointers
  public :: RE_fracture
  public :: RE_matrix
  
  contains

    subroutine RE_matrix()
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_por
      use Re_dual_reader
      use RE_constitutive
      use debug_tools
      use dual_tab
      use dual_coup
      use re_total
      !class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i

           
      pde(1)%getval => getval_retot_dual
      call Re_dual_readm(pde(1))
      call Re_dual_var() 
      pde(1)%initcond => dual_inicond_m

      pde(1)%flux => darcy_law_d
      if (drutes_config%fnc_method == 0) then
	    pde(1)%pde_fnc(1)%dispersion => dual_mualemm
	    select case(coup_model)
	     case(1:2)
	       pde(1)%pde_fnc(1)%reaction => dual_coupling
	       pde(1)%pde_fnc(2)%reaction => dual_coupling_f
	     case(3:4)
	       pde(1)%pde_fnc(1)%reaction =>dual_coup_min
	       pde(1)%pde_fnc(2)%reaction =>dual_coup_min_f
	     case default
	     stop
	    end select
	    
	    pde(1)%pde_fnc(1)%elasticity => dual_ret_capm
	    pde(1)%mass => vangen_d_m
      else
	    call dual_tabvalues(pde(1), Kfnc=dual_mualemm, Cfnc=dual_ret_capm,&
	     thetafnc=vangen_d_m,Kfnc_f=dual_mualemf, Cfnc_f=dual_ret_capf, &
	     thetafnc_f=vangen_d_f,ex_K_fnc=dual_coupling_K)
	    pde(1)%pde_fnc(1)%dispersion  => dual_mualem_m_tab		
	    pde(1)%pde_fnc(1)%reaction => dual_coupling_tab
	    pde(1)%pde_fnc(2)%reaction => dual_coupling_f_tab
	    pde(1)%pde_fnc(1)%elasticity => dual_ret_capm_tab
	    pde(1)%mass => vangen_d_m_tab
      end if
      
      ! boundary condition defined as different type boundary_vals
      do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
	select case(pde(1)%bc(i)%code)
	  case(-1)
	      pde(1)%bc(i)%value_fnc => retot_dirichlet_height_bc
	  case(0)
		pde(1)%bc(i)%value_fnc => re_null_bc
	  case(1)
		pde(1)%bc(i)%value_fnc => retot_dirichlet_bc
	  case(2)
		pde(1)%bc(i)%value_fnc => retot_neumann_bc_m
	  case default
		print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		print *, "the incorrect boundary code specified is:", pde(1)%bc(i)%code
		ERROR stop
	end select
      end do 

   
    end subroutine RE_matrix
    
    subroutine RE_fracture()
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_por
      use Re_dual_reader
      use RE_constitutive
      use debug_tools
      use dual_tab
      use dual_coup
      use re_total
      !class(pde_str), intent(in out) :: pde_loc  
      integer(kind=ikind) :: i
      
      pde(2)%getval => getval_retot_dual
      call Re_dual_readf(pde(2))
      pde(2)%initcond => dual_inicond_f

      
     if (drutes_config%fnc_method == 0) then
	    pde(2)%pde_fnc(2)%dispersion => dual_mualemf
	    select case(coup_model)
	     case(1:2)
	       pde(2)%pde_fnc(2)%reaction => dual_coupling_f
	       pde(2)%pde_fnc(1)%reaction => dual_coupling
	     case(3:4)
	       pde(2)%pde_fnc(2)%reaction =>dual_coup_min_f
	       pde(2)%pde_fnc(1)%reaction =>dual_coup_min
	     case default
	     stop
	    end select
	    pde(2)%pde_fnc(2)%elasticity => dual_ret_capf
	    pde(2)%mass => vangen_d_f
      else
	    pde(2)%pde_fnc(2)%dispersion  => dual_mualem_f_tab		
	    pde(2)%pde_fnc(2)%reaction => dual_coupling_f_tab
	    pde(2)%pde_fnc(1)%reaction => dual_coupling_tab
	    pde(2)%pde_fnc(2)%elasticity => dual_ret_capf_tab
	    pde(2)%mass => vangen_d_f_tab
      end if
      
      pde(2)%flux => darcy_law_d
      
      do i=lbound(pde(2)%bc,1), ubound(pde(2)%bc,1)
	select case(pde(2)%bc(i)%code)
	  case(-1)
	      pde(2)%bc(i)%value_fnc => retot_dirichlet_height_bc
	  case(0)
		pde(2)%bc(i)%value_fnc => re_null_bc
	  case(1)
		pde(2)%bc(i)%value_fnc => retot_dirichlet_bc
	  case(2)
		pde(2)%bc(i)%value_fnc => retot_neumann_bc_f
	  case default
		print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		print *, "the incorrect boundary code specified is:", pde(2)%bc(i)%code
		ERROR stop
	end select
      end do 
     

    end subroutine RE_fracture

 
 subroutine retot_neumann_bc_m(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals
      !use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
     

      integer(kind=ikind) :: i, edge_id, j,layer
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, weight
      integer :: i1
      
      

      if (present(value)) then
	edge_id = nodes%edge(elements%data(el_id, node_order))

	i = pde_loc%permut(elements%data(el_id, node_order))
	

	if (pde_loc%bc(edge_id)%file) then
	  do i=1, ubound(pde_loc%bc(edge_id)%series,1)
	    if (pde_loc%bc(edge_id)%series(i,1) > time) then
	      if (i > 1) then
		j = i-1
	      else
		j = i
	      end if
	      layer=pde_loc%bc(edge_id)%layer
	      bcval = pde_loc%bc(edge_id)%series(j,2)*exchange(layer)%weightm
	      EXIT
	    end if
	  end do
	else
	  layer=pde_loc%bc(edge_id)%layer
	  bcval = pde_loc%bc(edge_id)%value*exchange(layer)%weightm
	end if
	


	value = bcval

      end if
      
      if (present(code)) then
	code = 2
      end if


    end subroutine retot_neumann_bc_m
    
     subroutine retot_neumann_bc_f(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals
      !use re_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
     

      integer(kind=ikind) :: i, edge_id, j, layer
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, weight
      integer :: i1
      
      

      if (present(value)) then
	edge_id = nodes%edge(elements%data(el_id, node_order))

	i = pde_loc%permut(elements%data(el_id, node_order))
	
	

	if (pde_loc%bc(edge_id)%file) then
	  do i=1, ubound(pde_loc%bc(edge_id)%series,1)
	    if (pde_loc%bc(edge_id)%series(i,1) > time) then
	      if (i > 1) then
		j = i-1
	      else
		j = i
	      end if
	      layer=pde_loc%bc(edge_id)%layer
	      bcval = pde_loc%bc(edge_id)%series(j,2)*exchange(layer)%weightf
	      EXIT
	    end if
	  end do
	else
	 layer=pde_loc%bc(edge_id)%layer
	  bcval = pde_loc%bc(edge_id)%value*exchange(layer)%weightf
	end if
	


	value = bcval

      end if
      
      if (present(code)) then
	code = 2
      end if


    end subroutine retot_neumann_bc_f
 
end module Re_dual_pointers