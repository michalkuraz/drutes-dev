module Re_dual_pointers
  public :: RE_fracture
  public :: RE_matrix
  public :: dual_neumann_bc
  public :: dual_freedrainage
  public :: dual_atmospheric
  public :: inf_neumann_bc
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
      pde(1)%initcond => dual_inicond

      pde(1)%flux => darcy_law_d
      if (drutes_config%fnc_method == 0) then
	    pde(1)%pde_fnc(1)%dispersion => dual_mualem
	    select case(coup_model)
	     case(1:3)
	       pde(1)%pde_fnc(1)%reaction => dual_coupling_neg
	       pde(1)%pde_fnc(2)%reaction => dual_coupling
	     case(4:5)
	       pde(1)%pde_fnc(1)%reaction => dual_coup_min_neg
	       pde(1)%pde_fnc(2)%reaction => dual_coup_min
	     case default
	     stop
	    end select
	    
	    pde(1)%pde_fnc(1)%elasticity => dual_ret_cap
	    pde(1)%mass => vangen_d
      else
	    call dual_tabvalues(pde(1), Kfnc=dual_mualem, Cfnc=dual_ret_cap,&
            thetafnc=vangen_d,ex_K_fnc=dual_coupling_K)
	    pde(1)%pde_fnc(1)%dispersion  => dual_mualem_tab		
	    pde(1)%pde_fnc(1)%reaction => dual_coupling_neg_tab
	    pde(1)%pde_fnc(2)%reaction => dual_coupling_tab
	    pde(1)%pde_fnc(1)%elasticity => dual_ret_cap_tab
	    pde(1)%mass => vangen_d_tab
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
            pde(1)%bc(i)%value_fnc => dual_neumann_bc
	  case(3)
            pde(1)%bc(i)%value_fnc => dual_freedrainage
	  case(4)
            pde(1)%bc(i)%value_fnc => dual_atmospheric
	  case(5)
            pde(1)%bc(i)%value_fnc => inf_neumann_bc
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
      pde(2)%initcond => dual_inicond

      
     if (drutes_config%fnc_method == 0) then
	    pde(2)%pde_fnc(2)%dispersion => dual_mualem
	    select case(coup_model)
	     case(1:3)
	       pde(2)%pde_fnc(2)%reaction => dual_coupling_neg
	       pde(2)%pde_fnc(1)%reaction => dual_coupling
	     case(4:5)
	       pde(2)%pde_fnc(2)%reaction =>dual_coup_min_neg
	       pde(2)%pde_fnc(1)%reaction =>dual_coup_min
	     case default
	     stop
	    end select
	    pde(2)%pde_fnc(2)%elasticity => dual_ret_cap
	    pde(2)%mass => vangen_d
      else
      	    call dual_tabvalues(pde(2), Kfnc=dual_mualem, Cfnc=dual_ret_cap,&
            thetafnc=vangen_d,ex_K_fnc=dual_coupling_K)
	    pde(2)%pde_fnc(2)%dispersion  => dual_mualem_tab		
	    pde(2)%pde_fnc(2)%reaction => dual_coupling_neg_tab
	    pde(2)%pde_fnc(1)%reaction => dual_coupling_tab
	    pde(2)%pde_fnc(2)%elasticity => dual_ret_cap_tab
	    pde(2)%mass => vangen_d_tab
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
            pde(2)%bc(i)%value_fnc => dual_neumann_bc
	  case(3)
            pde(2)%bc(i)%value_fnc => dual_freedrainage
	  case(4)
            pde(2)%bc(i)%value_fnc => dual_atmospheric
	  case(5)
            pde(2)%bc(i)%value_fnc => inf_neumann_bc
	  case default
		print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		print *, "the incorrect boundary code specified is:", pde(2)%bc(i)%code
		ERROR stop
	end select
      end do 
     

    end subroutine RE_fracture

 
  subroutine dual_neumann_bc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals

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
	       
              select case (pde_loc%mfswitch)
                case("m")
                  bcval = pde_loc%bc(edge_id)%series(j,2)*exchange(layer)%weightm
                case("f")
                  bcval = pde_loc%bc(edge_id)%series(j,2)*exchange(layer)%weightf
              end select
	      EXIT
	    end if
	  end do
	else
	  layer=pde_loc%bc(edge_id)%layer
          select case (pde_loc%mfswitch)
            case("m")
              bcval = pde_loc%bc(edge_id)%value*exchange(layer)%weightm
            case("f")
              bcval = pde_loc%bc(edge_id)%value*exchange(layer)%weightf
          end select
	end if
	


	value = bcval
	print*, pde_loc%mfswitch, value

      end if
      
      if (present(code)) then
	code = 2
      end if


    end subroutine dual_neumann_bc
    
  subroutine inf_neumann_bc(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
     

      integer(kind=ikind) :: i, edge_id, j,layer
      real(kind=rkind), dimension(3) :: gravflux, bcflux
      real(kind=rkind) :: bcval, gfluxval, weight
      integer :: i1
      
      select case (pde_loc%mfswitch)
        case("m")
          weight=(1_rkind-infweight)
        case("f")
          weight=infweight
      end select

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
              bcval = pde_loc%bc(edge_id)%series(j,2)*weight
	      EXIT
	    end if
	  end do
	else
	  layer=pde_loc%bc(edge_id)%layer
          bcval = pde_loc%bc(edge_id)%value*weight
	end if

	value = bcval

      end if
      
      if (present(code)) then
	code = 2
      end if


    end subroutine inf_neumann_bc
    
 
  subroutine dual_freedrainage(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals
      use dual_por
      
      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      real(kind=rkind), dimension(3,3) :: K
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer, D
      real(kind=rkind), dimension(3) :: gravflux
      
      
      if (present(value)) then

	quadpnt%type_pnt = "ndpt"
	quadpnt%column = 2
	quadpnt%order = elements%data(el_id, node_order)
	layer = elements%material(el_id,1)
	D = drutes_config%dimen
        select case (pde_loc%mfswitch)
          case("m")
            call pde(1)%pde_fnc(1)%dispersion(pde_loc, layer, quadpnt, tensor=K(1:D,1:D))          
          case("f")
            call pde(2)%pde_fnc(2)%dispersion(pde_loc, layer, quadpnt, tensor=K(1:D,1:D))         
        end select
	
	
      	select case(D)
	  case(1)
	  
	    value = K(1,1) * elements%nvect_z(el_id, node_order)
	  
	  case(2)	  
	    gravflux(1) = sqrt(1-elements%nvect_z(el_id, node_order)*elements%nvect_z(el_id, node_order))*K(1,2)
	    
	    gravflux(2) = elements%nvect_z(el_id, node_order)*K(2,2)

	    value = sqrt(gravflux(1)*gravflux(1) + gravflux(2)*gravflux(2))

	end select
      end if
      
      if (present(code)) then
	code = 2
      end if
      
    end subroutine dual_freedrainage
     
  subroutine dual_atmospheric(pde_loc, el_id, node_order, value, code) 
      use typy
      use globals
      use global_objs
      use pde_objs
      use dual_globals

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in)  :: el_id, node_order
      real(kind=rkind), intent(out), optional    :: value
      integer(kind=ikind), intent(out), optional :: code
      
      
      
      type(integpnt_str) :: quadpnt
      integer(kind=ikind) :: layer
      real(kind=rkind) :: theta, rain, evap
      integer(kind=ikind) :: i, edge_id, j
      
      

      if (present(code)) then
	    code = 2
      end if
      
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
	      rain = pde_loc%bc(edge_id)%series(j,2)
	      evap = pde_loc%bc(edge_id)%series(j,3)
	      EXIT
	    end if
	  end do
	else
	  print *, "atmospheric boundary must be time dependent, check record for the boundary", edge_id
	  ERROR STOP
	end if

        quadpnt%type_pnt = "ndpt"
        quadpnt%column = 2
        quadpnt%order = elements%data(el_id,node_order)
        layer = elements%material(el_id,1)
        theta =  pde_loc%mass(layer, quadpnt)
        select case (pde_loc%mfswitch)
          case("m")
            value = (rain - evap*theta**(2.0_rkind/3.0_rkind))*(1_rkind-infweight)  
          case("f")
            value = (rain - evap*theta**(2.0_rkind/3.0_rkind))*(infweight)
        end select


      end if
      
    end subroutine dual_atmospheric
     
end module Re_dual_pointers