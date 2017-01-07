module RE_pointers
  public :: RE_std, RE_rot, RErotH, REstdH
  public :: RE_pressheadbc,  RE_totheadbc, allREpointers
  
  contains
  
  
     subroutine RErotH(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use re_total
      
      class(pde_str), intent(in out) :: pde_loc
      
      call allREpointers(pde_loc)
      call RE_totheadbc(pde_loc)
      
      pde_loc%pde_fnc(pde_loc%order)%convection => convection_rerot
      pde_loc%getval => getval_retot
      pde_loc%flux => darcy4totH
      pde_loc%initcond => retot_initcond

    end subroutine RErotH
    
    subroutine RE_rot(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      
      class(pde_str), intent(in out) :: pde_loc
      
      call allREpointers(pde_loc)
      call RE_pressheadbc(pde_loc)
      
      pde_loc%pde_fnc(pde_loc%order)%convection => convection_rerot   

      
      if (drutes_config%fnc_method == 0) then
      	pde_loc%pde_fnc(pde_loc%order)%convection => dmualem_dh
      else
      	pde_loc%pde_fnc(pde_loc%order)%convection => dmualem_dh_tab
      end if
      pde_loc%flux => darcy_law
      pde_loc%initcond => re_initcond
     
    
    end subroutine RE_rot
    
    subroutine RE_std(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      
      class(pde_str), intent(in out) :: pde_loc
      
      call allREpointers(pde_loc)
      call RE_pressheadbc(pde_loc)
       
      pde_loc%getval => getvalp1
      
      if (drutes_config%fnc_method == 0) then
      	pde_loc%pde_fnc(pde_loc%order)%convection => dmualem_dh
      else
      	pde_loc%pde_fnc(pde_loc%order)%convection => dmualem_dh_tab
      end if
	
      pde_loc%flux => darcy_law
      pde_loc%initcond => re_initcond
     
    
    end subroutine RE_std
    
    subroutine REstdH(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use re_total
      
      class(pde_str), intent(in out) :: pde_loc
      
      call allREpointers(pde_loc)
      call RE_totheadbc(pde_loc)

      pde_loc%getval => getval_retot
      
      pde_loc%flux => darcy4totH
      pde_loc%initcond => retot_initcond
    
    end subroutine REstdH
    
    
    subroutine RE_pressheadbc(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      
      class(pde_str), intent(in out) :: pde_loc
      
      integer(kind=ikind) :: i
      
      
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
	select case(pde_loc%bc(i)%code)
	  case(-1)
	      pde_loc%bc(i)%value_fnc => re_dirichlet_height_bc
	  case(0)
		pde_loc%bc(i)%value_fnc => re_null_bc
	  case(1)
		pde_loc%bc(i)%value_fnc => re_dirichlet_bc
	  case(2)
		pde_loc%bc(i)%value_fnc => re_neumann_bc
	  case(3)
		pde_loc%bc(i)%value_fnc => re_null_bc
	  case default
		print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		print *, "the incorrect boundary code specified is:", pde_loc%bc(i)%code
		ERROR stop
	end select
      end do
	

    end subroutine RE_pressheadbc
    
    subroutine RE_totheadbc(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use re_total
      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i
      
      do i=lbound(pde_loc%bc,1), ubound(pde_loc%bc,1)
	select case(pde_loc%bc(i)%code)
	  case(-1)
	      pde_loc%bc(i)%value_fnc => retot_dirichlet_height_bc
	  case(0)
		pde_loc%bc(i)%value_fnc => re_null_bc
	  case(1)
		pde_loc%bc(i)%value_fnc => retot_dirichlet_bc
	  case(2)
		pde_loc%bc(i)%value_fnc => retot_neumann_bc
	  case(3)
		pde_loc%bc(i)%value_fnc => retot_freedrainage
	  case default
		print *, "ERROR! You have specified an unsupported boundary type definition for the Richards equation"
		print *, "the incorrect boundary code specified is:", pde_loc%bc(i)%code
		ERROR stop
	end select
      end do
    
    
    end subroutine RE_totheadbc
 
    
    
    subroutine allREpointers(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive   
      use re_reader
      
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i	
      logical, save :: read=.false.
      
      ! read inputs
      if (.not. read) then
	call res_read(pde_loc)
	read = .true.
      end if
      
      
      call domainswitch("m")
      pde_common%nonlinear = .true.
      if (drutes_config%fnc_method == 0) then
	pde_loc%pde_fnc(pde_loc%order)%dispersion => mualem
	pde_loc%pde_fnc(pde_loc%order)%elasticity => vangen_elast
	pde_loc%mass => vangen
      else
	call tabvalues(pde_loc, Kfnc=mualem, dKdhfnc = dmualem_dh, Cfnc=vangen_elast, thetafnc=vangen)
	pde_loc%pde_fnc(pde_loc%order)%dispersion  => mualem_tab		
	pde_loc%pde_fnc(pde_loc%order)%elasticity => vangen_elast_tab
	pde_loc%mass => vangen_tab
      end if
      

      do i=1, ubound(vgmatrix,1)
	if (vgmatrix(i)%rcza_set%use) then
	  call init_zones(vgmatrix)
	  pde_loc%dt_check => rcza_check
	  EXIT
	end if
      end do
      
      pde_loc%pde_fnc(pde_loc%order)%zerord  => sinkterm
      

      
      
    
    end subroutine allREpointers
  

end module RE_pointers