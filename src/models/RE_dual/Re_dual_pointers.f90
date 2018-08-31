
! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher

! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.

!> \file Re_dual_pointers.f90
!! \brief This module contains the pointers used for dual permeability model
!<

module Re_dual_pointers
  public :: RE_fracture
  public :: RE_matrix
  contains
!> Links to releveant subroutines for matrix domain 
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
      use re_dual_bc

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
    
!> Links to releveant subroutines for fracture domain 
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
      use re_dual_bc

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

     
end module Re_dual_pointers
