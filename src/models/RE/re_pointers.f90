
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

!> \file re_pointers.f90
!! \brief Pointer linkers for the Richards equation
!<

module RE_pointers
  public :: RE_std, REstdH
  public :: RE_pressheadbc,  RE_totheadbc, allREpointers
  public :: RE_processes
  
  contains
  
  
    subroutine RE_processes(processes)
      use typy
      
      integer(kind=ikind), intent(out) :: processes
      
      processes = 1
      
    end subroutine RE_processes

    
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
      

      if (drutes_config%rotsym) then
        pde_loc%pde_fnc(pde_loc%order)%convection => convection_rerot
        if (drutes_config%fnc_method == 0) then
          pde_loc%pde_fnc(pde_loc%order)%hidden_convect => dmualem_dh
        else
          pde_loc%pde_fnc(pde_loc%order)%hidden_convect => dmualem_dh_tab
        end if
      else
        if (drutes_config%fnc_method == 0) then
          pde_loc%pde_fnc(pde_loc%order)%convection => dmualem_dh
        else
          pde_loc%pde_fnc(pde_loc%order)%convection => dmualem_dh_tab
        end if
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
      
      
      if (drutes_config%rotsym) then
        pde_loc%pde_fnc(pde_loc%order)%convection => convection_rerot
      end if

      pde_loc%getval => getval_retot
      
      pde_loc%flux => darcy4totH
      pde_loc%initcond => retot_initcond
      pde_loc%symmetric = .true.
    
    end subroutine REstdH
    
    
    subroutine RE_pressheadbc(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use re_globals
      use re_constitutive
      use core_tools

      
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
          case(4)
            print *, "seepage face is not created yet for pressure head RE"
            ERROR STOP
          case(5)
            if (cut(drutes_config%name) == "REevap") then
              CONTINUE
              ! to be linked in evappointers
            else
              print *, "if you want to use atmospheric bc, set problem type to REevap"
              ERROR stop
            end if
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
      use core_tools
      use re_evap_methods
      
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
          case(4)
            pde_loc%bc(i)%value_fnc => retot_seepage	
          case(5)
            if (cut(drutes_config%name) == "REevap") then
              CONTINUE
              ! to be linked in evappointers
            else
              pde_loc%bc(i)%value_fnc => retot_atmospheric
            end if
          case(6)
            pde_loc%bc(i)%value_fnc => retot_well
          case(7, 8, 9)
              CONTINUE
              ! Linked in freeze_pointers
          case(10)
            pde_loc%bc(i)%value_fnc => retot_atmosphericv2
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
      

      pde_common%nonlinear = .true.
      if (drutes_config%fnc_method == 0) then
        pde_loc%pde_fnc(pde_loc%order)%dispersion => mualem
        pde_loc%pde_fnc(pde_loc%order)%elasticity => vangen_elast
        pde_loc%mass(1)%val => vangen
      else
        call tabvalues(pde_loc, Kfnc=mualem, dKdhfnc = dmualem_dh, Cfnc=vangen_elast, thetafnc=vangen)
        pde_loc%pde_fnc(pde_loc%order)%dispersion  => mualem_tab
        pde_loc%pde_fnc(pde_loc%order)%elasticity => vangen_elast_tab
        pde_loc%mass(1)%val => vangen_tab
      end if
      
      pde_loc%mass(2)%val => satdegree
      

      do i=1, ubound(vgset,1)
        if (vgset(i)%rcza_set%use) then
          call init_zones(vgset)
          pde_loc%dt_check => rcza_check
          EXIT
        end if
      end do
      
      pde_loc%pde_fnc(pde_loc%order)%zerord  => sinkterm
      
    
    end subroutine allREpointers
  

end module RE_pointers
