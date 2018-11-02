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

!> \file bousspointers.f90
!! \brief Boussinesq equation pointer linker.
!<





module bousspointers
  public :: boussi, bouss_processes
  
  contains
  
   subroutine bouss_processes(processes)
      use typy
      
      integer(kind=ikind), intent(out) :: processes
      
      processes = 1
      
    end subroutine bouss_processes
  
    subroutine boussi(pde_loc)
      use typy
      use globals
      use global_objs
      use pde_objs
      use boussglob
      use boussfnc
      use boussread
	
      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind) :: i
      
      call boussreader(pde_loc)
      
      pde_loc%pde_fnc(pde_loc%order)%dispersion => bouss_cond
      pde_loc%pde_fnc(pde_loc%order)%convection => bouss_adv
      pde_loc%pde_fnc(pde_loc%order)%elasticity => bouss_elast

      pde_loc%pde_fnc(pde_loc%order)%zerord => boussreact
	  
      do i=lbound(pde(1)%bc,1), ubound(pde(1)%bc,1)
        select case(pde_loc%bc(i)%code)
          case(1,2)
            pde_loc%bc(i)%value_fnc => bouss_bc
        end select
      end do
	      
	
      pde_loc%flux => darcy4bouss
      pde_loc%initcond => boussicond     
       
    end subroutine boussi


end module bousspointers
