
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

module kinpointer

  contains
  
    subroutine kinwaveprocs(number) 
      use typy
      integer(kind=ikind), intent(out) :: number
      
      number = 1
      
    end subroutine kinwaveprocs
    
    subroutine kinwavelinker(pde_loc)
      use typy
      use kinfnc
      use kinreader
      use pde_objs
      
      class(pde_str), intent(in out) :: pde_loc
      
      call kininit(pde_loc)
      
      pde_loc%pde_fnc(pde_loc%order)%convection => kinconvect
      
      pde_loc%bc(101)%value_fnc => kinbor
      
      pde_loc%initcond => kinematixinit
      
      
    end subroutine kinwavelinker
    
    



end module kinpointer
