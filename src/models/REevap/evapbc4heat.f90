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


module evapbc4heat
  private :: RNterm, Hterm, Eterm
  public :: ebalance_flux
  
  contains
    function RNterm(quadpnt) result(Rn)
      use typy
      use global_objs
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: Rn
    
    end function RNterm
    
    
    function Hterm(quadpnt) result(H)
      use typy
      use global_objs
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: H
    end function Hterm
    
    function Eterm(quadpnt) result(E)
      use typy
      use global_objs
      type(integpnt_str), intent(in) :: quadpnt
      real(kind=rkind) :: E
    
    end function Eterm
    
    subroutine ebalance_flux()
    
    end subroutine ebalance_flux


end module evapbc4heat
