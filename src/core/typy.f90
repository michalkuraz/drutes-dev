
! Copyright 2008 Michal Kuraz, Petr Mayer


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

module typy

    
    integer, parameter, public :: lowikind=selected_int_kind(5)
    
    integer, parameter, public :: dprec=selected_real_kind(15,99) 
    
    integer, parameter, public :: sprec=selected_real_kind(8,9)

    !> kind pro realna cisla 18 cifer
    integer, parameter, public :: rkind = selected_real_kind(15,99)
    !> kind pro cela cisla 10 cifer
    integer, parameter, public :: ikind = selected_int_kind(10)
    !> velmi dlouha cela cisla
    integer, parameter, public :: likind = selected_int_kind(16)

    !> pocitadlo operaci
    type, public :: tcount
        !> pocet aditivnich operaci
        integer(kind=likind) :: ad = 0
        !> pocet nasobeni
        integer(kind=likind) :: mul = 0
        !> pocet deleni
        integer(kind=likind) :: div = 0
        !> doba behu
        real(kind=rkind) :: time = 0
    end type tcount

    public :: print_info
    public :: update_info
    contains
    !> vytiskne udaje o pocitani
    subroutine print_info(info)
        implicit none
        !> data o spotrebe prace
        type(tcount), intent(in) :: info
        print *, "pocty operaci aditivni:",info%ad," nasobeni:",info%mul,&
            " deleni:",info%div, " cas:",info%time
    end subroutine print_info

    !> pricte info2 k info1
    subroutine update_info(info1,info2)
        implicit none
        type(tcount), intent(inout) :: info1
        type(tcount), intent(in) :: info2

        info1%ad = info1%ad + info2%ad
        info1%mul = info1%mul + info2%mul
        info1%div = info1%div + info2%div
        info1%time = info1%time + info2%time

    end subroutine update_info
    
!     function norm2(a) result(cislo)
!       real(kind=rkind), dimension(:), intent(in) :: a
!       real(kind=rkind) :: cislo
!       
!       integer(kind=ikind) :: i
!       
!       cislo = 0
!       do i=1, ubound(a,1)
! 	cislo = cislo + a(1)*a(1)
!       end do
!       
!       cislo = sqrt(cislo)
!     
!     end function norm2


end module typy
