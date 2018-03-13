
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
    
    

    !> real number specification
    integer, parameter, public :: rkind = selected_real_kind(15,99)
    !> integer number specification
    integer, parameter, public :: ikind = selected_int_kind(10)
    !> long integers
    integer, parameter, public :: likind = selected_int_kind(16)
    !> unicode char type
    integer, parameter, public :: chkind=selected_char_kind("default")!"ASCII")!"ISO_10646")


    !---------petr mayer's definitions----------
    integer, parameter, public :: maxlen = 200
    !> seznam retezcu
    type, public :: StringList
        !> jednotlive retezce
        character(len=maxlen), dimension(:), allocatable :: name
        !> delky retezcu
        integer, dimension(:), allocatable :: delka
        !> skutecny pocet retezcu
        integer :: pocet = 0
    end type StringList

    !> 4byte real
    integer, parameter :: r4  = selected_real_kind(5,10)
    !> 8byte real
    integer, parameter :: r8  = selected_real_kind(12,20)
    !> 16byte real
    integer, parameter :: r16 = selected_real_kind(30)

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

    !> maximalni delka retezce jedne polozky beznych dat
    integer, parameter, public :: strmax = 1024

    !> typ pro retezec s promennou delkou
    type, public :: Tident
        !> skutecny obsah
        character(len=:), allocatable :: nm
        !prozatim jen
        ! dopsat metody
        !character(len=strmax), allocatable :: nm ! jen pro stary linux
    end type Tident


    !> koren pro funkce
    type, public :: funroot
        contains
        !> vycisleni funkce v komplexnim oboru
        procedure evalc
        !>  vycisleni funkce v realnem oboru
        procedure evalr
        !> pretizena verze vycisleni funkce
        generic :: eval => evalr, evalc
        !> vycisli derivaci funkce, komplexni
        procedure derevalc
        !> vycisli derivaci funkce, realna
        procedure derevalr
        !> vycisli derivaci funkce, pretizena
        generic :: dereval => derevalc, derevalr
    end type funroot

    public :: print_info
    public :: update_info
    public :: norm2 ! casem zrusit, duvodem je chyba v 16b implementaci
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
        type(tcount), intent(in)    :: info2

        info1%ad   = info1%ad   + info2%ad
        info1%mul  = info1%mul  + info2%mul
        info1%div  = info1%div  + info2%div
        info1%time = info1%time + info2%time

    end subroutine update_info


    function evalr(r,x) result(y)
        implicit none
        class(funroot), intent(in) :: r
        real(kind=rkind), intent(in) :: x
        real(kind=rkind) :: y
        y = 0
        stop "Neimplentovana funkcionalita"
    end function evalr

    function evalc(r,x) result(y)
        implicit none
        class(funroot), intent(in) :: r
        complex(kind=rkind), intent(in) :: x
        complex(kind=rkind) :: y
        y = 0
        stop "Neimplentovana funkcionalita"
    end function evalc

    function derevalr(r,x) result(y)
        implicit none
        class(funroot), intent(in) :: r
        real(kind=rkind), intent(in) :: x
        real(kind=rkind) :: y
        stop "Neimplentovana funkcionalita"
    end function derevalr

    function derevalc(r,x) result(y)
        implicit none
        class(funroot), intent(in) :: r
        complex(kind=rkind), intent(in) :: x
        complex(kind=rkind) :: y
        stop "Neimplentovana funkcionalita"
    end function derevalc


    !> spocte L2 normu vektoru
    function norm2(v) result(y)
        implicit none
        real(kind=rkind), dimension(:), intent(in) :: v
        real(kind=rkind) :: y
        integer :: i

        y=0
        do i = lbound(v,1), ubound(v,1)
            y = y + v(i)*v(i)
        end do
        y = sqrt(y)
    end function norm2

end module typy



