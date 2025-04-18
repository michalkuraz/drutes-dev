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

!> \file debug_tools.f90
!! \brief Debug tools for debugging DRUtES in terminal. 
!<

!> Contains generic procedure printmtx for printing different types of arrays, and 
!! substitution for pause (depricated procedure)
!<


module debug_tools
  private :: print_real_matrix, print_int_matrix, print_real_vector, print_int_vector, print_real_vector4
  public :: wait, print_filename
  private :: print_quadpnt
  
  !> generic procedure, can be overloaded with different vector/matrix types
  interface printmtx
    module procedure print_real_matrix
    module procedure print_int_matrix
    module procedure print_real_vector
    module procedure print_char_vector
    module procedure print_int_vector
    module procedure print_sparse_matrix
    module procedure print_smartmatrix_i
    module procedure print_smartarray_i
    module procedure print_smartmatrix_r
    module procedure print_smartarray_r
    module procedure print_logical_array
    module procedure print_real_vector4
    module procedure print_quadpnt
  end interface printmtx
 

  contains
  
    function print_filename(unt) result(answer)
      use typy

      integer, intent(in) :: unt
      integer :: ids
      character(len=512) :: filename
      character(:), allocatable :: answer
      logical :: openedq
      
      ids = unt
      
      inquire(unit=ids, opened=openedq, name=filename)
      
      if (openedq) then
        answer = trim(filename)
      else
        answer = "file not opened"
      end if
      
    end function print_filename

    !> prints integpnt_str
    subroutine print_quadpnt(quadpnt, filunit, name)
      use global_objs
      use core_tools
      use globals
      
      type(integpnt_str), intent(in) :: quadpnt
      integer, intent(in), optional :: filunit
      character(len=*), intent(in), optional :: name
      
      integer :: filloc
      integer :: ierr
      logical :: op   
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      
      write(unit=filloc, fmt=*) "type:", quadpnt%type_pnt
      
      if (quadpnt%type_pnt == "ndpt" .or. quadpnt%type_pnt == "obpt") then
        write(unit=filloc, fmt=*) "order of point:", quadpnt%order
      end if
      
      if (quadpnt%type_pnt == "gqnd") then
        write(unit=filloc, fmt=*) "order of element:", quadpnt%element
      end if
      
      if (quadpnt%type_pnt == "xypt") then
        write(unit=filloc, fmt=*) "coordinates:", quadpnt%xy
      end if
      
      write(unit=filloc, fmt=*) "column:", quadpnt%column
      
      write(unit=filloc, fmt=*) "element:", quadpnt%element
      
      if (quadpnt%ddlocal) then
        write(unit=filloc, fmt=*) "using subdomain local data"
        write(unit=filloc, fmt=*) "subdomain id", quadpnt%subdom
        if (quadpnt%extended) then
          write(unit=filloc, fmt=*) "node from extended subdomain (see subcycling man)"
        else
        write(unit=filloc, fmt=*) "node inside the subdomain"
        end if
      else
        write(unit=filloc, fmt=*) "using global data"
      end if
      
      if (.not. quadpnt%globtime) then
        write(unit=filloc, fmt=*) "not using the global time, solution will be interpolated"
        write(unit=filloc, fmt=*) "time 4 use:", quadpnt%time4eval
      else
        write(unit=filloc, fmt=*) "using the global time"
      end if
      
  
      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
      
    end subroutine print_quadpnt
  
    !> prints sparse matrix  
    subroutine print_sparse_matrix(A, filunit, name)
      use sparsematrix
      use mtxiotools
      use globals
      use core_tools
      class(smtx), intent(in out) :: A
      integer, intent(in), optional :: filunit
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
    
      if (present(name)) then
        print *, "cdswav"
        call a%write(name)
        return
      end if
      

      call a%print()

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if

    end subroutine print_sparse_matrix

    !> prints dense matrix of reals
    subroutine print_real_matrix(a, filunit, name)
      use typy
      use globals
      use core_tools

      real(kind=rkind), dimension(:,:), intent(in out) :: a
      integer, intent(in), optional :: filunit     
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      
      do i=lbound(a,1), ubound(a,1)
        write(unit=filloc, fmt=*)  i, "|",  a(i,:)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
      
    end subroutine print_real_matrix

    !> prints dense matrix of integers
    subroutine print_int_matrix(a, filunit, name)
      use typy
      use globals
      use core_tools
      
      integer(kind=ikind), dimension(:,:), intent(in) :: a
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      
      do i=lbound(a,1), ubound(a,1)
        write(unit=filloc, fmt=*)  i, "|",  a(i,:)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
      
    end subroutine print_int_matrix
    
        !> prints vector of characters
    subroutine print_char_vector(V, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      character(len=*), dimension(:), intent(in) :: V  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
     

      do i=lbound(V,1),ubound(V,1)
       write(unit=filloc, fmt=*)  i,  V(i)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
  
    end subroutine print_char_vector


    !> prints vector of reals
    subroutine print_real_vector(V, filunit, name)
    ! vytiskne vektor, pocet sloupcu tisku je nc
      use typy
      use globals
      use core_tools
      
      !parametry
      real(kind=rkind), dimension(:), intent(in) :: V  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
     

      do i=lbound(V,1),ubound(V,1)
       write(unit=filloc, fmt=*)  i,  V(i)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
  
    end subroutine print_real_vector


        !> prints vector of single reals
    subroutine print_real_vector4(V, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      real(4), dimension(:), intent(in) :: V  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
     

      do i=lbound(V,1),ubound(V,1)
       write(unit=filloc, fmt=*) "row:", i, "value:", V(i)
      end do

      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
  
    end subroutine print_real_vector4

    !> prints vector of integers
    subroutine print_int_vector(V, filunit, name)
    ! vytiskne vektor, pocet sloupcu tisku je nc
      use typy
      use globals
      use core_tools
      
      !parametry
      integer(kind=ikind), dimension(:), intent(in) :: V  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if

      do i=lbound(V,1),ubound(V,1)
         write(unit=filloc, fmt=*) "row:", i, "value:", V(i)
      end do
  
      if (terminal /= filloc) then
        close(filloc)
      else
        call flush(terminal)
      end if
  

    end subroutine print_int_vector
    
    
    !> prints smartarray vector of smartarray vectors of integers
    subroutine print_smartmatrix_i(array, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      class(smartarray_int), dimension(:), intent(in) :: array  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
            
      do i=lbound(array,1),ubound(array,1)
        if (.not. allocated(array(i)%data)) then
          print *, "no values to print at row:", i
        else
          write(unit=filloc, fmt=*) i, "|", array(i)%data(1:array(i)%pos)
          write(unit=filloc, fmt=*) "-----------------------------------------------"
        end if
      end do
      
      call flush(filloc)
      
    end subroutine print_smartmatrix_i

    
    !> prints smartarray vector of integers
    subroutine print_smartarray_i(array, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      class(smartarray_int), intent(in) :: array  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      do i=1, array%pos
        write(unit=filloc, fmt=*)  i, "|", array%data(i)
      end do

      
      call flush(filloc)
      
    end subroutine print_smartarray_i
    
    
    
        !> prints smartarray vector of smartarray vectors of integers
    subroutine print_smartmatrix_r(array, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      class(smartarray_real), dimension(:), intent(in) :: array  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
            
      do i=lbound(array,1),ubound(array,1)
        if (.not. allocated(array(i)%data)) then
          print *, "no values to print at row:", i
        else
          write(unit=filloc, fmt=*) i, "|", array(i)%data(1:array(i)%pos)
          write(unit=filloc, fmt=*) "-----------------------------------------------"
        end if
      end do
      
      call flush(filloc)
      
    end subroutine print_smartmatrix_r

    
    !> prints smartarray vector of integers
    subroutine print_smartarray_r(array, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      class(smartarray_real), intent(in) :: array  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      do i=1, array%pos
        write(unit=filloc, fmt=*)  i, "|", array%data(i)
      end do

      
      call flush(filloc)
      
    end subroutine print_smartarray_r
    
    
    
    !> prints vector of logicals
    subroutine print_logical_array(array, filunit, name)
      use typy
      use globals
      use core_tools
      
      !parametry
      logical, dimension(:), intent(in) :: array  !<vektor k tisknuti
      integer, intent(in), optional :: filunit   
      character(len=*), intent(in), optional :: name

      integer :: filloc
      integer :: ierr
      logical :: op
      integer(kind=ikind) :: i
      
      if (present(name)) then
        call find_unit(filloc)
        open(unit=filloc, file=name, action="write", status="replace", iostat=ierr)
        if (ierr /= 0) then
          print *, "unable to open dump file, called from debug_tools::printmtx"
          error stop
        end if
      else if (present(filunit)) then
        filloc = filunit
        inquire(unit=filloc, opened=op)
        if (.not. op) then
          print *, "file not opened, called from debug_tools::printmtx"
          error stop
        end if
      else
        filloc = terminal
      end if
      
      do i=1, ubound(array,1)
        write(unit=filloc, fmt=*) i, "|", array(i)
      end do
      
    end subroutine print_logical_array




    !> substitution of depricated function pause
    subroutine wait(ch) 
      use globals       
      character(len=*), optional, intent(in) :: ch
      
      if (present(ch)) then
        print *, "-------------"
        print *, trim(ch)
        print *, "-------------"
      end if

      print *, "press [ENTER] to continue"
  
      call flush(terminal)
  
      read(*,*)


    end subroutine wait
    
    
!   


end module debug_tools
