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


module core_tools
  public :: write_log

  !> dummy functions, the pointer procedure is set to this function in case if not being considered
  public :: avg
  public :: mesh_allocater
  public :: find_unit
  public :: cut


  contains

  function cut(ch) result(out_ch)
  
    character(len=*), intent(in) :: ch
    character(len=:), allocatable :: out_ch
    
    allocate(character(len=LEN(adjustl(TRIM(ch)))) :: out_ch)
    
    out_ch=adjustl(TRIM(ch))
    
  end function cut
  
  subroutine write_log(text, real1, int1, text2, real2, int2, text3, real3, int3, hidden)
    use typy
    use globals
    character(len=*), intent(in) :: text
    real(kind=rkind), intent(in), optional :: real1, real2, real3
    integer(kind=ikind), intent(in), optional :: int1, int2, int3
    character(len=*), intent(in), optional :: text2, text3
    !> character to store the time info
    character(len=10) :: timer
    !> the character to store the date info
    character(len=8) :: dater

    character(len=1024) :: ch1, ch2, ch3, ch4, ch5, ch6, ch7, ch8
    logical, optional :: hidden

    if (present(real1)) then
      write(unit=ch1, fmt=*) real1
    else
      write(unit=ch1, fmt=*) ""
    end if


    if (present(int1)) then
      write(unit=ch2, fmt=*) int1
    else
      write(unit=ch2, fmt=*) ""
    end if

    if (present(text2)) then
      write(unit=ch3, fmt=*) text2
    else
      write(unit=ch3, fmt=*) ""
    end if

    if (present(real2)) then
      write(unit=ch4, fmt=*) real2
    else
      write(unit=ch4, fmt=*) ""
    end if

    if (present(int2)) then
      write(unit=ch5, fmt=*) int2
    else
      write(unit=ch5, fmt=*) ""
    end if

    
    if (present(text3)) then
      write(unit=ch6, fmt=*) text3
    else
      write(unit=ch6, fmt=*) ""
    end if

    if (present(real3)) then
      write(unit=ch7, fmt=*) real3
    else
      write(unit=ch7, fmt=*) ""
    end if

    if (present(int3)) then
      write(unit=ch8, fmt=*) int3
    else
      write(unit=ch8, fmt=*) ""
    end if

    call date_and_time(dater, timer)
    write(unit=logfile, fmt = *) "---------------------------------------------------------------"
    write(unit=logfile, fmt = *)  timer(1:2), "hrs", " : ", timer(3:4), "min"," : ", timer(5:10), "s", "  -  ", dater(1:4), &
       "/", dater(5:6), "/", dater(7:8), "    ",  trim(text), trim(ch1), trim(ch2), trim(ch3), trim(ch4), trim(ch5), &
	    trim(ch6), trim(ch7), trim(ch8)
    if (.not. present(hidden)) then	    
      write(unit=terminal, fmt=*) trim(text), trim(ch1), trim(ch2), trim(ch3), trim(ch4), trim(ch5), &
		      trim(ch6), trim(ch7), trim(ch8)
    else 
      if (.not. hidden) then
	write(unit=terminal, fmt=*) trim(text), trim(ch1), trim(ch2), trim(ch3), trim(ch4), trim(ch5), &
	      trim(ch6), trim(ch7), trim(ch8)
      end if
    end if
      

 
!     if (present(time)) then
!       write(unit=logfile, fmt = *) "simulation time=", time
!     end if
    call flush(logfile)

  end subroutine write_log



    function avg(a,b) result(c)
      use typy
      real(kind=rkind), intent(in) :: a
      real(kind=rkind), intent(in) :: b
      real(kind=rkind) :: c

      c = (a+b)/2.0_rkind

    end function avg

    !> allocates structures for nodes and elements
    subroutine mesh_allocater()
      use typy
      use global_objs
      use pde_objs
      use globals

      integer :: n

      allocate(nodes%id(nodes%kolik))
      allocate(nodes%edge(nodes%kolik))
      allocate(nodes%data(nodes%kolik,drutes_config%dimen))
      allocate(nodes%boundary(nodes%kolik))
      allocate(nodes%el2integ(nodes%kolik))
      allocate(nodes%element(nodes%kolik))
      allocate(elements%length(elements%kolik,drutes_config%dimen + 1))
      allocate(elements%nvect_z(elements%kolik,drutes_config%dimen + 1))
      allocate(elements%id(elements%kolik))
      allocate(elements%data(elements%kolik, drutes_config%dimen+1))
      allocate(elements%areas(elements%kolik))
      allocate(elements%ders(elements%kolik, drutes_config%dimen+1, drutes_config%dimen))
      
      n = pde_common%processes
      allocate(elements%material(elements%kolik, n))
      
      allocate(elements%neighbours(elements%kolik, ubound(elements%data,2)))

      

    end subroutine mesh_allocater


    subroutine find_unit(iunit, start_id)
      integer, intent(out) :: iunit
      integer, intent(in), optional :: start_id
      logical :: op

      if (present(start_id)) then
        iunit = start_id*1 + 1
      else
        iunit = 100 ! the lower values are often reserved
      end if

      do
        inquire(unit=iunit,opened=op)
        if (.not.op) then
          exit
        end if 
        iunit = iunit + 1  
        if (iunit >= huge(iunit) - 1) then
          ERROR STOP "SYSTEM PANIC, unable to find any free unit"
          exit
        end if
      end do
  
  end subroutine find_unit





end module core_tools
