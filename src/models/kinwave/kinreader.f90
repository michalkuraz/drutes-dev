
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



module kinreader

  private :: read_catchment, gen_catchment
  public :: kininit

  contains
  
  
    subroutine read_catchment()
      use typy
      use kinglobs
      use readtools
      
      integer :: ierr, filenodes, fileel
      integer(kind=ikind) :: ndcounter, itmp, i
      
      open(newunit=filenodes, file="drutes.conf/kinwave/nodes.in", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "unable to open drutes.conf/kinwave/nodes.in, exiting..."
        error stop
      end if
      
      open(newunit=fileel, file="drutes.conf/kinwave/elements.in", status="old", action="read", iostat=ierr)
    
      if (ierr /= 0) then
        print *, "unable to open drutes.conf/kinwave/elements.in, exiting..."
        error stop
      end if
      
      ndcounter = 0
      do 
        call comment(filenodes)
        read(unit=filenodes, fmt=*, iostat=ierr) itmp
        if (ierr /= 0) then
          EXIT 
        end if
        if (itmp > ndcounter) then
          ndcounter = itmp
        end if
      end do
      
      allocate(watershed_nd(ndcounter))
      
      close(filenodes)
      
      open(newunit=filenodes, file="drutes.conf/kinwave/nodes.in", status="old", action="read", iostat=ierr)
      
      print *, ndcounter
      
      do i=1, ndcounter
        call fileread(itmp, filenodes) 
        print *, itmp
      end do
      stop
    
    end subroutine read_catchment
    
    
    subroutine gen_catchment()
      use typy
      use kinglobs
      use globals
      use geom_tools
      
      integer(kind=ikind) :: i
      real(kind=rkind), dimension(3) :: a,b,c
      
      allocate(watershed_nd(nodes%kolik))
      
      allocate(watershed_el(elements%kolik))
      
      
      do i=1, nodes%kolik
        watershed_nd(i)%xyz(1:2) = nodes%data(i,:)
        watershed_nd(i)%xyz(3) = nodes%data(i,1)**2*0.001 + nodes%data(i,2)**2*0.002
      end do
      
      watershed_el(:)%cover = 1
      
      do i=1, elements%kolik
        a = watershed_nd(elements%data(i,1))%xyz
        b = watershed_nd(elements%data(i,2))%xyz
        c = watershed_nd(elements%data(i,3))%xyz
        
        call plane_derivative(a,b,c, watershed_el(i)%sx, watershed_el(i)%sy)
        
      end do
      

    
    end subroutine gen_catchment
    
    

    
    subroutine kininit(pde_loc)
      use typy
      use pde_objs
      use kinglobs
      use readtools
      use core_tools
      
      class(pde_str), intent(in out) :: pde_loc
      integer :: file_kinematix, ierr, filerain
      integer(kind=ikind) :: n, i, counter
      character(len=512) :: msg
      real(kind=rkind) :: tmp
      
      select case(drutes_config%mesh_type)
        case(1) 
          call gen_catchment()
          
        case(4)
          call read_catchment()
          
        case default
          print *, "for kinematic wave equation use Arc GIS data option only"
          print *, "  this is option number 4"
          print *, "exiting..."
          ERROR STOP
          
      end select
      
      pde_loc%problem_name(1) = "runoff"
      pde_loc%problem_name(2) = "Kinematic wave equation for real catchments"

      pde_loc%solution_name(1) = "runoff_depth" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "[L] " !popisek grafu

      pde_loc%flux_name(1) = "runoff_flux"  
      pde_loc%flux_name(2) = "runoff flux [L3.T-1.L-2]"
    
      allocate(pde_loc%mass_name(1,2))
      
      pde_loc%mass_name(1,1) = "runoff_height"
      pde_loc%mass_name(1,2) = "runoff elevation [m a.m.s.l]"
      
      open(newunit=file_kinematix, file="drutes.conf/kinwave/kinwave.conf", action="read", status="old", iostat = ierr)
      
      n = maxval(elements%material)
      
      write(msg, fmt=*) "You have specified incorrect number of subregions with different Manning cofficients,", new_line("a"),  &
      "according to mesh definitions, you should define", n, "different Manning values." 
      
      
      allocate(manning(n))
      
      call fileread(i, file_kinematix, ranges=(/n,n/), errmsg=cut(msg))
      
      do i=1, n
        call fileread(manning(i), file_kinematix, ranges=(/epsilon(manning(1)), huge(manning(1))/))
      end do
      
      open(newunit=filerain, file="drutes.conf/kinwave/rain.in", status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "unable to open file with rainfall data drutes.conf/kinwave/rain.in"
        error stop
      end if
      
      
      counter = 0
      do
        call comment(filerain)
        read(unit=filerain, fmt=*, iostat=ierr) tmp
        
        if (ierr /= 0) then
          EXIT
        else
          counter = counter + 1
        end if
      end do
      
      close(filerain)
      
      open(newunit=filerain, file="drutes.conf/kinwave/rain.in", status="old", action="read", iostat=ierr)
      
      allocate( 
    
    end subroutine kininit

end module kinreader
